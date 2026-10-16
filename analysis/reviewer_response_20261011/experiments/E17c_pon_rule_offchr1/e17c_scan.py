"""E17c scan: PoN features of every labelled tensor and every somatic truth record of one chromosome.

    python e17c_scan.py chr22 [--workers 4] [--check 300]      -> DATA/scan/chr22.json

Nothing here decides a rule: the scan stores, per GRCh38 allele, what each PoN says about it, and
e17c_tables.py turns that into the rule grid. Inputs are read only.

Tensors. Every tensor of label 0 (non), 1 (somatic) or 2 (germline) in <SET>/tensors/<KIND>/<chrom>_labels.ndjson
of the six HPRC v1.1 d9 sets and their six HPRC v2.1 d46 sets (label -1 = ignore is skipped). The allele the PoN
sees is the one a call of that tensor would have in the GRCh38 VCF: the labels' `grch38` projection (the same
projection merge stored in <chrom>_variant_summary.ndjson), checked against the FASTA, INS / DEL left-aligned and
padded with pansoma_net_v2.linear_vcf.Genome (imported, as linear_vcf does it). No projection = `unplaced`: no PoN
can reach the call.

Truth. The somatic truth VCF of each sample (labels.manifest.json): PASS or '.' records, kind as
pansoma_net_v2.vcfeval.record_kinds (the truth vcfeval scores) and as filter_panel_of_normals.record_kind (the
rule the PoN applies to the record); `in_bed` = the record's REF span overlaps the somatic benchmark BED (vcfeval
-e). Records as they are, as vcfeval tags its truth (truth_pon_tags).

PoN features (scripts/filter_panel_of_normals.py, imported, not modified). Each PoN VCF is streamed once over the
chromosome; its records at the positions of interest are kept in memory and handed to the script's own
record_matches() through an object with the same fetch() interface (record_matches keeps only records starting at
POS, so this is the same comparison as its tabix fetch; --check repeats a random sample against the real tabix
files). Per allele (POS, REF, ALT):
  g_any, c_any     gnomAD / CoLoRSdb have the exact allele at any AF (record_matches, min_af=None)
  g_lvl, c_lvl     highest floor of FLOORS at which record_matches(min_af=floor) matches; -1 = none
  d_any, d_nonsom  dbSNP has the allele; ... in a record not flagged SAO=2 (record_matches non_somatic=True)
  k_any            the 1000G PoN has the allele
  hp               homopolymer_run() of the script (INDEL alleles; 0 otherwise), capped at 50
  rc               indexed_gam_pipeline_v4 truth_labels.RepeatContext class of the INDEL event (HP>=7, HP4-6, STR,
                   VNTR, none, complex); SNV / MNP for equal-length alleles
"""
from __future__ import annotations

import argparse
import gc
import json
import multiprocessing
import random
import sys
import time
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

REPO = Path("/scratch/jshen/Github/Pansoma")
sys.path[:0] = [str(REPO), str(REPO / "scripts"), str(REPO / "machine_learning")]

import pysam  # noqa: E402

from filter_panel_of_normals import (homopolymer_run, record_kind, record_matches,  # noqa: E402
                                     resolve_contig)
from indexed_gam_pipeline_v4.tensor_postprocessing.truth_labels import RepeatContext, split_alleles  # noqa: E402
from pansoma_net_v2.linear_vcf import Genome  # noqa: E402
from pansoma_net_v2.vcfeval import record_kinds  # noqa: E402

pysam.set_verbosity(0)
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E17c_pon_rule_offchr1")
PON_DIR = Path("/scratch/jshen/data/Pansoma/panel_of_normal_VCFs")
PON_FILES = ("af-only-gnomad.hg38.vcf.gz", "Homo_sapiens_assembly38.dbsnp138.vcf.gz", "1000g_pon.hg38.vcf.gz",
             "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")  # the script's PON_NAMES order
GNOMAD, DBSNP, KG, COLORS = range(4)
FLOORS = (1e-5, 1e-4, 1e-3, 1e-2, 5e-2)  # every AF floor of the rule grid
D9 = ("HG008T_PacBio", "HG008T_ONT", "HG008T_Illumina", "COLO829T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT")
SETS = D9 + tuple(f"{s}_d46" for s in D9)
SAMPLE_SET = {"HG008T": "HG008T_PacBio", "COLO829T": "COLO829T_ONT"}  # manifest that names the sample's truth
FEATURES = ("g_any", "g_lvl", "c_any", "c_lvl", "d_any", "d_nonsom", "k_any", "hp", "rc")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ tensors
def read_labels(task):
    """Counter (label, reason, (chrom, pos0, ref, alt, event) | None) of the labelled tensors of one file."""
    set_name, kind, chrom = task
    path = TENSORS / set_name / "tensors" / kind / f"{chrom}_labels.ndjson"
    out = Counter()
    if not path.exists():
        return task, None
    with open(path) as stream:
        for line in stream:
            if '"label": -1,' in line:
                continue
            r = json.loads(line)
            if r["label"] not in (0, 1, 2):
                continue
            g = r.get("grch38")
            lin = None
            if g:
                lin = (g["chrom"], g["pos0"], g["ref"].upper(), g["alt"].upper(), r["candidate_id"].split(":")[2])
            out[(r["label"], r["reason"], lin)] += 1
    return task, out


# ------------------------------------------------------------------ the chromosome in memory
class MemFasta:
    """One chromosome of the FASTA in memory with pysam.FastaFile.fetch()'s interface (end clipped at the contig
    end), handed to Genome, homopolymer_run and RepeatContext instead of random-access FASTA reads."""

    def __init__(self, path, chrom):
        with pysam.FastaFile(path) as fa:
            self.seq = fa.fetch(chrom)
        self.chrom = chrom

    def fetch(self, contig, start, end):
        if contig != self.chrom:
            raise KeyError(contig)
        return self.seq[max(0, start):end]


# ------------------------------------------------------------------ PoNs in memory for record_matches
class Rec:
    __slots__ = ("pos", "ref", "alts", "info")

    def __init__(self, pos, ref, alts, info):
        self.pos, self.ref, self.alts, self.info = pos, ref, alts, info


class MemPoN:
    """The records of one PoN that start at the positions of interest, with the fetch() record_matches uses."""

    def __init__(self, filename, records):
        self.filename = filename
        self.at = defaultdict(list)
        for pos, ref, alts, info in records:
            self.at[pos].append(Rec(pos, ref, alts, info))

    def fetch(self, contig, start, stop):
        # record_matches asks for (POS - 1, POS) and keeps candidates with .pos == POS
        return self.at.get(stop, ())


def stream_pon(task):
    """(pos, REF, ALTs, info) of every record of PoN p on the chromosome at a position of interest."""
    p, chrom, positions = task
    path = PON_DIR / PON_FILES[p]
    out = []
    n = 0
    with pysam.VariantFile(str(path)) as vcf:
        contig = resolve_contig(vcf, chrom)
        if contig is None:
            return p, contig, out, 0
        has_af, has_sao = "AF" in vcf.header.info, "SAO" in vcf.header.info
        for rec in vcf.fetch(contig):
            n += 1
            if rec.pos not in positions:
                continue
            info = {}
            if has_af:
                info["AF"] = rec.info.get("AF")
            if has_sao:
                info["SAO"] = rec.info.get("SAO")
            out.append((rec.pos, rec.ref, rec.alts, info))
    return p, contig, out, n


def level(pon, contig, pos, ref, alt):
    """Highest FLOORS index at which record_matches(min_af) matches; -1 = none."""
    best = -1
    for i, floor in enumerate(FLOORS):
        if record_matches(pon, contig, pos, ref, (alt,), min_af=floor):
            best = i
        else:
            break
    return best


def features(pons, contigs, key, fasta, repeats, in_memory=True):
    """Feature tuple of one allele; in_memory=False queries the PoN objects without the position shortcut."""
    chrom, pos, ref, alt = key
    f = dict.fromkeys(FEATURES, 0)
    f["g_lvl"] = f["c_lvl"] = -1
    if not in_memory or any(pos in pon.at for pon in pons):
        g, d, k, c = pons
        f["g_any"] = int(record_matches(g, contigs[GNOMAD], pos, ref, (alt,)))
        f["c_any"] = int(record_matches(c, contigs[COLORS], pos, ref, (alt,)))
        f["g_lvl"] = level(g, contigs[GNOMAD], pos, ref, alt) if f["g_any"] else -1
        f["c_lvl"] = level(c, contigs[COLORS], pos, ref, alt) if f["c_any"] else -1
        f["d_any"] = int(record_matches(d, contigs[DBSNP], pos, ref, (alt,)))
        f["d_nonsom"] = int(record_matches(d, contigs[DBSNP], pos, ref, (alt,), non_somatic=True))
        f["k_any"] = int(record_matches(k, contigs[KG], pos, ref, (alt,)))
    if len(ref) != len(alt):
        f["hp"] = min(50, homopolymer_run(fasta, chrom, pos, ref, alt))
        events = list(split_alleles(pos, ref, (alt,)))
        if len(events) == 1 and events[0][3] in ("INS", "DEL"):
            p0, r, a, _, _ = events[0]
            f["rc"] = repeats.classify(chrom, p0, r, a)
        else:
            f["rc"] = "complex"
    else:
        f["rc"] = "SNV/MNP"
    return tuple(f[k] for k in FEATURES)


# ------------------------------------------------------------------ truth
def bed_intervals(path, chrom):
    iv = []
    with open(path) as fh:
        for line in fh:
            c = line.split("\t")
            if c[0] == chrom:
                iv.append((int(c[1]), int(c[2])))
    iv.sort()
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return merged


def overlaps(merged, starts, a, b):
    import bisect
    i = bisect.bisect_right(starts, a) - 1
    for j in (i, i + 1):
        if 0 <= j < len(merged) and merged[j][0] < b and a < merged[j][1]:
            return True
    return False


def read_truth(chrom):
    """[(sample, pon_kind, vcfeval kinds, in_bed, keys)] of the PASS / '.' truth records on the chromosome."""
    out = []
    for sample, set_name in SAMPLE_SET.items():
        m = json.loads((TENSORS / set_name / "tensors" / "SNV" / "labels.manifest.json").read_text())
        vcf_path, bed_path = m["truth"]["somatic"]["vcf"], m["truth"]["somatic"]["bed"]
        merged = bed_intervals(bed_path, chrom)
        starts = [a for a, _ in merged]
        with pysam.VariantFile(vcf_path) as vcf:
            contig = resolve_contig(vcf, chrom)
            if contig is None:
                continue
            for v in vcf.fetch(contig):
                if set(v.filter.keys()) - {"PASS"}:
                    continue
                alts = tuple(a for a in (v.alts or ()) if a)
                kinds = record_kinds(v.ref, ",".join(alts))
                if not kinds:
                    continue
                in_bed = overlaps(merged, starts, v.pos - 1, v.pos - 1 + max(1, len(v.ref)))
                keys = tuple((chrom, v.pos, v.ref.upper(), a.upper()) for a in alts
                             if not a.startswith("<") and a not in ("*", "."))
                out.append((sample, record_kind(v.ref, alts), tuple(sorted(kinds)), in_bed, keys))
    return out


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("chrom")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--check", type=int, default=300, help="alleles re-checked against the tabix PoN files")
    ap.add_argument("--sets", nargs="*", default=list(SETS))
    args = ap.parse_args()
    chrom = args.chrom
    t0 = time.time()
    m = json.loads((TENSORS / D9[0] / "tensors" / "SNV" / "labels.manifest.json").read_text())
    fasta_path = m["fasta"]
    genome, repeats, fasta = Genome(fasta_path), RepeatContext(fasta_path), MemFasta(fasta_path, chrom)
    genome.fasta = repeats.genome = fasta  # same reads, from memory

    # 1. labelled tensors
    tasks = [(s, k, chrom) for s in args.sets for k in ("SNV", "INDEL")]
    with Pool(args.workers) as pool:
        label_counts = dict(pool.map(read_labels, tasks))
    log(chrom, "labels read", sum(sum(c.values()) for c in label_counts.values() if c), "tensors")

    norm = {}       # (chrom, pos0, ref, alt, ev) -> key | "refmismatch"

    def normalize(lin):
        if lin in norm:
            return norm[lin]
        c, pos0, ref, alt, ev = lin
        if c not in genome.lengths or genome.get(c, pos0, pos0 + len(ref)) != ref:
            norm[lin] = "refmismatch"
            return norm[lin]
        if ev in ("INS", "DEL"):
            pos0, ref, alt, _ = genome.left_align(c, pos0, ref, alt, ev)
        pos, vref, valt = genome.vcf_allele(c, pos0, ref, alt, ev)
        norm[lin] = (c, pos, vref, valt)
        return norm[lin]

    keys = {}

    def key_id(key):
        return keys.setdefault(key, len(keys))

    tensor_rows = Counter()   # (set, kind, label, reason, key id | 'unplaced' | 'refmismatch' | 'otherchrom')
    for (set_name, kind, _), counts in label_counts.items():
        if counts is None:
            continue
        for (label, reason, lin), n in counts.items():
            if lin is None:
                k = "unplaced"
            elif lin[0] != chrom:
                k = "otherchrom"
            else:
                key = normalize(lin)
                k = key if key == "refmismatch" else key_id(key)
            tensor_rows[(set_name, kind, label, reason, k)] += n
    missing = [list(t) for t, c in label_counts.items() if c is None]
    del label_counts
    norm.clear()
    gc.collect()
    truth = read_truth(chrom)
    truth_rows = Counter()
    for sample, pon_kind, kinds, in_bed, tkeys in truth:
        truth_rows[(sample, pon_kind, kinds, in_bed, tuple(key_id(k) for k in tkeys))] += 1
    log(chrom, f"{len(keys):,} distinct alleles; {len(truth):,} truth records")

    # 2. PoN records at the positions of interest
    positions = {k[1] for k in keys}
    with Pool(min(4, args.workers)) as pool:
        streamed = pool.map(stream_pon, [(p, chrom, positions) for p in range(4)])
    pons, contigs, n_streamed = [None] * 4, [None] * 4, {}
    for p, contig, recs, n in streamed:
        pons[p] = MemPoN(PON_FILES[p], recs)
        contigs[p] = contig
        n_streamed[PON_FILES[p]] = n
    log(chrom, "PoNs streamed", n_streamed)

    # 3. features per allele
    feat_ids, feats = {}, []
    allele_feat = [0] * len(keys)
    for key, i in keys.items():
        f = features(pons, contigs, key, fasta, repeats)
        allele_feat[i] = feat_ids.setdefault(f, len(feats))
        if allele_feat[i] == len(feats):
            feats.append(f)
    log(chrom, f"features done: {len(feats):,} distinct feature tuples")

    # 4. check: a random sample of alleles against the tabix PoN files, as the script queries them
    rng = random.Random(f"E17c:{chrom}")
    sample_keys = rng.sample(sorted(keys), min(args.check, len(keys)))
    # half of the sample among alleles some PoN has at POS, so the check covers matches too
    with_records = sorted(k for k in keys if any(k[1] in pon.at for pon in pons))
    sample_keys += rng.sample(with_records, min(args.check, len(with_records)))
    files = [pysam.VariantFile(str(PON_DIR / f)) for f in PON_FILES]
    mismatches = []
    for key in sample_keys:
        mine = feats[allele_feat[keys[key]]]
        theirs = features(files, contigs, key, fasta, repeats, in_memory=False)
        if mine != theirs:
            mismatches.append([list(key), list(mine), list(theirs)])
    for f in files:
        f.close()
    log(chrom, f"check: {len(sample_keys)} alleles vs tabix, {len(mismatches)} differ")

    out = dict(chrom=chrom, floors=FLOORS, features=FEATURES, sets=args.sets,
               feature_tuples=[list(f) for f in feats],
               tensors=[[s, k, lab, reason, (allele_feat[i] if isinstance(i, int) else i), n]
                        for (s, k, lab, reason, i), n in sorted(tensor_rows.items(), key=str)],
               truth=[[sample, pk, list(kinds), in_bed, [allele_feat[i] for i in ids], n]
                      for (sample, pk, kinds, in_bed, ids), n in sorted(truth_rows.items(), key=str)],
               pon_records_streamed=n_streamed, distinct_alleles=len(keys),
               missing_label_files=missing,
               check=dict(alleles=len(sample_keys), mismatches=mismatches[:50], n_mismatches=len(mismatches)),
               seconds=round(time.time() - t0, 1))
    (DATA / "scan").mkdir(parents=True, exist_ok=True)
    tmp = DATA / "scan" / f"{chrom}.json.tmp"
    tmp.write_text(json.dumps(out))
    tmp.rename(DATA / "scan" / f"{chrom}.json")
    log(chrom, f"done in {out['seconds']} s")


if __name__ == "__main__":
    multiprocessing.set_start_method("forkserver")  # workers do not inherit (and copy) the parent's tables
    main()

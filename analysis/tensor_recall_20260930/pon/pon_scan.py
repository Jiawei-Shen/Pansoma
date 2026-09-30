"""Somatic truth alleles and tensors of HG008T / COLO829T against the four PoN VCFs, one chromosome per run.

    python pon_scan.py chr1      -> scan/chr1.json

Rule = scripts/filter_panel_of_normals.py (the repo PoN used by pansoma_net_v2.vcfeval): a PoN record at the same
POS with the same REF and one ALT equal; gnomAD and CoLoRSdb only if that ALT's AF >= 1e-4, dbSNP only if not
SAO=2, 1000G any. The same comparison is done here by streaming each PoN over the chromosome once (the script's
tabix fetch keeps only records starting at POS, so the result is the same), and three looser overlaps are
recorded beside it:
  position  any PoN record starts at POS (any REF/ALT)
  allele    same POS, REF, ALT at any AF / SAO
  trimmed   the same after removing the common suffix of REF and ALT on both sides (a multi-allelic PoN record
            pads REF for its longest allele, e.g. REF ATTT ALT ATT = the deletion AT>A)
Truth alleles: the VCF records as they are (as vcfeval tags them), one per ALT.
Tensors: the GRCh38 allele of the label record (the summaries' projection) left-aligned and padded exactly as
pansoma_net_v2.linear_vcf does; tensors without a GRCh38 allele (off-reference, no unique GRCh38 node) are
`unplaced` and a PoN cannot reach them. Tensors are counted per (set, kind, label, reason).
"""
import json, os, sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "machine_learning"))
from filter_panel_of_normals import PON_NAMES, PON_RULES  # noqa: E402
from pansoma_net_v2.linear_vcf import Genome  # noqa: E402

pysam.set_verbosity(0)
HERE = Path(__file__).resolve().parent
DATA = Path("/scratch/jshen/data/pansoma_v2_tensors")
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
PON_DIR = Path("/scratch/jshen/data/Pansoma/panel_of_normal_VCFs")
PON_FILES = ["af-only-gnomad.hg38.vcf.gz", "Homo_sapiens_assembly38.dbsnp138.vcf.gz", "1000g_pon.hg38.vcf.gz",
             "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"]  # PON_NAMES order
TRUTH_VCF = {"HG008T": "/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz",
             "COLO829T": str(DATA / "COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz")}
SETS = ["HG008T_PacBio", "HG008T_ONT", "HG008T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT", "COLO829T_Illumina"]


def trimmed(pos, ref, alt):
    """Common suffix removed (at least one base kept on each side); position unchanged."""
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return pos, ref, alt


def main(chrom):
    genome = Genome(FASTA)
    (HERE / "scan").mkdir(exist_ok=True)
    keys = {}            # (pos, REF, ALT) -> index
    trim_keys = defaultdict(set)  # trimmed (pos, REF, ALT) -> indexes (truth only)

    def key_index(pos, ref, alt):
        k = (pos, ref.upper(), alt.upper())
        if k not in keys:
            keys[k] = len(keys)
        return keys[k]

    # ---- truth alleles, as the VCF records are
    truth = []
    for sample, path in TRUTH_VCF.items():
        with pysam.VariantFile(path) as vcf:
            if chrom not in vcf.header.contigs:
                continue
            for v in vcf.fetch(chrom):
                for alt in v.alts:
                    i = key_index(v.pos, v.ref, alt)
                    trim_keys[trimmed(v.pos, v.ref, alt)].add(i)
                    truth.append(dict(sample=sample, chrom=chrom, vcf_pos=v.pos, vcf_ref=v.ref, vcf_alt=alt, key=i))

    # ---- tensors: GRCh38 allele as linear_vcf builds it
    tensors = Counter()   # (set, kind, label, reason, key or -1) -> count
    ref_mismatch = Counter()
    for name in SETS:
        for kind in ("SNV", "INDEL"):
            f = DATA / name / "tensors" / kind / f"{chrom}_labels.ndjson"
            if not f.exists():
                continue
            with f.open() as stream:
                for line in stream:
                    r = json.loads(line)
                    lin = r.get("grch38")
                    group = (name, kind, r["label"], r["reason"])
                    if not lin:
                        tensors[group + (-1,)] += 1
                        continue
                    ev = r["candidate_id"].split(":")[2]
                    c, pos0, ref, alt = lin["chrom"], lin["pos0"], lin["ref"].upper(), lin["alt"].upper()
                    if genome.get(c, pos0, pos0 + len(ref)) != ref:
                        ref_mismatch[name] += 1
                        tensors[group + (-1,)] += 1
                        continue
                    if ev in ("INS", "DEL"):
                        pos0, ref, alt, _ = genome.left_align(c, pos0, ref, alt, ev)
                    pos, vref, valt = genome.vcf_allele(c, pos0, ref, alt, ev)
                    tensors[group + (key_index(pos, vref, valt),)] += 1

    # ---- stream every PoN once over the chromosome
    positions = {k[0] for k in keys} | {k[0] for k in trim_keys}
    n = len(keys)
    hit = {m: [0] * n for m in ("position", "allele", "rule")}  # bit per PoN
    af = [[None] * 4 for _ in range(n)]                        # best AF of an allele match (gnomAD, CoLoRSdb)
    trim_hit = defaultdict(int)                                # truth key index -> bit per PoN (trimmed allele)
    by_pos = defaultdict(list)
    for (pos, ref, alt), i in keys.items():
        by_pos[pos].append((ref, alt, i))
    streamed = {}
    for p, (pname, fname, rule) in enumerate(zip(PON_NAMES, PON_FILES, PON_RULES)):
        bit = 1 << p
        with pysam.VariantFile(str(PON_DIR / fname)) as pon:
            contig = chrom if chrom in pon.header.contigs else chrom[3:] if chrom[3:] in pon.header.contigs else None
            if contig is None:
                streamed[pname] = None
                continue
            count = 0
            for rec in pon.fetch(contig):
                count += 1
                if rec.pos not in positions:
                    continue
                rref = (rec.ref or "").upper()
                alts = [(a or "").upper() for a in (rec.alts or ())]
                values = rec.info.get("AF") if "AF" in rec.header.info else None
                values = values if isinstance(values, tuple) else (values,) * len(alts)
                somatic = rule.get("non_somatic") and rec.info.get("SAO") == 2
                for ref, alt, i in by_pos.get(rec.pos, ()):
                    hit["position"][i] |= bit
                    if ref != rref or alt not in alts:
                        continue
                    k = alts.index(alt)
                    hit["allele"][i] |= bit
                    a = values[k]
                    if a is not None and (af[i][p] is None or a > af[i][p]):
                        af[i][p] = a
                    min_af = rule.get("min_af")
                    if somatic or (min_af is not None and (a is None or a < min_af)):
                        continue
                    hit["rule"][i] |= bit
                for alt in alts:
                    for i in trim_keys.get(trimmed(rec.pos, rref, alt), ()):
                        trim_hit[i] |= bit
            streamed[pname] = count

    if os.environ.get("DUMP_KEYS"):  # per distinct allele (check against a vcfeval PoN output)
        with open(HERE / "scan" / f"{chrom}.keys.tsv", "w") as dump:
            dump.write("pos\tref\talt\tposition\tallele\trule\n")
            for (pos, ref, alt), i in keys.items():
                dump.write(f"{pos}\t{ref}\t{alt}\t{hit['position'][i]}\t{hit['allele'][i]}\t{hit['rule'][i]}\n")
    out = dict(chrom=chrom, pon_names=PON_NAMES, pon_records_streamed=streamed, ref_mismatch=dict(ref_mismatch),
               truth=[dict(t, position=hit["position"][t["key"]], allele=hit["allele"][t["key"]],
                           trimmed=trim_hit.get(t["key"], 0), rule=hit["rule"][t["key"]], af=af[t["key"]]) for t in truth],
               tensors=[[*g[:4], (hit["position"][g[4]], hit["allele"][g[4]], hit["rule"][g[4]]) if g[4] >= 0 else None, c]
                        for g, c in tensors.items()])
    # tensors: aggregate by (set, kind, label, reason, bits) -> count; None = unplaced
    agg = Counter()
    for name, kind, label, reason, bits, c in out["tensors"]:
        agg[(name, kind, label, reason, "unplaced" if bits is None else f"{bits[0]}:{bits[1]}:{bits[2]}")] += c
    out["tensors"] = [[*k, v] for k, v in sorted(agg.items())]
    for t in out["truth"]:
        t.pop("key")
    (HERE / "scan" / f"{chrom}.json").write_text(json.dumps(out))
    print(chrom, f"{len(truth):,} truth alleles, {sum(c for *_, c in out['tensors']):,} tensors, {n:,} distinct alleles,",
          "PoN records streamed", streamed, "REF mismatches", dict(ref_mismatch), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])

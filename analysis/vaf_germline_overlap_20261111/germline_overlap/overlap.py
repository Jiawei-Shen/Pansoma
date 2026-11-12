#!/usr/bin/env python3
"""Germline alleles of the three benchmark normals vs the alleles of the GRCh38-only, HPRC v1.1 (d9) and HPRC v2.1 (d46)
graphs, per chromosome (chr1-22), SNV and INDEL: normal-only, graph-only and shared allele counts.

Inputs are the allele tables of E7 (analysis/reviewer_response_20261011/experiments/E7_graph_composition, WORK =
/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition), where every allele is left-normalized against
GRCh38 and decomposed the same way (MNPs into SNVs; INS / DEL with one anchor base; key chrom:pos:REF:ALT):
  normal   WORK/sample_alleles/{HG008,COLO829,HCC1395}.germline.tsv.gz: the germline truth the labels use, HG008-N dipcall,
           COLO829BL dipcall, HCC1395BL DeepVariant hybrid QUAL >= 10; passed / in_bed from the label tables.
  graph    HPRC v1.1: WORK/hprc_alleles (all 89 haplotypes incl. CHM13, carrier count per allele); the d9 graph = alleles
           with >= 9 carriers. HPRC v2.1: WORK/v2.1/hprc_alleles (463 haplotypes); the d46 graph = >= 46 carriers.
           GRCh38-only: no non-reference allele (Linear_Reference_VG_Indexes/hg38_linear.giraffe.gbz is vg autoindex of
           the GRCh38 FASTA alone: one linear path per contig, its path covers run over the same nodes).
           Check for v1.1: the alleles of the real d9 graph's own deconstruct VCF (WORK/d9_alleles), reference
           'HPRC_v1.1_d9_VCF' (a lower bound: complex d9 ALTs are not decomposed there).
Strata: SNV; INDEL = INS / DEL of <= 50 bp. Not counted: alleles > 50 bp, COMPLEX, and normal alleles that decompose
into several primitives (listed in excluded.tsv). Counts are distinct allele keys, so the shared count is the same seen
from either side (as bcftools isec's 0002 / 0003).
Scopes: 'all' = every chr1-22 allele of both sides (as the earlier HG008N isec figure, whole VCFs);
        'bed' = normal alleles with FILTER PASS inside the normal's germline BED, graph alleles whose REF span lies
        inside that BED (one merged interval).
Outputs (next to this script): overlap_long.tsv; per (normal, graph, scope) R input
isec_<graph>_vs_<normal>.<scope>.per_chr.snv_indel.tsv with the columns of the earlier HG008N table
(A = graph, B = normal: chr, variant_type, A_only_count, B_only_count, overlap_count_Arep, overlap_count_Brep,
A_total_count, B_total_count), read by overlapping_HPRC_germline_HG008.R as AF-HPRC = A, HG008N = B; excluded.tsv.
usage: overlap.py [THREADS]   (Slurm: -c 8 --mem 48G)
"""
import csv
import gzip
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = Path("/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition")
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
NORMALS = {  # E7 sample -> (normal name, germline BED of the labels)
    "HG008": ("HG008N", "/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed"),
    "COLO829": ("COLO829BL", "/scratch/qfu/COLO829BL_DSA/dipcall_hg38/dipcall_hg38.dip.bed"),
    "HCC1395": ("HCC1395BL", "/scratch/jshen/data/CASTLE/HCC1395/truth/HCC1395_germline_high-confidence-regions_QUAL10.bed"),
}
GRAPHS = ("GRCh38", "HPRC_v1.1_d9", "HPRC_v2.1_d46", "HPRC_v1.1_d9_VCF")
SCOPES = ("all", "bed")
LONG = 50


def read_bed(path):
    iv = defaultdict(list)
    for line in open(path):
        if line.strip() and not line.startswith(("#", "track", "browser")):
            c, s, e = line.split("\t")[:3]
            iv[c].append((int(s), int(e)))
    out = {}
    for c, items in iv.items():
        merged = []
        for s, e in sorted(items):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = (np.array([m[0] for m in merged]), np.array([m[1] for m in merged]))
    return out


def key_span(k):
    """0-based half-open REF span of a key chrom:pos:REF:ALT (REF length; md5 keys are > 200 bp and never counted)."""
    _, pos, ref, _ = k.split(":", 3)
    return int(pos) - 1, int(pos) - 1 + len(ref)


def inside(bed, chrom, keys):
    """Boolean array: key's REF span inside one merged interval."""
    if chrom not in bed or not keys:
        return np.zeros(len(keys), bool)
    starts, ends = bed[chrom]
    sp = np.array([key_span(k) for k in keys], dtype=np.int64).reshape(-1, 2)
    i = np.searchsorted(starts, sp[:, 0], side="right") - 1
    ok = i >= 0
    return ok & (ends[np.maximum(i, 0)] >= sp[:, 1])


def stratum_from_key(kind, k):
    if kind == "SNV":
        return "SNV"
    if kind in ("INS", "DEL"):
        if "md5:" in k:  # an allele > 200 bp (E7 common.seqkey)
            return None
        _, _, ref, alt = k.split(":", 3)
        return None if abs(len(alt) - len(ref)) > LONG else "INDEL"
    return None


def load_normals():
    """NORM[sample][chrom][(scope, stratum)] = set of keys; EXCL counts."""
    norm, excl = {}, Counter()
    for s in NORMALS:
        d = defaultdict(lambda: defaultdict(set))
        with gzip.open(WORK / "sample_alleles" / f"{s}.germline.tsv.gz", "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                c = r["chrom"]
                if c not in AUTOSOMES:
                    continue
                keys = r["keys"].split(";")
                if len(keys) != 1:
                    excl[(s, "normal", "several primitives")] += 1
                    continue
                st = stratum_from_key(r["nkind"], keys[0])
                if st is None:
                    excl[(s, "normal", f"{r['nkind']} / > {LONG} bp")] += 1
                    continue
                d[c][("all", st)].add(keys[0])
                if r["passed"] == "T" and r["in_bed"] == "T":
                    d[c][("bed", st)].add(keys[0])
        norm[s] = d
        print(f"[INFO] {s}: " + ", ".join(f"{sc}/{st} {sum(len(d[c][(sc, st)]) for c in d):,}"
                                           for sc in SCOPES for st in ("SNV", "INDEL")), flush=True)
    return norm, excl


def graph_keys(chrom):
    """{graph: {stratum: list of keys}} for one chromosome; Counter of excluded graph alleles."""
    out, excl = {g: {"SNV": [], "INDEL": []} for g in GRAPHS}, Counter()
    for g, path, d in (("HPRC_v1.1_d9", WORK / "hprc_alleles" / f"{chrom}.tsv.gz", 9),
                       ("HPRC_v2.1_d46", WORK / "v2.1" / "hprc_alleles" / f"{chrom}.tsv.gz", 46)):
        with gzip.open(path, "rt") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            ik, ikind, irc, inc = (hdr.index(x) for x in ("key", "kind", "rep_class", "n_carriers"))
            for line in f:
                x = line.split("\t", inc + 1)
                if int(x[inc]) < d:
                    continue
                kind = x[ikind]
                if kind == "SNV":
                    out[g]["SNV"].append(x[ik])
                elif kind in ("INS", "DEL") and x[irc] != "SV":
                    out[g]["INDEL"].append(x[ik])
                else:
                    excl[(g, "graph", "COMPLEX" if kind == "COMPLEX" else f"> {LONG} bp")] += 1
    with gzip.open(WORK / "d9_alleles" / f"{chrom}.tsv.gz", "rt") as f:
        next(f)
        for line in f:
            k, kind = line.split("\t")[:2]
            st = stratum_from_key(kind, k)
            if st is None:
                excl[("HPRC_v1.1_d9_VCF", "graph", kind if kind == "COMPLEX" else f"> {LONG} bp")] += 1
            else:
                out["HPRC_v1.1_d9_VCF"][st].append(k)
    for g in GRAPHS:
        for st in out[g]:
            out[g][st] = sorted(set(out[g][st]))
    return out, excl


def work(chrom):
    gk, excl = graph_keys(chrom)
    rows = []
    for s, (name, _) in NORMALS.items():
        bed = BEDS[s]
        for g in GRAPHS:
            for st in ("SNV", "INDEL"):
                keys = gk[g][st]
                in_b = inside(bed, chrom, keys)
                for sc in SCOPES:
                    gset = set(keys) if sc == "all" else {k for k, b in zip(keys, in_b) if b}
                    nset = NORM[s][chrom][(sc, st)]
                    ov = len(nset & gset)
                    rows.append(dict(normal=name, graph=g, scope=sc, chr=chrom, variant_type=st, normal_total=len(nset),
                                     graph_total=len(gset), overlap=ov, normal_only=len(nset) - ov,
                                     graph_only=len(gset) - ov))
    print(f"[INFO] {chrom} done", flush=True)
    return rows, excl


def main():
    global NORM, BEDS
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    NORM, excl = load_normals()
    BEDS = {s: read_bed(b) for s, (_, b) in NORMALS.items()}
    rows = []
    with Pool(threads) as pool:
        for r, e in pool.imap_unordered(work, AUTOSOMES):
            rows += r
            excl.update(e)
    order = {c: i for i, c in enumerate(AUTOSOMES)}
    rows.sort(key=lambda r: (r["normal"], r["graph"], r["scope"], order[r["chr"]], r["variant_type"] != "SNV"))
    cols = ["normal", "graph", "scope", "chr", "variant_type", "normal_total", "graph_total", "overlap", "normal_only",
            "graph_only"]
    with open(HERE / "overlap_long.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    rdir = HERE / "R_input"
    rdir.mkdir(exist_ok=True)
    by = defaultdict(list)
    for r in rows:
        by[(r["normal"], r["graph"], r["scope"])].append(r)
    for (name, g, sc), rr in by.items():
        with open(rdir / f"isec_{g}_vs_{name}.{sc}.per_chr.snv_indel.tsv", "w") as f:
            f.write("chr\tvariant_type\tA_only_count\tB_only_count\toverlap_count_Arep\toverlap_count_Brep\t"
                    "A_total_count\tB_total_count\n")
            for r in rr:
                f.write(f"{r['chr']}\t{r['variant_type']}\t{r['graph_only']}\t{r['normal_only']}\t{r['overlap']}\t"
                        f"{r['overlap']}\t{r['graph_total']}\t{r['normal_total']}\n")
    with open(HERE / "excluded.tsv", "w") as f:
        f.write("source\tside\treason\tcount (chr1-22, distinct keys for graphs, records for normals)\n")
        for (src, side, why), n in sorted(excl.items()):
            f.write(f"{src}\t{side}\t{why}\t{n}\n")
    print(f"[INFO] wrote {len(rows)} rows, {len(by)} R input files", flush=True)


if __name__ == "__main__":
    main()

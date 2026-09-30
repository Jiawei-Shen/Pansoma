"""Shared helpers for the v1 vs v2 head-to-head (read-only on all inputs)."""
import gzip
import json
import bisect
from collections import defaultdict

V2 = "/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors"
TRUTH_TSV = "/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/truth/somatic.graph.tsv"
RECALL_TSV = V2 + "/somatic.recall.tsv"
SOM_BED = "/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed"
GERM_BED = "/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed"
COMP = str.maketrans("ACGTN", "TGCAN")


def read_bed(path, chrom):
    iv = []
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            p = line.split("\t")
            if p[0] == chrom:
                iv.append((int(p[1]), int(p[2])))
    iv.sort()
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def intersect(a, b):
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append([s, e])
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


class Region:
    def __init__(self, iv):
        self.starts = [s for s, _ in iv]
        self.ends = [e for _, e in iv]

    def contains0(self, pos0):
        k = bisect.bisect_right(self.starts, pos0) - 1
        return k >= 0 and pos0 < self.ends[k]


def confident_chr1():
    return Region(intersect(read_bed(SOM_BED, "chr1"), read_bed(GERM_BED, "chr1")))


def somatic_bed_chr1():
    return Region(read_bed(SOM_BED, "chr1"))


def truth_697():
    """chr1 SNP truth alleles that PASS and are in the somatic BED (somatic.recall.tsv), with pos0/ref/alt."""
    rec = {}
    with open(RECALL_TSV) as f:
        next(f)
        for line in f:
            tid, chrom, vpos, kind, passed, in_bed, status, *_ = line.rstrip("\n").split("\t")
            if chrom == "chr1" and kind == "SNP" and passed == "True" and in_bed == "True":
                rec[int(tid)] = dict(vcf_pos=int(vpos), status=status)
    with open(TRUTH_TSV) as f:
        hdr = next(f).rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            t = int(r["truth_id"])
            if t in rec:
                rec[t].update(pos0=int(r["pos0"]), ref=r["ref"], alt=r["alt"])
    return rec

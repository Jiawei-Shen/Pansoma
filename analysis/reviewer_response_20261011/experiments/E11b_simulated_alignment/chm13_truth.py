#!/usr/bin/env python3
"""E11b CHM13 control: lift every simulated read's haplotype template to T2T-CHM13.

Needs build_segments.py run on the haplotype -> CHM13 alignments (chm13_paf.sbatch outputs):
  python build_segments.py --dipcall-prefix <D>/chm13_truth/HG008N_CHM13_dipcall --out <D>/chm13_truth
A read gets a CHM13 truth (ok = 1) when its whole template lies in one CHM13 one-to-one segment of the same haplotype
contig. Output: <D>/chm13_truth/<platform>.chm13_truth.npz (chrom, ref_start, ref_end, ref_strand, ok; rows = truth.npz).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e11b_common as C  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment")
    ap.add_argument("--platform", required=True)
    a = ap.parse_args()
    T = np.load(f"{a.base}/sim/{a.platform}/truth.npz")
    g38 = C.read_tsv(f"{a.base}/sim/segments.tsv")
    contig_of_seg = np.array([r["contig"] for r in g38])
    al = C.Alignments(f"{a.base}/chm13_truth")
    n = len(T["read"])
    hs, he, hstr, hap = T["hap_start"], T["hap_end"], T["hap_strand"], T["hap"]
    contig = contig_of_seg[T["seg"]]
    aln = np.full(n, -1, np.int64)
    by_contig = {}
    for s in al.segs:
        by_contig.setdefault(s["contig"], []).append((s["qstart"], s["qend"], s["aln"]))
    for c, lst in by_contig.items():
        lst.sort()
        st = np.array([x[0] for x in lst])
        en = np.array([x[1] for x in lst])
        ai = np.array([x[2] for x in lst])
        m = np.flatnonzero(contig == c)
        k = np.searchsorted(st, hs[m], side="right") - 1
        okk = (k >= 0)
        kk = np.maximum(k, 0)
        inside = okk & (he[m] <= en[kk])
        aln[m[inside]] = ai[kk[inside]]
    ok = aln >= 0
    chrom = np.full(n, -1, np.int8)
    rs = np.full(n, -1, np.int64)
    re_ = np.full(n, -1, np.int64)
    rstr = np.zeros(n, np.int8)
    if ok.any():
        c2, s2, e2, st2, _ = al.lift_reads(aln[ok], hs[ok], he[ok], hstr[ok])
        chrom[ok], rs[ok], re_[ok], rstr[ok] = c2, s2, e2, st2
    np.savez_compressed(f"{a.base}/chm13_truth/{a.platform}.chm13_truth.npz", chrom=chrom, ref_start=rs, ref_end=re_,
                        ref_strand=rstr, ok=ok)
    print(a.platform, "reads", n, "with CHM13 truth", int(ok.sum()), f"({ok.mean():.4%})")


if __name__ == "__main__":
    main()

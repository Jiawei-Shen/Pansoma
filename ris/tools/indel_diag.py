"""INDEL call diagnostic on a run's rtg curve (raw_curve: every PASS + LowQual call, each TP or FP).

For each call: p_somatic (QUAL), FORMAT AF, the max population AF of an exact allele match in gnomAD and CoLoRSdb, and
the homopolymer run at the event on GRCh38. Then the F1 that a decision (p threshold, PoN AF floor, read-AF floor)
would give, counted on calls (a TP call = one truth found; truth = the baseline of the curve run).
For exploration only: on the test chromosome this is tuning on the test set.

    py.sh tools/indel_diag.py RUN_DIR [RUN_DIR ...] --out TSV_PREFIX
"""
import argparse
import gzip
import math
import os
import sys

import numpy as np
import pysam

PON = "/storage2/fs1/epigenome/Active/jiawei.s/Dataset/panel_of_normal_VCFs"
FASTA = "/storage2/fs1/epigenome/Active/jiawei.s/Dataset/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"


def pon_af(vcf, chrom, pos, ref, alt):
    """Max AF of an exact (pos, REF, ALT) match, 0 if none."""
    best = 0.0
    for c in vcf.fetch(chrom, max(0, pos - 1), pos):
        if c.pos != pos or (c.ref or "").upper() != ref.upper():
            continue
        alts = [(a or "").upper() for a in (c.alts or ())]
        af = c.info.get("AF")
        af = af if isinstance(af, tuple) else (af,) * len(alts)
        for k, a in enumerate(alts):
            if a == alt.upper() and af[k] is not None:  # the PoN script ignores AF=None
                best = max(best, af[k])
    return best


def homopolymer(fa, chrom, pos, ref, alt):
    """Longest run of the event's base right after the padding base (VCF pos is the padding base)."""
    ev = alt[1:] if len(alt) > len(ref) else ref[1:]
    if not ev:
        return 0
    b = ev[0]
    if any(x != b for x in ev):
        return 0
    seq = fa.fetch(chrom, pos, pos + 60).upper()
    n = 0
    for x in seq:
        if x != b:
            break
        n += 1
    return n


def calls(path):
    out = []
    with pysam.VariantFile(path) as v:
        for r in v:
            p = r.info.get("P_SOMATIC")
            s = r.samples[0]
            af = s.get("AF")
            af = af[0] if isinstance(af, tuple) else af
            out.append((r.chrom, r.pos, r.ref, r.alts[0], float(p) if p is not None else float("nan"),
                        float(af) if af is not None else float("nan")))
    return out


def count(path):
    with pysam.VariantFile(path) as v:
        return sum(1 for _ in v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--eval", default="eval_bed")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    gn = pysam.VariantFile(f"{PON}/af-only-gnomad.hg38.vcf.gz")
    co = pysam.VariantFile(f"{PON}/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")
    fa = pysam.FastaFile(FASTA)
    for run in a.runs:
        name = os.path.basename(run.rstrip("/"))
        d = f"{run}/vcf_chr1/{a.eval}/vcfeval/raw_curve"
        n_truth = count(f"{d}/tp-baseline.vcf.gz") + count(f"{d}/fn.vcf.gz")
        rows = [(c, 1) for c in calls(f"{d}/tp.vcf.gz")] + [(c, 0) for c in calls(f"{d}/fp.vcf.gz")]
        tab = []
        for (chrom, pos, ref, alt, p, af), tp in rows:
            g, c = pon_af(gn, chrom, pos, ref, alt), pon_af(co, chrom, pos, ref, alt)
            tab.append((p, af, max(g, c), homopolymer(fa, chrom, pos, ref, alt), tp))
        t = np.array(tab, dtype=float)
        with gzip.open(f"{a.out}.{name}.tsv.gz", "wt") as f:
            f.write("p_somatic\tread_af\tpon_af\thomopolymer\ttp\n")
            for r in t:
                f.write("\t".join(f"{x:.6g}" for x in r) + "\n")
        p, raf, paf, hp, tp = t.T
        print(f"=== {name}: truth {n_truth}, curve calls {len(t)} (TP {int(tp.sum())}, FP {int((1 - tp).sum())})")
        # where TP and FP sit
        for lab, m in (("TP", tp == 1), ("FP", tp == 0)):
            print(f"  {lab}: pon_af=0 {np.mean(paf[m] == 0):.2f}  <0.01 {np.mean((paf[m] > 0) & (paf[m] < 0.01)):.2f}"
                  f"  0.01-0.05 {np.mean((paf[m] >= 0.01) & (paf[m] < 0.05)):.2f}  0.05-0.2 {np.mean((paf[m] >= 0.05) & (paf[m] < 0.2)):.2f}"
                  f"  >=0.2 {np.mean(paf[m] >= 0.2):.2f} | read AF median {np.nanmedian(raf[m]):.3f}"
                  f" | homopolymer>=6 {np.mean(hp[m] >= 6):.2f} | p median {np.nanmedian(p[m]):.3f}")
        best = []
        for floor in (0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, math.inf):
            for rmin in (0.0, 0.05, 0.1, 0.15, 0.2):
                keep0 = (paf < floor) & (raf >= rmin)
                for thr in np.unique(np.round(p[keep0], 3)):
                    keep = keep0 & (p >= thr)
                    n = keep.sum()
                    if n == 0:
                        continue
                    tpk = tp[keep].sum()
                    prec, rec = tpk / n, tpk / n_truth
                    f1 = 2 * prec * rec / (prec + rec) if tpk else 0.0
                    best.append((f1, floor, rmin, thr, prec, rec, int(n)))
        best.sort(reverse=True)
        print("  best (F1, PoN floor, read-AF min, p thr, P, R, calls):")
        for b in best[:3]:
            print("   ", " ".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in b))
        for floor in (0.01, 0.05, 0.1, math.inf):
            bf = max((b for b in best if b[1] == floor and b[2] == 0.0), default=None)
            print(f"  PoN floor {floor}: best F1 {bf[0]:.3f} (thr {bf[3]:.3f}, P {bf[4]:.3f}, R {bf[5]:.3f}, {bf[6]} calls)" if bf else "")
        sys.stdout.flush()


if __name__ == "__main__":
    main()

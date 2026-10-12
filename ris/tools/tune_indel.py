"""Choose the INDEL decision on a tuning chromosome (not the test one): p_somatic threshold, PoN AF floor (gnomAD /
CoLoRSdb exact allele) and read-AF floor, by the F1 of the calls against the truth.

    python ris/tools/tune_indel.py EVAL_DIR OUT.json [--floors F ...] [--read-af A ...]
EVAL_DIR = a vcfeval.py --output of the tuning chromosome. The general decision tunes the threshold only:
--floors 0.01 --read-af 0 (one flat 0.01 floor for every INDEL, no read-AF floor; the PoN script's own INDEL rule
splits 0.01 / 0.001 by homopolymer, which this tool does not apply).

Uses EVAL_DIR/vcfeval/raw_curve (every PASS or LowQual call, TP or FP, from graph_vcf --min-score up): a call is kept
when p >= threshold, read AF >= floor and its population AF < PoN floor; precision = kept TP / kept, recall = kept TP /
truth alleles (one TP call per truth allele: rtg's call-side TP count)."""
import argparse
import json
import math
import sys

import numpy as np
import pysam

from indel_diag import PON, calls, count, pon_af

FLOORS = (0.005, 0.01, 0.02, 0.05)
READ_AF = (0.0, 0.1, 0.15, 0.2, 0.25)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("eval_dir")
    ap.add_argument("out")
    ap.add_argument("--floors", nargs="+", type=float, default=list(FLOORS))
    ap.add_argument("--read-af", nargs="+", type=float, default=list(READ_AF))
    a = ap.parse_args()
    eval_dir, out = a.eval_dir, a.out
    d = f"{eval_dir}/vcfeval/raw_curve"
    n_truth = count(f"{d}/tp-baseline.vcf.gz") + count(f"{d}/fn.vcf.gz")
    gn = pysam.VariantFile(f"{PON}/af-only-gnomad.hg38.vcf.gz")
    co = pysam.VariantFile(f"{PON}/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")
    rows = [(c, 1) for c in calls(f"{d}/tp.vcf.gz")] + [(c, 0) for c in calls(f"{d}/fp.vcf.gz")]
    t = np.array([(p, af, max(pon_af(gn, ch, pos, ref, alt), pon_af(co, ch, pos, ref, alt)), tp)
                  for (ch, pos, ref, alt, p, af), tp in rows], dtype=float)
    p, raf, paf, tp = t.T
    best = None
    for floor in a.floors:
        for rmin in a.read_af:
            keep0 = ~(paf >= floor * (1 - 1e-6)) & (raf >= rmin)  # the PoN script's float32 tolerance
            for thr in np.unique(np.round(p[keep0], 3)):
                keep = keep0 & (p >= thr)
                n, k = int(keep.sum()), float(tp[keep].sum())
                if not k:
                    continue
                prec, rec = k / n, k / n_truth
                f1 = 2 * prec * rec / (prec + rec)
                cand = (f1, -n, floor, rmin, float(thr), prec, rec, n)
                if best is None or cand[:2] > best[:2]:
                    best = cand
    f1, _, floor, rmin, thr, prec, rec, n = best
    res = dict(eval_dir=eval_dir, truth=n_truth, curve_calls=len(t), threshold=thr, pon_floor=floor, read_af_min=rmin,
               tuning_f1=f1, tuning_precision=prec, tuning_recall=rec, tuning_calls=n)
    json.dump(res, open(out, "w"), indent=2)
    print(json.dumps(res))


if __name__ == "__main__":
    main()

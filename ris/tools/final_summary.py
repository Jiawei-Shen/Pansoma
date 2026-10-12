"""Final chr1 table of the INDEL rescorer runs (exp_rescore.sh without --chr2-only), with a 5 Mb block bootstrap 95%
interval of the PoN F1 at the threshold. Standard library only (login node).

    python3 ris/tools/final_summary.py RESCORE_RUN [RESCORE_RUN ...] [--boot 2000]

For each run, set and variant (pooled, loso_<sample>): rtg's pon row in the BED (calls, TP, FP, FN, P, R, F1) and
without the BED (F1). The bootstrap resamples 5 Mb blocks of chr1 with replacement and recomputes
F1 = 2PR/(P+R), P = TP_call/(TP_call+FP), R = TP_base/(TP_base+FN), from vcfeval/pon_calls/{tp,fp,tp-baseline,fn}.vcf.gz."""
import glob
import gzip
import os
import random
import re
import sys

RUNS = "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"
ORDER = ["HG008T_Illumina", "HG008T_PacBio", "HG008T_ONT", "COLO829T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT"]
BLOCK = 5_000_000


def positions(path):
    out = []
    if not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                out.append(int(line.split("\t", 2)[1]))
    return out


def block_counts(d):
    c = {}
    for name in ("tp", "fp", "tp-baseline", "fn"):
        for p in positions(f"{d}/{name}.vcf.gz"):
            b = p // BLOCK
            c.setdefault(b, {"tp": 0, "fp": 0, "tp-baseline": 0, "fn": 0})[name] += 1
    return c


def f1(t):
    p = t["tp"] / (t["tp"] + t["fp"]) if t["tp"] + t["fp"] else 0.0
    r = t["tp-baseline"] / (t["tp-baseline"] + t["fn"]) if t["tp-baseline"] + t["fn"] else 0.0
    return 2 * p * r / (p + r) if p + r else 0.0


def ci(d, n, rng):
    c = block_counts(d)
    keys = list(c)
    if not keys:
        return None
    vals = []
    for _ in range(n):
        t = {"tp": 0, "fp": 0, "tp-baseline": 0, "fn": 0}
        for _ in keys:
            b = c[rng.choice(keys)]
            for k in t:
                t[k] += b[k]
        vals.append(f1(t))
    vals.sort()
    return vals[int(0.025 * n)], vals[int(0.975 * n) - 1]


def pon_row(report):
    if not os.path.exists(report):
        return None
    m = re.search(r"^pon\s+(.*)$", open(report).read(), re.M)
    if not m:
        return None
    v = m.group(1).split()
    return dict(calls=v[0], tp=v[1], fp=v[2], fn=v[3], p=float(v[4]), r=float(v[5]), f1=float(v[6]))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n = int(sys.argv[sys.argv.index("--boot") + 1]) if "--boot" in sys.argv else 2000
    args = [a for a in args if not a.isdigit()]
    rng = random.Random(0)
    print(f"{'run':14s} {'variant':14s} {'set':18s} {'calls':>6s} {'TP':>4s} {'FP':>5s} {'FN':>4s} {'P':>6s} {'R':>6s} "
          f"{'F1':>6s} {'95% CI':>14s} {'noBED F1':>8s}")
    for run in args:
        base = f"{RUNS}/{run}"
        def key(x):
            v = x.rsplit("_chr1_", 1)[1]
            si = next((i for i, st in enumerate(ORDER) if f"__{st}_chr1_" in x), 9)
            return (v != "pooled", v, si)
        for d in sorted(glob.glob(f"{base}/vcf_*_chr1_*"), key=key):
            variant = d.rsplit("_chr1_", 1)[1]
            s = next((x for x in ORDER if f"__{x}_chr1_" in d), "?")
            r = pon_row(f"{d}/eval_bed/report.txt")
            nb = pon_row(f"{d}/eval_nobed/report.txt")
            if not r:
                print(f"{run:14s} {variant:14s} {s:18s} (missing)")
                continue
            nob = f"{nb['f1']:.3f}" if nb else "-"
            lo_hi = ci(f"{d}/eval_bed/vcfeval/pon_calls", n, rng)
            cis = f"[{lo_hi[0]:.3f}, {lo_hi[1]:.3f}]" if lo_hi else "-"
            print(f"{run:14s} {variant:14s} {s:18s} {r['calls']:>6s} {r['tp']:>4s} {r['fp']:>5s} {r['fn']:>4s} "
                  f"{r['p']:6.3f} {r['r']:6.3f} {r['f1']:6.3f} {cis:>14s} {nob:>8s}")


if __name__ == "__main__":
    main()

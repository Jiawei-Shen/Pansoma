import gzip, sys
from collections import Counter
p = sys.argv[1]
R = f"/scratch/jshen/data/pansoma_net_v2_runs/multicancer6_{p}_SNV/vcf_chr1_HCC1395T_{p}_d46_normfix/eval_bed"
E8 = "/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score"
def recs(f):
    out = {}
    for l in gzip.open(f, "rt"):
        if l[0] == "#": continue
        x = l.rstrip("\n").split("\t")
        if x[0] == "chr1" and len(x[3]) == 1 and len(x[4]) == 1: out[(int(x[1]), x[3], x[4])] = x
    return out
calls = recs(f"{R}/HCC1395T_{p}_d46.SNV.linear.pon.vcf.gz")
af = {k: float(dict(zip(x[8].split(":"), x[9].split(":")))["AF"]) for k, x in calls.items()}
truth = set(recs(f"{R}/vcfeval/pon_calls/tp-baseline.vcf.gz")) | set(recs(f"{R}/vcfeval/pon_calls/fn.vcf.gz"))
found = {"Pansoma": set(recs(f"{R}/vcfeval/pon_calls/tp-baseline.vcf.gz"))}
for c, n in (("deepsomatic", "DeepSomatic-TO"), ("clairs_to", "ClairS-TO")):
    found[n] = set(recs(f"{E8}/{c}/HCC1395T_{p}/SNV/eval_bed/pass_pon/tp-baseline.vcf.gz"))
bins = [(0.08, 0.25), (0.25, 0.4), (0.4, 0.6), (0.6, 0.85), (0.85, 1.01)]
print(f"{p}: truth SNVs with a Pansoma tensor in the VCF, by Pansoma read AF; recall per caller")
print("| AF | truth | " + " | ".join(found) + " |")
for lo, hi in bins:
    t = [k for k in truth if k in af and lo <= af[k] < hi]
    print(f"| [{lo}, {hi}) | {len(t)} | " + " | ".join(f"{sum(k in f for k in t)/len(t):.3f}" for f in found.values()) + " |")
t = [k for k in truth if k not in af]
print(f"| no tensor in the VCF | {len(t)} | " + " | ".join(f"{sum(k in f for k in t)/len(t):.3f}" for f in found.values()) + " |")

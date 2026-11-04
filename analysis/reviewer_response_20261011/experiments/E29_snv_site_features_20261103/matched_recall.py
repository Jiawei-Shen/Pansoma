# usage: python matched_recall.py PLAT  -- diagnostic, not part of the pre-specified evaluation: the site-feature model's
# HCC1395T recall by read AF at the p_somatic cutoff where its overall recall equals the normfix model's PASS recall
# (vcfeval pon_curve: every scored call after the PoN, BED), and the precision there.
import gzip, sys
p = sys.argv[1]; R = "/scratch/jshen/data/pansoma_net_v2_runs"
def recs(f, chrom):
    out = {}
    for l in gzip.open(f, "rt"):
        if l[0] == "#": continue
        x = l.rstrip("\n").split("\t")
        if x[0] == chrom and len(x[3]) == 1 and len(x[4]) == 1: out[(int(x[1]), x[3], x[4])] = x
    return out
ps = lambda x: float(dict(kv.split("=") for kv in x[7].split(";") if "=" in kv)["P_SOMATIC"])
for chrom in ("chr20", "chr1"):
    s = f"HCC1395T_{p}_d46"
    O, N = f"{R}/multicancer6_{p}_SNV/vcf_{chrom}_{s}_normfix/eval_bed", f"{R}/multicancer6_{p}_SNV_sitef_ppr9/vcf_{chrom}_{s}/eval_bed"
    calls = recs(f"{O}/{s}.SNV.linear.pon.vcf.gz", chrom)
    af = {k: float(dict(zip(x[8].split(":"), x[9].split(":")))["AF"]) for k, x in calls.items()}
    truth = set(recs(f"{O}/vcfeval/pon_calls/tp-baseline.vcf.gz", chrom)) | set(recs(f"{O}/vcfeval/pon_calls/fn.vcf.gz", chrom))
    old = set(recs(f"{O}/vcfeval/pon_calls/tp-baseline.vcf.gz", chrom))
    tp = {k: ps(x) for k, x in recs(f"{N}/vcfeval/pon_curve/tp.vcf.gz", chrom).items()}
    fp = [ps(x) for x in recs(f"{N}/vcfeval/pon_curve/fp.vcf.gz", chrom).values()]
    otp = {k: ps(x) for k, x in recs(f"{O}/vcfeval/pon_curve/tp.vcf.gz", chrom).items()}
    ofp = [ps(x) for x in recs(f"{O}/vcfeval/pon_curve/fp.vcf.gz", chrom).values()]
    target = len(old & truth) / len(truth)
    sc = sorted((v for k, v in tp.items() if k in truth), reverse=True)
    t = sc[min(len(sc), round(target * len(truth))) - 1]
    new = {k for k, v in tp.items() if k in truth and v >= t}
    prec = lambda tps, fps, c: sum(v >= c for v in tps) / max(1, sum(v >= c for v in tps) + sum(v >= c for v in fps))
    print(f"## {p} HCC1395T {chrom}: recall {target:.3f} (normfix PASS, p >= 0.622) vs site features at p >= {t:.3f} (recall {len(new)/len(truth):.3f})")
    print(f"precision there: normfix {prec(list(otp.values()), ofp, 0.6224303245544434):.3f}, site features {prec(list(tp.values()), fp, t):.3f}\n")
    print("| AF | truth | site features, matched recall | normfix PASS |"); print("|---|---|---|---|")
    for lo, hi in [(0.08, 0.25), (0.25, 0.4), (0.4, 0.6), (0.6, 0.85), (0.85, 1.01)]:
        tt = [k for k in truth if k in af and lo <= af[k] < hi]
        print(f"| [{lo}, {min(hi, 1.0)}] | {len(tt)} | {sum(k in new for k in tt)/len(tt):.3f} | {sum(k in old for k in tt)/len(tt):.3f} |")
    print()

# usage: python at_cutoffs.py PLAT -- test P / R / F1 after the PoN (BED) of the site-feature model at the PASS cutoffs its
# checkpoint's validation gives (post-PoN callable recall 0.9 / 0.95, post-PoN best F1), from the vcfeval pon_curve;
# then HCC1395T recall by read AF at each cutoff (bins and truth of compare.py).
import gzip, json, sys
p = sys.argv[1]; R = "/scratch/jshen/data/pansoma_net_v2_runs"; RUN = f"{R}/multicancer6_{p}_SNV_sitef_ppr9"
CO = "COLO829T_fiberseq" if p == "PacBio" else f"COLO829T_{p}"
L = [json.loads(l) for l in open(f"{RUN}/metrics.jsonl")]
best = max(L, key=lambda r: r["val"]["postpon"]["f1"]); v = best["val"]
assert f"best epoch {best['epoch']} " in open(f"{RUN}/train.log").read(), "best epoch differs from train.log"
cut = {"recall 0.9 (ppr9)": v["somatic_at_recall_postpon_callable"]["0.9"]["threshold"],
       "recall 0.95": v["somatic_at_recall_postpon_callable"]["0.95"]["threshold"],
       "best F1": v["somatic_threshold_postpon"]}
print(f"{p} epoch {best['epoch']}: validation cutoffs " + ", ".join(f"{k} p >= {t:.3f}" for k, t in cut.items()))
ps = lambda l: float(l.split("P_SOMATIC=")[1].split(";")[0])
def scores(f): return [ps(l) for l in gzip.open(f, "rt") if l[0] != "#"]
print("| test | " + " | ".join(cut) + " |"); print("|---" * (len(cut) + 1) + "|")
for name, s, chrom in [("HCC1395T chr20", f"HCC1395T_{p}", "chr20"), ("HCC1395T chr1", f"HCC1395T_{p}", "chr1"),
                       ("HG008T chr1", f"HG008T_{p}", "chr1"), ("COLO829T chr1", CO, "chr1")]:
    E = f"{RUN}/vcf_{chrom}_{s}_d46/eval_bed"; a = json.load(open(f"{E}/report.json"))["pon"]["at_threshold"]
    n = a["tp_baseline"] + a["fn"]; tp, fp = scores(f"{E}/vcfeval/pon_curve/tp.vcf.gz"), scores(f"{E}/vcfeval/pon_curve/fp.vcf.gz")
    row = []
    for t in cut.values():
        T, F = sum(x >= t for x in tp), sum(x >= t for x in fp); P, Rc = T / (T + F), T / n
        row.append(f"{2*P*Rc/(P+Rc):.3f} ({P:.3f} / {Rc:.3f})")
    print(f"| {name} | " + " | ".join(row) + " |")

def recs(f, chrom):
    out = {}
    for l in gzip.open(f, "rt"):
        if l[0] == "#": continue
        x = l.rstrip("\n").split("\t")
        if x[0] == chrom and len(x[3]) == 1 and len(x[4]) == 1: out[(int(x[1]), x[3], x[4])] = x
    return out
for chrom in ("chr20", "chr1"):
    s = f"HCC1395T_{p}_d46"
    O, N = f"{R}/multicancer6_{p}_SNV/vcf_{chrom}_{s}_normfix/eval_bed", f"{RUN}/vcf_{chrom}_{s}/eval_bed"
    af = {k: float(dict(zip(x[8].split(":"), x[9].split(":")))["AF"]) for k, x in recs(f"{O}/{s}.SNV.linear.pon.vcf.gz", chrom).items()}
    truth = set(recs(f"{O}/vcfeval/pon_calls/tp-baseline.vcf.gz", chrom)) | set(recs(f"{O}/vcfeval/pon_calls/fn.vcf.gz", chrom))
    old = set(recs(f"{O}/vcfeval/pon_calls/tp-baseline.vcf.gz", chrom))
    tp = {k: ps(x[7] + ";") for k, x in recs(f"{N}/vcfeval/pon_curve/tp.vcf.gz", chrom).items()}
    print(f"\nHCC1395T {chrom} recall by read AF\n\n| AF | truth | " + " | ".join(cut) + " | normfix PASS |"); print("|---" * (len(cut) + 3) + "|")
    for lo, hi in [(0.08, 0.25), (0.25, 0.4), (0.4, 0.6), (0.6, 0.85), (0.85, 1.01)]:
        tt = [k for k in truth if k in af and lo <= af[k] < hi]
        if tt: print(f"| [{lo}, {min(hi, 1.0)}] | {len(tt)} | " + " | ".join(f"{sum(tp.get(k, -1) >= t for k in tt)/len(tt):.3f}" for t in cut.values()) + f" | {sum(k in old for k in tt)/len(tt):.3f} |")

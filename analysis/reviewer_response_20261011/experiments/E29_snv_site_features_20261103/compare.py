# usage: python compare.py PLAT   (Illumina | PacBio | ONT)
# The pre-specified evaluation of README.md: the site-feature model (multicancer6_<PLAT>_SNV_sitef_ppr9) against
# multicancer6_<PLAT>_SNV (normfix predictions) and the comparators (E8), SNVs after the PoN in the benchmark BED, PASS at
# each model's own rule; then HCC1395T recall by read AF (bins of E27 v3/hcc1395_snv_gap/recall_by_af.py).
import gzip, json, os, sys
p = sys.argv[1]
R = "/scratch/jshen/data/pansoma_net_v2_runs"
E8 = "/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score"
NEW, OLD = f"{R}/multicancer6_{p}_SNV_sitef_ppr9", f"{R}/multicancer6_{p}_SNV"
CO = "COLO829T_fiberseq" if p == "PacBio" else f"COLO829T_{p}"
TESTS = [("HCC1395T chr20", f"HCC1395T_{p}", "chr20"), ("HCC1395T chr1", f"HCC1395T_{p}", "chr1"),
         ("HG008T chr1", f"HG008T_{p}", "chr1"), ("COLO829T chr1", CO, "chr1")]

def pansoma(run, s, chrom, suffix):
    f = f"{run}/vcf_{chrom}_{s}_d46{suffix}/eval_bed/report.json"
    if not os.path.exists(f): return None
    d = json.load(open(f)); a = d["pon"]["at_threshold"]
    return a["precision"], a["recall"], a["f1"], d["threshold"]

def comparator(c, s, chrom):
    f = f"{E8}/{c}/{s}{'_chr20' if chrom == 'chr20' else ''}/SNV/report.json"
    if not os.path.exists(f): return None
    a = json.load(open(f))["eval_bed"]["pass_pon"]["at_threshold"]
    return a["precision"], a["recall"], a["f1"], None

fmt = lambda r: "-" if r is None else f"{r[2]:.3f} ({r[0]:.3f} / {r[1]:.3f})"
print(f"## {p}: SNV F1 (P / R) after the PoN, benchmark BED, PASS at each caller's own rule\n")
print("| test | Pansoma + site features (ppr9) | Pansoma (normfix) | DeepSomatic-TO | ClairS-TO |")
print("|---|---|---|---|---|")
for name, s, chrom in TESTS:
    new, old = pansoma(NEW, s, chrom, ""), pansoma(OLD, s, chrom, "_normfix")
    print(f"| {name} | {fmt(new)} | {fmt(old)} | {fmt(comparator('deepsomatic', s, chrom))} | {fmt(comparator('clairs_to', s, chrom))} |")
for name, s, chrom in TESTS:
    new, old = pansoma(NEW, s, chrom, ""), pansoma(OLD, s, chrom, "_normfix")
    if new and old: print(f"\nPASS threshold {name}: site features {new[3]:.3f}, normfix {old[3]:.3f}", end="")
print("\n")

def recs(f, chrom):
    out = {}
    for l in gzip.open(f, "rt"):
        if l[0] == "#": continue
        x = l.rstrip("\n").split("\t")
        if x[0] == chrom and len(x[3]) == 1 and len(x[4]) == 1: out[(int(x[1]), x[3], x[4])] = x
    return out

for chrom in ("chr20", "chr1"):
    s = f"HCC1395T_{p}"
    ON, NN = (f"{OLD}/vcf_{chrom}_{s}_d46_normfix/eval_bed", f"{NEW}/vcf_{chrom}_{s}_d46/eval_bed")
    if not os.path.exists(f"{NN}/vcfeval/pon_calls/tp-baseline.vcf.gz"): continue
    calls = recs(f"{ON}/{s}_d46.SNV.linear.pon.vcf.gz", chrom)
    af = {k: float(dict(zip(x[8].split(":"), x[9].split(":")))["AF"]) for k, x in calls.items()}
    truth = set(recs(f"{ON}/vcfeval/pon_calls/tp-baseline.vcf.gz", chrom)) | set(recs(f"{ON}/vcfeval/pon_calls/fn.vcf.gz", chrom))
    found = {"site features": set(recs(f"{NN}/vcfeval/pon_calls/tp-baseline.vcf.gz", chrom)),
             "normfix": set(recs(f"{ON}/vcfeval/pon_calls/tp-baseline.vcf.gz", chrom))}
    for c, n in (("deepsomatic", "DeepSomatic-TO"), ("clairs_to", "ClairS-TO")):
        found[n] = set(recs(f"{E8}/{c}/{s}{'_chr20' if chrom == 'chr20' else ''}/SNV/eval_bed/pass_pon/tp-baseline.vcf.gz", chrom))
    print(f"## {p}: HCC1395T {chrom} recall by read AF (truth SNVs with a Pansoma tensor in the VCF; after the PoN, BED)\n")
    print("| AF | truth | " + " | ".join(found) + " |"); print("|---" * (len(found) + 2) + "|")
    for lo, hi in [(0.08, 0.25), (0.25, 0.4), (0.4, 0.6), (0.6, 0.85), (0.85, 1.01)]:
        t = [k for k in truth if k in af and lo <= af[k] < hi]
        if t: print(f"| [{lo}, {min(hi, 1.0)}] | {len(t)} | " + " | ".join(f"{sum(k in f for k in t) / len(t):.3f}" for f in found.values()) + " |")
    t = [k for k in truth if k not in af]
    if t: print(f"| no tensor in the VCF | {len(t)} | " + " | ".join(f"{sum(k in f for k in t) / len(t):.3f}" for f in found.values()) + " |")
    print()

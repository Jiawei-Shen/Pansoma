# python compare.py [PLAT] -- the evaluation of README.md (fixed before results): SNVs after the PoN in the benchmark BED.
# Models: the AF 0.05 model, the control (floor 0.08, same negative sampling), sitef_ppr9 (floor 0.08, every negative),
# normfix (no site features), DeepSomatic-TO and ClairS-TO (E8). Site-feature models at rule 1 (their PASS: post-PoN
# validation callable recall 0.9) and rule 2 (the same at 0.95, the best epoch's validation threshold, from the
# vcfeval pon_curve of the same VCF). Then HCC1395T recall by read AF (Pansoma's tensor AF; 0.05-0.08 = the new band).
import glob, gzip, json, os, sys
R = "/scratch/jshen/data/pansoma_net_v2_runs"
E8 = "/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score"
# usage: python compare.py [PLAT [AF05_RUN CONTROL_RUN]]   (PLAT PacBio (default) | Illumina | ONT)
PLAT = sys.argv[1] if len(sys.argv) > 1 else "PacBio"
AF05, CTRL = sys.argv[2:4] if len(sys.argv) > 3 else (f"multicancer6_{PLAT}_SNV_sitef_af05_pnf25", f"multicancer6_{PLAT}_SNV_sitef_pnf25")
SITEF = [("AF 0.05 + sampling", AF05), ("0.08 + sampling (control)", CTRL), ("0.08, all negatives (ppr9)", f"multicancer6_{PLAT}_SNV_sitef_ppr9")]
NORMFIX = f"multicancer6_{PLAT}_SNV"
CO = "COLO829T_fiberseq" if PLAT == "PacBio" else f"COLO829T_{PLAT}"
TESTS = [("HCC1395T chr20", f"HCC1395T_{PLAT}", "chr20"), ("HCC1395T chr1", f"HCC1395T_{PLAT}", "chr1"),
         (f"HG008T chr1", f"HG008T_{PLAT}", "chr1"), ("COLO829T chr1", CO, "chr1")]
ps = lambda info: float(info.split("P_SOMATIC=")[1].split(";")[0])

def recs(f, chrom):
    out = {}
    for l in gzip.open(f, "rt"):
        if l[0] == "#": continue
        x = l.rstrip("\n").split("\t")
        if x[0] == chrom and len(x[3]) == 1 and len(x[4]) == 1: out[(int(x[1]), x[3], x[4])] = x
    return out

def cutoffs(run):
    L = [json.loads(l) for l in open(f"{R}/{run}/metrics.jsonl")]
    best = max(L, key=lambda r: r["val"]["postpon"]["f1"])
    assert f"best epoch {best['epoch']} " in open(f"{R}/{run}/train.log").read(), run
    v = best["val"]["somatic_at_recall_postpon_callable"]
    return best["epoch"], v["0.9"]["threshold"], v["0.95"]["threshold"]

def at(run, s, chrom, cut):
    """P, R, F1 and the true calls (keys) at p_somatic >= cut, or at PASS (cut None, report.json)."""
    E = f"{R}/{run}/vcf_{chrom}_{s}_d46{'_normfix' if run == NORMFIX else ''}/eval_bed"
    if not os.path.exists(f"{E}/report.json"): return None
    a = json.load(open(f"{E}/report.json"))["pon"]["at_threshold"]; n = a["tp_baseline"] + a["fn"]
    if cut is None:
        return a["precision"], a["recall"], a["f1"], set(recs(f"{E}/vcfeval/pon_calls/tp-baseline.vcf.gz", chrom))
    tp = {k: ps(x[7]) for k, x in recs(f"{E}/vcfeval/pon_curve/tp.vcf.gz", chrom).items()}
    fp = [ps(x[7]) for x in recs(f"{E}/vcfeval/pon_curve/fp.vcf.gz", chrom).values()]
    T, F = sum(v >= cut for v in tp.values()), sum(v >= cut for v in fp)
    P, Rc = T / max(1, T + F), T / n
    return P, Rc, 2 * P * Rc / max(1e-12, P + Rc), {k for k, v in tp.items() if v >= cut}

def comparator(c, s, chrom):
    d = f"{E8}/{c}/{s}{'_chr20' if chrom == 'chr20' else ''}/SNV"
    if not os.path.exists(f"{d}/report.json"): return None
    a = json.load(open(f"{d}/report.json"))["eval_bed"]["pass_pon"]["at_threshold"]
    return a["precision"], a["recall"], a["f1"], set(recs(f"{d}/eval_bed/pass_pon/tp-baseline.vcf.gz", chrom))

cols = []
for name, run in SITEF:
    try: ep, c1, c2 = cutoffs(run)
    except FileNotFoundError: continue
    print(f"{name} ({run}): best epoch {ep}, rule 1 p >= {c1:.3f}, rule 2 p >= {c2:.3f}")
    cols += [(f"{name}, rule 1", lambda s, ch, run=run: at(run, s, ch, None)),
             (f"{name}, rule 2", lambda s, ch, run=run, c2=c2: at(run, s, ch, c2))]
cols += [("normfix", lambda s, ch: at(NORMFIX, s, ch, None)),
         ("DeepSomatic-TO", lambda s, ch: comparator("deepsomatic", s, ch)), ("ClairS-TO", lambda s, ch: comparator("clairs_to", s, ch))]
fmt = lambda r: "-" if r is None else f"{r[2]:.3f} ({r[0]:.3f} / {r[1]:.3f})"
res = {(t, c): f(s, ch) for t, s, ch in TESTS for c, f in cols}
print(f"\n## {PLAT} SNV F1 (P / R) after the PoN, benchmark BED\n")
print("| model | " + " | ".join(t for t, _, _ in TESTS) + " |"); print("|---" * (len(TESTS) + 1) + "|")
for c, _ in cols: print(f"| {c} | " + " | ".join(fmt(res[(t, c)]) for t, _, _ in TESTS) + " |")

BINS = [(0.05, 0.08), (0.08, 0.25), (0.25, 0.4), (0.4, 0.6), (0.6, 0.85), (0.85, 1.01)]
for t, s, chrom in TESTS[:2]:
    base = f"{R}/{AF05}/vcf_{chrom}_{s}_d46/eval_bed"; old = f"{R}/{NORMFIX}/vcf_{chrom}_{s}_d46_normfix/eval_bed"
    if not os.path.exists(base): continue
    af = {}
    for f in [f"{old}/{s}_d46.SNV.linear.pon.vcf.gz"] + glob.glob(f"{base}/*.SNV.linear.pon.vcf.gz"):
        for k, x in recs(f, chrom).items(): af[k] = float(dict(zip(x[8].split(":"), x[9].split(":")))["AF"])
    truth = set(recs(f"{base}/vcfeval/pon_calls/tp-baseline.vcf.gz", chrom)) | set(recs(f"{base}/vcfeval/pon_calls/fn.vcf.gz", chrom))
    print(f"\n## {t} recall by read AF (Pansoma tensor AF; truth {len(truth):,})\n")
    print("| AF | truth | " + " | ".join(c for c, _ in cols) + " |"); print("|---" * (len(cols) + 2) + "|")
    for lo, hi in BINS + [(None, None)]:
        tt = [k for k in truth if k not in af] if lo is None else [k for k in truth if k in af and lo <= af[k] < hi]
        if not tt: continue
        lab = "no tensor" if lo is None else f"[{lo}, {min(hi, 1.0)}]"
        print(f"| {lab} | {len(tt)} | " + " | ".join("-" if res[(t, c)] is None else f"{sum(k in res[(t, c)][3] for k in tt) / len(tt):.3f}" for c, _ in cols) + " |")

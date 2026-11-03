import gzip, json, sys
from collections import Counter
p = sys.argv[1] if len(sys.argv) > 1 else "Illumina"
run = sys.argv[2] if len(sys.argv) > 2 else "multicancer6"
R = f"/scratch/jshen/data/pansoma_net_v2_runs/{run}_{p}_SNV/vcf_chr1_HCC1395T_{p}_d46_normfix/eval_bed"
thr = json.load(open(f"{R}/report.json"))["threshold"]
def recs(f):
    out = {}
    for l in gzip.open(f, "rt"):
        if l[0] == "#": continue
        x = l.rstrip("\n").split("\t")
        if len(x[3]) == 1 and len(x[4]) == 1: out[(int(x[1]), x[3], x[4])] = x
    return out
calls = recs(f"{R}/HCC1395T_{p}_d46.SNV.linear.pon.vcf.gz")
tpb = recs(f"{R}/vcfeval/pon_calls/tp-baseline.vcf.gz"); fn = recs(f"{R}/vcfeval/pon_calls/fn.vcf.gz")
st = {}
for l in open(f"/scratch/jshen/data/pansoma_v2_tensors/HCC1395T_{p}_d46/tensors/somatic.recall.tsv"):
    x = l.rstrip("\n").split("\t")
    if x[1] == "chr1" and x[3] == "SNP": st[int(x[2])] = x[6] + (":" + x[7] if x[6] == "filtered" else "")
def info(x):
    d = dict(kv.split("=", 1) for kv in x[7].split(";") if "=" in kv)
    fmt = dict(zip(x[8].split(":"), x[9].split(":")))
    return float(d["P_SOMATIC"]), float(d["P_GERMLINE"]), float(d["P_NON"]), float(fmt["AF"]), int(fmt["DP"]), x[6]
cls = Counter(); rows = []
for k in fn:
    s = st.get(k[0], "?")
    if s != "tensor_representative":
        cls["no tensor: " + s] += 1; continue
    c = calls.get(k)
    if c is None:
        cls["tensor, P_SOMATIC < 0.01 (not in VCF)"] += 1; continue
    ps, pg, pn, af, dp, flt = info(c)
    if "PanelOfNormals" in flt:
        cls["tensor, population allele (PoN)"] += 1; continue
    top = max((ps, "somatic"), (pg, "germline"), (pn, "non"))[1]
    cls[f"tensor, below threshold ({top} is the top class)"] += 1
    rows.append((ps, pg, pn, af, dp))
print(p, run, "truth", len(tpb) + len(fn), "TP", len(tpb), "FN", len(fn), "threshold", round(thr, 3))
for k, v in cls.most_common(): print(f"  {k}: {v}")
import statistics as S
tp_af = [info(calls[k])[3] for k in tpb if k in calls]; tp_dp = [info(calls[k])[4] for k in tpb if k in calls]
def q(v): v = sorted(v); return f"median {S.median(v):.2f} [IQR {v[len(v)//4]:.2f}-{v[3*len(v)//4]:.2f}]"
print("  TP AF", q(tp_af), "DP", q(tp_dp))
print("  FN below thr AF", q([r[3] for r in rows]), "DP", q([r[4] for r in rows]), "P_SOMATIC", q([r[0] for r in rows]))
bins = [(0, 0.15), (0.15, 0.25), (0.25, 0.4), (0.4, 0.6), (0.6, 0.85), (0.85, 1.01)]
for lo, hi in bins:
    a = sum(lo <= x < hi for x in tp_af); b = [r for r in rows if lo <= r[3] < hi]
    g = sum(r[1] > r[0] for r in b)
    print(f"  AF [{lo},{hi}): TP {a}, FN below thr {len(b)} (germline > somatic in {g}); recall of tensor truth {a/(a+len(b)) if a+len(b) else 0:.2f}")

"""v2 chr1 SNV per-tensor numbers on v1's definition (read-only inputs).

v1 definition (pasted log, best epoch 67 = published HG008_WGS_SNV.pth): chr1 val = first 8 v1 shards (262,144 tensors,
524 true); every linear v1 candidate with AF >= 0.08, inside or outside any BED; true = PASS somatic SNP truth allele,
false = everything else (germline included). Here: same rule applied to v2 tensors.
"""
import gzip
import json
import pickle
import sys
from collections import Counter

import numpy as np

HERE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"
sys.path.insert(0, HERE)
from common import V2, TRUTH_TSV, confident_chr1  # noqa

RUNS = "/scratch/jshen/data/pansoma_net_v2_runs"
NAME = "Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz"
conf = confident_chr1()

# PASS chr1 SNP somatic truth alleles, any BED
pass_snp = {}
with open(TRUTH_TSV) as f:
    hdr = next(f).rstrip("\n").split("\t")
    for line in f:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        if r["chrom"] == "chr1" and r["kind"] == "SNP" and r["passed"] == "True":
            pass_snp[(int(r["pos0"]), r["alt"])] = int(r["truth_id"])
print("PASS chr1 SNP truth alleles (any BED)", len(pass_snp))

# v1 range of the pasted log: the first 8 shards (rows < 262144) of the v1 chr1 val set
recs = pickle.load(open(f"{HERE}/v1_val_recs.pkl", "rb"))
v1_last_pos = max(r[0] for r in recs[:262144])
print("v1 log range: chr1 pos0 <=", v1_last_pos, "true in first 8 shards", sum(r[2] for r in recs[:262144]))

af = []
with open(f"{V2}/SNV/chr1_variant_summary.ndjson") as f:
    for line in f:
        i = line.find('"af": ')
        af.append(float(line[i + 6:line.index(",", i)]))
af = np.array(af)
lin = []
with open(f"{V2}/SNV/chr1_labels.ndjson") as f:
    for line in f:
        r = json.loads(line)
        g = r.get("grch38")
        lin.append((g["pos0"], g["alt"]) if g else None)
print("tensors", len(lin), "af parsed", len(af))

pos = np.array([k[0] if k else -1 for k in lin])
is_lin = pos >= 0
true = np.array([k is not None and k in pass_snp for k in lin])
inconf = np.array([k is not None and conf.contains0(k[0]) for k in lin])
print("v2 linear tensors", is_lin.sum(), "AF>=0.08", (is_lin & (af >= 0.08)).sum(), "true", true.sum())


def pr_at_recall(p, y, r):
    o = np.argsort(-p, kind="stable")
    ys = y[o]
    tp = np.cumsum(ys)
    k = np.searchsorted(tp, np.ceil(r * y.sum()))
    return tp[k] / (k + 1), int(tp[k]), int(k + 1 - tp[k]), float(p[o][k])


out = {}
for run in sys.argv[1:]:
    ps, pred = [], []
    with gzip.open(f"{RUNS}/{run}/test_chr1/{NAME}", "rt") as f:
        for line in f:
            p = json.loads(line)
            ps.append(p["p_somatic"])
            pred.append(p["pred"] == "somatic")
    ps, pred = np.array(ps), np.array(pred)
    np.save(f"{HERE}/v2_{run}_psom.npy", ps)
    res = {}
    for name, m in (("v1_def_range262k", is_lin & (af >= 0.08) & (pos <= v1_last_pos)),
                    ("v1_def_all_chr1", is_lin & (af >= 0.08)),
                    ("v1_def_all_chr1_inconf", is_lin & (af >= 0.08) & inconf)):
        y = true[m]
        p = ps[m]
        tp = int((pred[m] & y).sum()); fp = int((pred[m] & ~y).sum())
        res[name] = dict(tensors=int(m.sum()), true=int(y.sum()),
                         argmax=dict(tp=tp, fp=fp, P=round(tp / max(tp + fp, 1), 4), R=round(tp / y.sum(), 4)),
                         P_at_R0935=[round(float(v), 4) if isinstance(v, float) else v for v in pr_at_recall(p, y, 0.935)],
                         P_at_R08=[round(float(v), 4) if isinstance(v, float) else v for v in pr_at_recall(p, y, 0.8)])
    out[run] = res
    print(run, json.dumps(res))
json.dump(out, open(f"{HERE}/v2_tensor_level_v1def.json", "w"), indent=1)
np.save(f"{HERE}/v2_af.npy", af)
pickle.dump(dict(lin=lin, true=true, inconf=inconf), open(f"{HERE}/v2_lin.pkl", "wb"))

"""Composition of v2 false calls (label 0 vs 2, AF bins) at fixed truth recalls (read-only)."""
import gzip
import json
import pickle
import sys
from collections import Counter

import numpy as np

HERE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"
sys.path.insert(0, HERE)
from common import truth_697  # noqa

truth = truth_697()
v2 = pickle.load(open(f"{HERE}/v2_lin.pkl", "rb"))
lin = v2["lin"]
af = np.load(f"{HERE}/v2_af.npy")
tkey = {(t["pos0"], t["alt"]) for t in truth.values()}
for run in sys.argv[1:]:
    ps = np.load(f"{HERE}/v2_{run}_psom.npy")
    tl, intest, tids = [], [], []
    with gzip.open(f"/scratch/jshen/data/pansoma_net_v2_runs/{run}/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz", "rt") as f:
        for line in f:
            p = json.loads(line)
            tl.append(-9 if p["test_label"] is None else p["test_label"]); intest.append(p["in_test"])
            tids.append([t for t in p.get("truth_ids") or [] if t in truth])
    tl, intest = np.array(tl), np.array(intest)
    best = {}
    for i in np.nonzero(intest)[0]:
        for t in tids[i]:
            best[t] = max(best.get(t, 0), ps[i])
    tp_s = np.sort(np.array(list(best.values())))[::-1]
    fpmask = intest & np.isin(tl, [0, 2]) & np.array([k is not None and k not in tkey for k in lin])
    print(run, "false in-test GRCh38 tensors", fpmask.sum(), "of which germline", (fpmask & (tl == 2)).sum())
    for r in (0.5, 0.8, 0.9):
        thr = tp_s[int(np.ceil(r * len(truth))) - 1]
        called = fpmask & (ps >= thr)
        g = called & (tl == 2)
        n = called & (tl == 0)
        bins = [0.07, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 1.01]
        hg = np.histogram(af[g], bins)[0]; hn = np.histogram(af[n], bins)[0]
        tpafs = []
        print(f"  R={r} thr={thr:.4f} FP={called.sum()} germline={g.sum()} ({g.sum() / max(called.sum(), 1):.2%}) non={n.sum()}")
        print("    AF bins", bins[:-1], "\n    germline", hg.tolist(), "\n    non     ", hn.tolist())
    # germline tensors overall: how many get p_somatic > 0.5
    gm = intest & (tl == 2)
    print("  germline test tensors", gm.sum(), "p_som>=0.5:", (ps[gm] >= 0.5).sum(), " AF median", float(np.median(af[gm])))
    sm = np.array([bool(x) for x in tids]) & intest
    print("  truth tensors AF median", float(np.median(af[sm])), "AF quantiles", np.quantile(af[sm], [0.1, 0.25, 0.5, 0.75, 0.9]).round(3).tolist())
    print("  germline tensors AF quantiles", np.quantile(af[gm], [0.1, 0.25, 0.5, 0.75, 0.9]).round(3).tolist())

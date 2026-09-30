"""Bootstrap CI of the v1 sample-estimated truth-level points (resample truths and each FP stratum)."""
import glob
import pickle
import sys

import numpy as np

HERE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"
sys.path.insert(0, HERE)
from common import truth_697  # noqa

truth = truth_697()
tkey = {(t["pos0"], t["alt"]) for t in truth.values()}
recs = pickle.load(open(f"{HERE}/v1_val_recs.pkl", "rb"))
v2 = pickle.load(open(f"{HERE}/v2_lin.pkl", "rb"))
import json
labels = np.array([json.loads(l)["label"] for l in open("/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors/SNV/chr1_labels.ndjson")])
som_sites = {v2["lin"][i] for i in np.nonzero(labels == 1)[0] if v2["lin"][i] is not None}
rows, strata, pt = [], [], []
for p in sorted(glob.glob(f"{HERE}/v1_sample_part*.npz")):
    z = np.load(p)
    rows.append(z["rows"]); strata.append(z["strata"]); pt.append(z["p_true"])
    nG_all, nN_all = int(z["n_G_all"]), int(z["n_N_all"])
rows, strata, pt = np.concatenate(rows), np.concatenate(strata), np.concatenate(pt)
tp = np.array([p for r, s, p in zip(rows, strata, pt) if s == "T" and (recs[r][0], recs[r][1]) in tkey])
keep = np.array([s != "T" and (recs[r][0], recs[r][1]) not in som_sites and (recs[r][0], recs[r][1]) not in tkey
                 for r, s in zip(rows, strata)])
G = pt[keep & (strata == "G")]
Nn = pt[keep & (strata == "N")]
wG, wN = nG_all / (strata == "G").sum(), nN_all / (strata == "N").sum()
N = len(truth)


def points(tp, G, Nn):
    thr = np.unique(np.concatenate([tp, G, Nn]))[::-1]
    t = len(tp) - np.searchsorted(np.sort(tp), thr)
    f = wG * (len(G) - np.searchsorted(np.sort(G), thr)) + wN * (len(Nn) - np.searchsorted(np.sort(Nn), thr))
    rec, prec = t / N, t / np.maximum(t + f, 1e-9)
    f1 = 2 * prec * rec / np.maximum(prec + rec, 1e-12)
    out = [f1.max()]
    for r in (0.5, 0.8, 0.9):
        m = rec >= r
        out.append(prec[m].max() if m.any() else np.nan)
    m = prec >= 0.10
    out.append(rec[m].max() if m.any() else 0)
    return out


est = points(tp, G, Nn)
rng = np.random.default_rng(1)
bs = np.array([points(rng.choice(tp, len(tp)), rng.choice(G, len(G)), rng.choice(Nn, len(Nn))) for _ in range(300)])
names = ["maxF1", "P@R0.5", "P@R0.8", "P@R0.9", "R@P0.10"]
for n, e, lo, hi in zip(names, est, np.nanpercentile(bs, 2.5, 0), np.nanpercentile(bs, 97.5, 0)):
    print(f"v1 {n}: {e:.4f}  95% CI [{lo:.4f}, {hi:.4f}]")
print("sampled G", len(G), "N", len(Nn), "w", round(wG, 2), round(wN, 2),
      "G p>=0.5:", int((G >= 0.5).sum()), "N p>=0.5:", int((Nn >= 0.5).sum()),
      "est. FP at p>=0.5:", round(wG * (G >= 0.5).sum() + wN * (Nn >= 0.5).sum()), "TP at p>=0.5:", int((tp >= 0.5).sum()))

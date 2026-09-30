"""Head-to-head truth-level curves: v1 (stratified-sample estimate) vs v2 (exact), same 697 truths, same region.

FP units = distinct (pos0, alt) sites inside somatic BED ∩ germline BED; tensors whose site is a v2 label-1 tensor
(somatic, incl. partial INDEL matches) are neutral for both models, as in the v2 metrics.
"""
import glob
import json
import pickle
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"
sys.path.insert(0, HERE)
from common import V2, truth_697, confident_chr1  # noqa
from curve import summary  # noqa

truth = truth_697()
N = len(truth)
tkey = {(t["pos0"], t["alt"]): tid for tid, t in truth.items()}

# ---------------- v2 (exact) ----------------
v2 = pickle.load(open(f"{HERE}/v2_lin.pkl", "rb"))
lin = v2["lin"]
af2 = np.load(f"{HERE}/v2_af.npy")
labels = []
with open(f"{V2}/SNV/chr1_labels.ndjson") as f:
    for line in f:
        r = json.loads(line)
        labels.append(r["label"])
labels = np.array(labels)
v2_som_sites = {lin[i] for i in np.nonzero(labels == 1)[0] if lin[i] is not None}
v2_site_row = {}
for i, k in enumerate(lin):
    if k is not None:
        v2_site_row.setdefault(k, i)


def weighted_curve(tp_scores, fp_scores, fp_weights, n_truth):
    tp_scores = np.asarray(tp_scores, float)
    fp_scores = np.asarray(fp_scores, float)
    fp_weights = np.asarray(fp_weights, float)
    thr = np.unique(np.concatenate([tp_scores, fp_scores]))[::-1]
    tps = np.sort(tp_scores)
    tp = len(tps) - np.searchsorted(tps, thr, side="left")
    o = np.argsort(fp_scores)
    fs, fw = fp_scores[o], fp_weights[o]
    cw = np.concatenate([[0], np.cumsum(fw)])
    fp = cw[-1] - cw[np.searchsorted(fs, thr, side="left")]
    rec = tp / n_truth
    prec = np.where(tp + fp > 0, tp / np.maximum(tp + fp, 1e-9), 1.0)
    f1 = np.where(prec + rec > 0, 2 * prec * rec / np.maximum(prec + rec, 1e-12), 0)
    return dict(thr=thr, tp=tp, fp=np.round(fp).astype(int), rec=rec, prec=prec, f1=f1)


def v2_units(run, af_min=0.0):
    import gzip
    ps = np.load(f"{HERE}/v2_{run}_psom.npy")
    tp = defaultdict(float)
    fp = defaultdict(float)
    with gzip.open(f"/scratch/jshen/data/pansoma_net_v2_runs/{run}/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz", "rt") as f:
        for i, line in enumerate(f):
            if not line.startswith('{"chrom"'):
                continue
            if af2[i] < af_min:
                continue
            p = json.loads(line)
            if not p["in_test"]:
                continue
            for t in p.get("truth_ids") or []:
                if t in truth:
                    tp[t] = max(tp[t], ps[i])
            if p["test_label"] in (0, 2) and lin[i] is not None and lin[i] not in tkey:
                fp[lin[i]] = max(fp[lin[i]], ps[i])
    return tp, fp


# ---------------- v1 (sample estimate) ----------------
recs = pickle.load(open(f"{HERE}/v1_val_recs.pkl", "rb"))  # (pos0, alt, lab, inconf, isgerm, sidx, iws, af)
if "--v1-full" in sys.argv:  # exact run (GPU output of v1_full_chr1_gpu.py): every in-confident row, weight 1
    j = sys.argv.index("--v1-full")
    full = np.load(sys.argv[j + 1])
    del sys.argv[j:j + 2]
    rows = np.array([i for i, r in enumerate(recs) if r[2] == 1 or r[3]])
    strata = np.array(["T" if recs[i][2] == 1 else ("G" if recs[i][4] else "N") for i in rows])
    pt = full[rows]
    nG_all, nN_all = int((strata == "G").sum()), int((strata == "N").sum())
    parts = ["full"]
else:
    parts = sorted(glob.glob(f"{HERE}/v1_sample_part*.npz"))
    rows, strata, pt = [], [], []
    for p in parts:
        z = np.load(p)
        rows.append(z["rows"]); strata.append(z["strata"]); pt.append(z["p_true"])
        nG_all, nN_all, nG, nN = int(z["n_G_all"]), int(z["n_N_all"]), int(z["n_G"]), int(z["n_N"])
    rows, strata, pt = np.concatenate(rows), np.concatenate(strata), np.concatenate(pt)
nparts = len(parts)
sG = (strata == "G").sum(); sN = (strata == "N").sum()
wG, wN = nG_all / sG, nN_all / sN
print(f"v1 sample parts={nparts} T={int((strata == 'T').sum())} G={sG} (w={wG:.2f}) N={sN} (w={wN:.2f})")
v1_tp, v1_fp_s, v1_fp_w, v1_fp_lab, v1_neutral = {}, [], [], [], 0
pair_tp, pair_fp = [], []
for r, s, p in zip(rows, strata, pt):
    pos0, alt, lab, inconf, isgerm, *_ = recs[r]
    k = (pos0, alt)
    if s == "T":
        if k in tkey:
            v1_tp[tkey[k]] = max(v1_tp.get(tkey[k], 0.0), float(p))
        continue
    if k in v2_som_sites or k in tkey:
        v1_neutral += 1
        continue
    v1_fp_s.append(float(p)); v1_fp_w.append(wG if s == "G" else wN); v1_fp_lab.append(s)
    j = v2_site_row.get(k)
    pair_fp.append((s, float(p), j))
print("v1 truths with tensor scored", len(v1_tp), "neutral sampled", v1_neutral)

results = {}
c1 = weighted_curve(list(v1_tp.values()), v1_fp_s, v1_fp_w, N)
results["v1_published_HG008_WGS_SNV (sample est.)"] = summary(c1, N, len(v1_tp))
# v1, germline-stratum FPs only / non-stratum only
m = np.array(v1_fp_lab) == "N"
results["v1 FP=non-germline sites only"] = summary(weighted_curve(list(v1_tp.values()), np.array(v1_fp_s)[m], np.array(v1_fp_w)[m], N), N, len(v1_tp))
results["v1 FP=germline sites only"] = summary(weighted_curve(list(v1_tp.values()), np.array(v1_fp_s)[~m], np.array(v1_fp_w)[~m], N), N, len(v1_tp))
for run in sys.argv[1:]:
    for af_min in (0.0, 0.08):
        tp, fp = v2_units(run, af_min)
        c = weighted_curve(list(tp.values()), list(fp.values()), np.ones(len(fp)), N)
        results[f"v2 {run} AF>={af_min}"] = summary(c, N, len(tp))
    if run == sys.argv[1]:
        # v2 scores on exactly the v1-sampled sites (same weights): isolates the model from candidate-set size
        ps = np.load(f"{HERE}/v2_{run}_psom.npy")
        tp2 = {}
        for tid, t in truth.items():
            j = v2_site_row.get((t["pos0"], t["alt"]))
            if j is not None and tid in v1_tp:
                tp2[tid] = ps[j]
        fs, fw, missing = [], [], Counter()
        for (s, p, j), w in zip(pair_fp, v1_fp_w):
            if j is None:
                missing[s] += 1
                fs.append(0.0)  # site has no v2 tensor: v2 never calls it
            else:
                fs.append(ps[j])
            fw.append(w)
        results[f"v2 {run} scored on v1-sampled sites"] = summary(weighted_curve(list(tp2.values()), fs, fw, N), N, len(tp2))
        results[f"v2 {run} scored on v1-sampled sites"]["v1_sites_without_v2_tensor"] = dict(missing)
        # per-truth paired comparison
        both = [(v1_tp[t], tp2[t]) for t in tp2]
        a = np.array(both)
        results["paired_truth_scores"] = dict(n=len(a), v1_median=float(np.median(a[:, 0])), v2_median=float(np.median(a[:, 1])),
                                              v1_ge_0p5=int((a[:, 0] >= 0.5).sum()), v2_ge_0p5=int((a[:, 1] >= 0.5).sum()))
for k, v in results.items():
    print(k, json.dumps(v))
json.dump(results, open(f"{HERE}/head_to_head_results.json", "w"), indent=1)

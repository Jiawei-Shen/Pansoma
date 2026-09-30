"""Bootstrap 95% CI of the v1 e067 sample-estimated truth-level points (resample truths T and each FP stratum G, N),
no PoN and with both PoN rules. Same unit definitions as h2h_final.py. Read-only."""
import glob, json, pickle, sys
import numpy as np
HERE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"
sys.path.insert(0, HERE)
from common import V2, truth_697  # noqa
from pon_annot import annotate, rule_repo, rule_af001  # noqa
import pysam
truth = truth_697(); NT = len(truth)
tkey = {(t["pos0"], t["alt"]): tid for tid, t in truth.items()}
neutral = set()
for line in open(f"{V2}/SNV/chr1_labels.ndjson"):
    r = json.loads(line)
    if r["label"] == 1 and r.get("grch38") and (r["grch38"]["pos0"], r["grch38"]["alt"]) not in tkey:
        neutral.add((r["grch38"]["pos0"], r["grch38"]["alt"]))
fa = pysam.FastaFile(json.load(open(f"{V2}/SNV/labels.manifest.json"))["fasta"])
recs = pickle.load(open(f"{HERE}/v1_val_recs.pkl", "rb"))
R, S, P = [], [], []
for pth in sorted(glob.glob(f"{HERE}/v1_sample_part*.npz")):
    z = np.load(pth); R.append(z["rows"]); S.append(z["strata"]); P.append(z["p_true"]); nG, nN = int(z["n_G_all"]), int(z["n_N_all"])
R, S, P = np.concatenate(R), np.concatenate(S), np.concatenate(P)
wG, wN = nG / (S == "G").sum(), nN / (S == "N").sum()
sites = [(recs[r][0], fa.fetch("chr1", recs[r][0], recs[r][0] + 1).upper(), recs[r][1]) for r in R]
ann = annotate(sites)
print(f"sample T={(S == 'T').sum()} G={(S == 'G').sum()} (w {wG:.2f}) N={(S == 'N').sum()} (w {wN:.2f})")


def points(tp, G, Nn):
    thr = np.unique(np.concatenate([tp, G, Nn]))[::-1]
    t = len(tp) - np.searchsorted(np.sort(tp), thr)
    f = wG * (len(G) - np.searchsorted(np.sort(G), thr)) + wN * (len(Nn) - np.searchsorted(np.sort(Nn), thr))
    rec, prec = t / NT, t / np.maximum(t + f, 1e-9)
    f1 = 2 * prec * rec / np.maximum(prec + rec, 1e-12)
    o = [f1.max()]
    for x in (0.05, 0.10, 0.15, 0.20):
        m = prec >= x; o.append(rec[m].max() if m.any() else 0.0)
    for x in (0.5, 0.7, 0.8, 0.9):
        m = rec >= x; o.append(prec[m].max() if m.any() else np.nan)
    return o


names = ["maxF1", "R@P0.05", "R@P0.10", "R@P0.15", "R@P0.20", "P@R0.5", "P@R0.7", "P@R0.8", "P@R0.9"]
out = {}
for rn, drop in (("noPoN", lambda a: False), ("PoN_repo", rule_repo), ("PoN_af001", rule_af001)):
    tp, G, Nn = {}, [], []
    for r, s, p, k in zip(R, S, P, sites):
        key = (k[0], k[2])
        if drop(ann[k]):
            continue
        if s == "T":
            if key in tkey:
                tp[tkey[key]] = max(tp.get(tkey[key], 0.0), float(p))
            continue
        if key in tkey or key in neutral:
            continue
        (G if s == "G" else Nn).append(float(p))
    tp, G, Nn = np.array(list(tp.values())), np.array(G), np.array(Nn)
    est = points(tp, G, Nn)
    rng = np.random.default_rng(1)
    bs = np.array([points(rng.choice(tp, len(tp)), rng.choice(G, len(G)), rng.choice(Nn, len(Nn))) for _ in range(400)])
    lo, hi = np.nanpercentile(bs, 2.5, 0), np.nanpercentile(bs, 97.5, 0)
    out[rn] = {n: [round(float(e), 4), round(float(a), 4), round(float(b), 4)] for n, e, a, b in zip(names, est, lo, hi)}
    print(rn, "  ".join(f"{n} {e:.3f} [{a:.3f},{b:.3f}]" for n, e, a, b in zip(names, est, lo, hi)))
json.dump(out, open(f"{HERE}/h2h_boot.json", "w"), indent=1)

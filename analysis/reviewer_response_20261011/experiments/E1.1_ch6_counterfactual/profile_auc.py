"""How well do simple ch6 summaries of the chr1 test tensors separate somatic from germline tensors (no model)?

For every E1.1 run of collect.RUNS whose base mode is done (profile.npz and positions.npy, written by
ch6_counterfactual.py or ch6_counterfactual_sampled.py): the base pass's per-tensor ch6 profile
  site      the site code (median ch6 of column 50 over the rows that have one)
  below     the fraction of the tensor's valid ch6 cells below 100 (not every haplotype)
  distinct  the number of distinct ch6 codes in columns 40-60 (a proxy for node changes near the site)
and the AUC of somatic (label 1) against germline (label 2) test tensors, scoring "more somatic" as a higher site
code, a lower fraction below 100 and fewer distinct codes (AUC 0.5 = no separation; ties count one half).
Writes results/profile_auc.json. usage: python profile_auc.py"""
import json
from pathlib import Path

import numpy as np
import torch

import ch6_counterfactual as cf
from collect import RUNS

DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E1.1_ch6_counterfactual")


def auc(pos_scores, neg_scores):
    """P(score of a random positive > that of a random negative) + 0.5 P(tie)."""
    a = np.concatenate([pos_scores, neg_scores])
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a))
    s = a[order]
    i = 0
    while i < len(s):  # average ranks of ties
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    n1, n0 = len(pos_scores), len(neg_scores)
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main():
    out = {}
    for run, set_name, _ in RUNS:
        d = next((x for x in (DATA / run, DATA / "linear_sampled" / run) if (x / "profile.npz").exists()), None)
        if d is None or not (d / "positions.npy").exists():
            continue
        rep = dict(run=run, set=set_name)
        ck = torch.load(cf.RUNS / run / "best.pth", map_location="cpu", weights_only=False)
        index = cf.open_index(cf.TENSORS / set_name / "tensors" / "SNV", DATA / "index_cache")
        a, reasons = index.arrays, index.meta["reasons"]
        ignore = set(ck.get("args", {}).get("ignore_reasons") or [])
        lab = a["eval_label"].astype(np.int64).copy()
        lab[np.isin(a["reason"], [k for k, r in enumerate(reasons) if r in ignore])] = -1
        pos = np.load(d / "positions.npy")
        chr1 = index.meta["chroms"].index("chr1")
        test = (lab[pos] >= 0) & (a["chrom"][pos] == chr1)
        z = np.load(d / "profile.npz")
        som, ger = test & (lab[pos] == 1), test & (lab[pos] == 2)
        r = dict(set=rep["set"], somatic=int(som.sum()), germline=int(ger.sum()))
        for name, sign in (("site", 1), ("below", -1), ("distinct", -1)):
            v = sign * z[name].astype(np.float64)
            r[f"auc_{name}"] = auc(v[som], v[ger])
            r[f"median_{name}_somatic"] = float(np.median(z[name][som]))
            r[f"median_{name}_germline"] = float(np.median(z[name][ger]))
        out[rep["run"]] = r
        print(rep["run"], json.dumps(r))
    res = Path(__file__).resolve().parent / "results"
    res.mkdir(exist_ok=True)
    (res / "profile_auc.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()

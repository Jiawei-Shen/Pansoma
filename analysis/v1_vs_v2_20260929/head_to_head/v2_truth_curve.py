"""v2 chr1 SNV truth-level curve against the 697 PASS in-somatic-BED SNP truths (read-only)."""
import gzip
import json
import sys
from collections import defaultdict, Counter

import numpy as np

sys.path.insert(0, "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head")
from common import V2, truth_697, confident_chr1  # noqa
from curve import curve, summary  # noqa

RUNS = "/scratch/jshen/data/pansoma_net_v2_runs"
NAME = "Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz"

truth = truth_697()
N = len(truth)
conf = confident_chr1()
truth_key = {(t["pos0"], t["alt"]): tid for tid, t in truth.items()}
print("truth", N, "in confident (som∩germ BED):", sum(conf.contains0(t["pos0"]) for t in truth.values()))

labels = []
with open(f"{V2}/SNV/chr1_labels.ndjson") as f:
    for line in f:
        r = json.loads(line)
        g = r.get("grch38")
        labels.append((r["candidate_id"], r["label"], r["reason"], (g["pos0"], g["alt"], g["ref"]) if g else None))
print("label records", len(labels))

results = {}
for run in sys.argv[1:]:
    path = f"{RUNS}/{run}/test_chr1/{NAME}"
    tp_any = defaultdict(float)        # truth -> max p over in_test tensors
    tp_all = defaultdict(float)        # truth -> max p over all tensors (incl. not in test)
    fp_site = defaultdict(float)       # (pos0, alt) -> max p, label 0/2, in test, on GRCh38
    fp_off = defaultdict(float)        # candidate_id -> p, off-reference false, in test
    fp_site_lab = {}
    cnt = Counter()
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            p = json.loads(line)
            cid, lab, reason, g = labels[i]
            assert p["candidate_id"] == cid, (i, p["candidate_id"], cid)
            s = p["p_somatic"]
            tids = [t for t in p.get("truth_ids") or [] if t in truth]
            for t in tids:
                tp_all[t] = max(tp_all[t], s)
                if p["in_test"]:
                    tp_any[t] = max(tp_any[t], s)
            if not p["in_test"]:
                cnt["not_in_test"] += 1
                continue
            tl = p["test_label"]
            cnt[f"test_label_{tl}"] += 1
            if tl in (0, 2):
                if g is None:
                    fp_off[cid] = max(fp_off.get(cid, 0.0), s)
                    cnt["off_ref_false"] += 1
                else:
                    key = (g[0], g[1])
                    if key in truth_key:
                        cnt["false_tensor_on_truth_key"] += 1
                        continue
                    fp_site[key] = max(fp_site[key], s)
                    fp_site_lab[key] = fp_site_lab.get(key, tl) if fp_site_lab.get(key, tl) == tl else "mixed"
            elif tl == 1 and not tids:
                cnt["somatic_label_not_in_697"] += 1
    print(run, dict(cnt), "fp sites", len(fp_site), "off-ref false tensors", len(fp_off),
          "truth with in-test tensor", len(tp_any), "truth with any tensor", len(tp_all))
    lab_counts = Counter(fp_site_lab.values())
    tp_s = list(tp_any.values())
    fp_s = list(fp_site.values())
    c1 = curve(tp_s, fp_s, N)
    c2 = curve(tp_s, fp_s + list(fp_off.values()), N)
    # FP restricted to label-0 sites only (germline sites excluded), for a v1-like "non only" view
    fp_non = [v for k, v in fp_site.items() if fp_site_lab[k] == 0]
    c3 = curve(tp_s, fp_non, N)
    results[run] = dict(
        fp_sites=len(fp_site), fp_site_labels={str(k): v for k, v in lab_counts.items()}, off_ref_false=len(fp_off),
        with_in_test_tensor=len(tp_any), with_any_tensor=len(tp_all),
        grch38_only=summary(c1, N, len(tp_any)),
        with_off_reference=summary(c2, N, len(tp_any)),
        non_only_fp=summary(c3, N, len(tp_any)))
    for k in ("grch38_only", "with_off_reference", "non_only_fp"):
        print(" ", k, json.dumps(results[run][k]))
    np.savez_compressed(f"/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head/v2_{run}_units.npz",
                        tp=np.array(tp_s), fp=np.array(fp_s), fp_off=np.array(list(fp_off.values())))
json.dump(results, open("/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head/v2_truth_curves.json", "w"), indent=1)

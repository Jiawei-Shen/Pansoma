"""E15b (HapMap 200-row SNV model on COLO829T Illumina P100 low-AF rebuild, chr1): tensor-level precision of the
somatic calls by the call's own read AF (top-level 'af' of the tensor's summary record), at the two thresholds of
the predictions' metrics.json. Calls = tensors in the evaluation (in BED, not ignored / excluded, label 0 / 1 / 2)
with p_somatic >= threshold; TP = label 1. Also the same for the tensors with AF < 0.08 only."""
import gzip, json, re, sys
from pathlib import Path
R = Path("/scratch/jshen/data/pansoma_net_v2_runs/HapMap_Illumina_SNV_full/test_chr1_colo829t_illumina")
pred = R / "COLO829T_Illumina_lowaf_chr1.tensors.SNV.predictions.ndjson.gz"
summ = Path("/scratch/jshen/data/pansoma_v2_tensors/COLO829T_Illumina_lowaf_chr1/tensors/SNV/chr1_variant_summary.ndjson")
met = json.load(open(R / "COLO829T_Illumina_lowaf_chr1.tensors.SNV.metrics.json"))
thr = {"pipeline": met["rules"]["pipeline"]["threshold"], "threshold_t": met["rules"]["threshold_t"]["threshold"]}
tmin = min(thr.values())
called = {}
with gzip.open(pred, "rt") as f:
    for line in f:
        r = json.loads(line)
        if r["p_somatic"] < tmin or r["ignored"] or r["excluded"]:
            continue
        lab = r["test_label"] if r["test_label"] is not None else r["label"]
        if lab is None or lab < 0:
            continue
        called[r["candidate_id"]] = (lab, r["p_somatic"])
cid_re, af_re = re.compile(rb'"candidate_id": "([^"]+)"'), re.compile(rb'"af": ([0-9.eE+-]+)')
af = {}
with open(summ, "rb") as f:
    for line in f:
        m = cid_re.search(line, 0, 200)
        c = m.group(1).decode()
        if c in called:
            af[c] = float(af_re.search(line).group(1))
bins = ((0, 0.02), (0.02, 0.05), (0.05, 0.08), (0.08, 0.10), (0.10, 0.20), (0.20, 1.01))
print(f"calls considered {len(called)}, with AF {len(af)}")
print("| threshold | " + " | ".join(f"call AF {a}-{b}" for a, b in bins) + " | all |")
print("|---|" + "---|" * (len(bins) + 1))
for name, t in thr.items():
    cells, TP, N = [], 0, 0
    for a, b in bins:
        tp = n = 0
        for c, (lab, p) in called.items():
            if p >= t and a <= af[c] < b:
                n += 1; tp += lab == 1
        TP += tp; N += n
        cells.append(f"{tp / n:.3f} ({tp}/{n})" if n else "-")
    print(f"| {name} {t:.4f} | " + " | ".join(cells) + f" | {TP / N:.3f} ({TP}/{N}) |")

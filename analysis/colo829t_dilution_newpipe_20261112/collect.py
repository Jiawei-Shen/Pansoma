#!/usr/bin/env python3
"""COLO829T dilutions, SNV chr1: the new-pipeline joint models against ClairS-TO, DeepSomatic-TO and the earlier joint
models. Reads <run>/vcf_chr1_<SET><suffix>/eval_bed/report.json of the best checkpoint of epochs 1-10 (FINAL below:
Illumina and fiberseq = COLO829_<P>_SNV_joint_newpipe_e5, the snapshot after epoch 5, whose best.pth was still the
run's best.pth when the trainings were stopped after epochs 10 and 15; ONT = COLO829_ONT_SNV_joint_newpipe, ended after
epoch 10, its best.pth identical to the epoch-5 snapshot) and the comparator scores of
../indel_ppr8_9919d4a_20261026/comparators/summary.json (the callers' PASS records after our PoN, the same truth, BED
and rtg vcfeval; its "pansoma" SNV entry is COLO829_<P>_SNV_joint_7fc7396_gpu4).
Final PASS thresholds (the user, 2026-11-12): Illumina and fiberseq the validation call threshold - 0.10, written into
the checkpoints (set_call_threshold.sh), so vcf_chr1_<SET> holds those calls and vcf_chr1_<SET>_cr095 the calls at the
recipe's threshold (the post-PoN validation recall 0.95 of the callable positives); ONT the recipe's threshold. The
tests before the change: vcf_chr1_<SET>_tm010 / _tm015 (the recipe's threshold - 0.10 / - 0.15 by vcf.sh THRESHOLD).
Standard library only.
Usage: python3 analysis/colo829t_dilution_newpipe_20261112/collect.py [> tables.md]"""
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
CMP = json.loads((HERE.parent / "indel_ppr8_9919d4a_20261026/comparators/summary.json").read_text())
PLATFORMS = ("Illumina", "fiberseq", "ONT")
PURITIES = (("P100", ""), ("P75", "_P75"), ("P50", "_P50"), ("P25", "_P25"))
FINAL = {"Illumina": "COLO829_Illumina_SNV_joint_newpipe_e5", "fiberseq": "COLO829_fiberseq_SNV_joint_newpipe_e5",
         "ONT": "COLO829_ONT_SNV_joint_newpipe"}  # the run that holds the test of the best of epochs 1-10
RECIPE = {"Illumina": "_cr095", "fiberseq": "_cr095", "ONT": ""}  # the directory suffix of the recipe's threshold
BEST_EPOCH = {"Illumina": 4, "fiberseq": 3, "ONT": 3}  # train.log / snapshot.txt


def report(run, s, suffix=""):
    f = RUNS / run / f"vcf_chr1_{s}{suffix}" / "eval_bed" / "report.json"
    if not f.exists():
        return None
    j = json.loads(f.read_text())
    a = j["pon"]["at_threshold"]
    return dict(p=a["precision"], r=a["recall"], f1=a["f1"], tp=a["tp_baseline"], fp=a["fp"], fn=a["fn"],
                best=j["pon"]["best"]["f1"], ceiling=j["pon"]["ceiling"], threshold=j["threshold"],
                unplaced=j["unplaced_pass"])


def cell(d):
    return "-" if d is None else f"{d['p']:.3f} / {d['r']:.3f} / **{d['f1']:.3f}**"


def counts(d, truth=None):
    fn = d["fn"] if "fn" in d else truth - d["tp"]
    return f"{d['tp']:,} | {d['fp']:,} | {fn:,} | {d['p']:.3f} | {d['r']:.3f} | {d['f1']:.3f}"


rows = []
for p in PLATFORMS:
    for name, suf in PURITIES:
        s = f"COLO829T_{p}{suf}"
        final, recipe = report(FINAL[p], s), report(FINAL[p], s, RECIPE[p])
        variants = {"": final, "_cr095": recipe if RECIPE[p] else None, "_tm010": report(FINAL[p], s, "_tm010"),
                    "_tm015": report(FINAL[p], s, "_tm015")}
        rows.append((p, name, CMP[f"SNV/{s}"], final, recipe, variants))
truth = {(p, name): final["tp"] + final["fn"] for p, name, _, final, _, _ in rows if final}

print("# COLO829T dilutions, SNV, chr1: new pipeline vs ClairS-TO and DeepSomatic-TO\n")
print("(`collect.py`; benchmark BED, after our PoN, rtg vcfeval; cells P / R / **F1** of the PASS calls. One threshold "
      "per Pansoma model for its four purities. Final PASS thresholds (the user, 2026-11-12): Illumina and fiberseq "
      "the validation call threshold - 0.10, ONT the recipe's threshold, the post-PoN validation recall 0.95 of the "
      "callable positives.)\n")
print("| platform | purity | ClairS-TO | DeepSomatic-TO | Pansoma new pipeline, final threshold | Pansoma new pipeline, "
      "the recipe's threshold | Pansoma earlier joint (7fc7396) |")
print("|---|---|---|---|---|---|---|")
for p, name, c, final, recipe, _ in rows:
    print(f"| {p} | {name} | {cell(c['clairs_to'])} | {cell(c['deepsomatic'])} | {cell(final)} | "
          f"{cell(recipe) if RECIPE[p] else 'the final threshold'} | {cell(c['pansoma'])} |")

print("\n## Every caller and threshold (PASS after the PoN; truth 2,443 SNVs on chr1; `final` marks the Pansoma "
      "threshold in use)\n")
print("| platform | purity | caller | PASS threshold | TP | FP | FN | P | R | F1 | TP, FP, F1 against the recipe's "
      "threshold | best F1 of the chr1 curve | recall ceiling |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for p, name, c, final, recipe, variants in rows:
    n = truth.get((p, name))
    for label, d in (("ClairS-TO", c["clairs_to"]), ("DeepSomatic-TO", c["deepsomatic"]), ("earlier joint", c["pansoma"])):
        print(f"| {p} | {name} | {label} | | {counts(d, n)} | | {d['best']:.3f} | {d['ceiling']:.3f} |")
    if RECIPE[p]:
        new = [("new pipeline, the recipe's threshold", recipe, False), ("new pipeline, - 0.10, final", final, True),
               ("new pipeline, - 0.15", variants["_tm015"], True)]
        if final is None:  # before set_call_threshold.sh: the test by vcf.sh THRESHOLD
            new[1] = ("new pipeline, - 0.10", variants["_tm010"], True)
    else:
        new = [("new pipeline, the recipe's threshold, final", final, False)]
    for label, d, lowered in new:
        if d is None:
            continue
        delta = f"{d['tp'] - recipe['tp']:+,}, {d['fp'] - recipe['fp']:+,}, {d['f1'] - recipe['f1']:+.3f}" \
            if lowered and recipe else ""
        print(f"| {p} | {name} | {label} | {d['threshold']:.4f} | {counts(d)} | {delta} | {d['best']:.3f} | "
              f"{d['ceiling']:.3f} |")

print("\n## Models\n")
print("| run | best epoch | final PASS threshold (p_somatic) | the recipe's threshold | PASS calls without a GRCh38 "
      "position at the final threshold (P100 / P75 / P50 / P25) |")
print("|---|---|---|---|---|")
for p in PLATFORMS:
    ds = [(f, r) for q, _, _, f, r, _ in rows if q == p]
    final_t = next((f"{f['threshold']:.4f}" for f, _ in ds if f), "-")
    recipe_t = next((f"{r['threshold']:.4f}" for _, r in ds if r), "-")
    print(f"| `{FINAL[p]}` | {BEST_EPOCH[p]} | {final_t} | {recipe_t} | "
          + " / ".join("-" if f is None else f"{f['unplaced']:,}" for f, _ in ds) + " |")

# the reports of the lowered thresholds by vcf.sh THRESHOLD must equal those of the changed checkpoints
for p, name, _, final, _, variants in rows:
    t = variants["_tm010"]
    if RECIPE[p] and final and t and any(final[k] != t[k] for k in ("tp", "fp", "fn", "threshold")):
        raise SystemExit(f"{p} {name}: vcf_chr1 (checkpoint threshold) differs from _tm010: {final} / {t}")

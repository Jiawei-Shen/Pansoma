#!/usr/bin/env python3
"""E23 at the tensor level (predict's chr1 output, before any VCF): somatic precision / recall / F1 of the chr1 test
tensors at the validation recall-0.9 threshold (the VCF PASS rule), with and without the tensors whose somatic truth
is graph-carried (exclude_graph_carried.py), overall and per VAF range: recall by the truth VAF of the somatic tensors,
precision by the tensor's own allele fraction.
usage: python tensor_metrics_chr1.py --run RUN_DIR --kind SNV|INDEL --exclude exclude.<KIND>.tsv --output OUT.md
"""
import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

import torch

from tensor_vaf import read_vcf_af

SET = Path("/scratch/jshen/data/pansoma_v2_tensors/HapMap_Illumina_nonHG005")
VCF = "/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/HapMap_somatic_snv_indel.vcf.gz"
CACHE = "/scratch/jshen/data/pansoma_net_v2_runs/index_cache"
RANGES = [(0.0099, 0.0199, "1-2 %"), (0.0199, 0.0499, "2-5 %"), (0.0499, 0.0999, "5-10 %"), (0.0999, 0.1651, "10-16.5 %"),
          (0.0099, 0.1651, "1-16.5 %")]


def f1(p, r):
    return 2 * p * r / (p + r) if p + r else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--kind", required=True)
    ap.add_argument("--exclude", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    import sys
    sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning")
    from pansoma_net_v2.data import KindIndex, SCALARS, read_candidate_ids
    run = Path(a.run)
    thr = torch.load(run / "best.pth", map_location="cpu", weights_only=False)["val"]["somatic_at_recall"]["0.9"]["threshold"]
    excluded = read_candidate_ids(a.exclude)
    index = KindIndex(SET / "tensors" / a.kind, CACHE)
    chr1 = index.select(["chr1"], labelled=False)
    cands = index.candidates()
    afi = SCALARS.index("af")
    tensor_af = {cands[int(k)]: float(index.arrays["scalars"][int(k), afi]) for k in chr1}
    del cands
    truth_key = {}
    with open(SET / "truth" / "somatic.graph.tsv") as f:
        col = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
        last = max(col[k] for k in ("truth_id", "chrom", "vcf_pos", "vcf_ref", "vcf_alt")) + 1
        for line in f:
            r = line.split("\t", last)
            if r[col["chrom"]] == "chr1":
                truth_key[int(r[col["truth_id"]])] = ("chr1", int(r[col["vcf_pos"]]), r[col["vcf_ref"]], r[col["vcf_alt"]])
    vaf = read_vcf_af(VCF)
    rows = []  # (excluded, is_somatic, called, truth VAF of a somatic tensor, tensor AF)
    pred = next((run / "test_chr1").glob(f"*.{a.kind}.predictions.ndjson.gz"))
    with gzip.open(pred, "rt") as f:
        for line in f:
            r = json.loads(line)
            if not r["in_test"] or r["test_label"] is None:
                continue
            som = r["test_label"] == 1
            tv = max((vaf.get(truth_key.get(t)) or 0 for t in r["truth_ids"]), default=None) if som else None
            rows.append((r["candidate_id"] in excluded, som, r["p_somatic"] >= thr, tv, tensor_af.get(r["candidate_id"])))
    lines = [f"# E23 tensor level, chr1 test: {run.name}, {a.kind}, threshold {thr:.5f} (validation recall 0.9)", "",
             "| tensors | somatic | called | TP | precision | recall | F1 |", "|---|---|---|---|---|---|---|"]
    for name, keep in (("all", lambda e: True), ("graph-carried left out", lambda e: not e)):
        s = [x for x in rows if keep(x[0])]
        pos, called = sum(x[1] for x in s), sum(x[2] for x in s)
        tp = sum(x[1] and x[2] for x in s)
        p, r = tp / called, tp / pos
        lines.append(f"| {name} ({len(s):,}) | {pos:,} | {called:,} | {tp:,} | {p:.3f} | {r:.3f} | {f1(p, r):.3f} |")
    lines += ["", "Graph-carried left out; recall by truth VAF, precision by tensor allele fraction:", "",
              "| VAF | recall | precision | F1 |", "|---|---|---|---|"]
    s = [x for x in rows if not x[0]]
    for lo, hi, name in RANGES:
        pos = [x for x in s if x[1] and x[3] is not None and lo <= x[3] < hi]
        rec = sum(x[2] for x in pos) / len(pos) if pos else 0.0
        cal = [x for x in s if x[2] and x[4] is not None and lo <= x[4] < hi]
        pre = sum(x[1] for x in cal) / len(cal) if cal else 0.0
        lines.append(f"| {name} | {rec:.3f} ({sum(x[2] for x in pos):,} / {len(pos):,}) | "
                     f"{pre:.3f} ({sum(x[1] for x in cal):,} / {len(cal):,}) | {f1(pre, rec):.3f} |")
    Path(a.output).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

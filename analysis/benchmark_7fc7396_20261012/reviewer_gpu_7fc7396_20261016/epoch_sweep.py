"""The epoch sweep of submit_epoch_sweep.sh: COLO829T Illumina SNV, v1.1 d9, h90 recipe (1 GPU, seed 0), every epoch's
checkpoint predicted on chr1 and scored (chr1, in the BED, after the PoN). Standard library only (login node):

    PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 epoch_sweep.py

Per epoch: validation somatic F1 (the checkpoint selection score) and AP, the PASS threshold (validation recall 0.9),
chr1 tensors at or above it (all, and labelled -1 outside the confident region), and PASS F1 / precision / false
calls after the PoN with those in chr1:121-126 and 143-144 Mb. Writes epoch_sweep.md next to this script.
"""
import gzip
import json
import os
from pathlib import Path

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
HERE = Path(__file__).resolve().parent
NAME = "COLO829_Illumina_SNV_h90_sweep"
WINDOWS = [(121_000_000, 126_000_000), (143_000_000, 144_000_000)]


def main():
    rows = [json.loads(line) for line in open(RUNS / NAME / "metrics.jsonl")]
    best_epoch = max(rows, key=lambda r: r["val"]["thresholded"]["somatic"]["f1"])["epoch"]
    L = [f"# Epoch sweep of {NAME} (chr1, in the BED, after the PoN)\n",
         "(`epoch_sweep.py`; the run trained as COLO829_Illumina_SNV_h90: 1 GPU, seed 0, --norm train. * = the epoch the "
         "pipeline keeps as best.pth, highest validation somatic F1.)\n",
         "| epoch | val F1 | val AP | PASS threshold | chr1 tensors >= threshold (unlabelled) | PASS F1 (P / R) | best F1 | "
         "false calls (in the windows) |", "|---:|---:|---:|---:|---:|---|---:|---:|"]
    for r in rows:
        e = r["epoch"]
        d = RUNS / f"{NAME}_e{e:02d}"
        val = r["val"]
        f = d / "vcf_chr1" / "eval_bed" / "report.json"
        if not f.exists():
            L.append(f"| {e}{'*' if e == best_epoch else ''} | {val['thresholded']['somatic']['f1']:.4f} | "
                     f"{val['somatic_ap']:.4f} | - | - | not scored yet | - | - |")
            continue
        rep = json.load(open(f))
        thr = rep["threshold"]
        n = unl = 0
        for p in map(json.loads, gzip.open(next((d / "test_chr1").glob("*.predictions.ndjson.gz")), "rt")):
            if p["p_somatic"] >= thr:
                n += 1
                unl += p["reason"] == "outside_confident_region"
        fps = [int(line.split("\t", 2)[1]) for line in gzip.open(d / "vcf_chr1" / "eval_bed" / "vcfeval" / "pon_calls" /
                                                                    "fp.vcf.gz", "rt") if line[0] != "#"]
        win = sum(any(a <= x < b for a, b in WINDOWS) for x in fps)
        a = rep["pon"]["at_threshold"]
        L.append(f"| {e}{'*' if e == best_epoch else ''} | {val['thresholded']['somatic']['f1']:.4f} | "
                 f"{val['somatic_ap']:.4f} | {thr:.3f} | {n:,} ({unl:,}) | {a['f1']:.3f} ({a['precision']:.3f} / "
                 f"{a['recall']:.3f}) | {rep['pon']['best']['f1']:.3f} | {len(fps):,} ({win:,}) |")
    (HERE / "epoch_sweep.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

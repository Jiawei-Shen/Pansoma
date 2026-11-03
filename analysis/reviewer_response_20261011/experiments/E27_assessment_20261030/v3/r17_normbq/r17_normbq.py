#!/usr/bin/env python3
"""R1.7 / Supplementary Table R1.7b with the single-tumor SNV models trained and predicted with `--norm dataset`
(base quality z-scored per dataset; the RIS `*_SNV_d46_normbq` seed-0 runs of ris_7fc7396_20261012 on their own chr1,
and their cross-tumor predictions `reviewer_cpu_20261013/xs__*_normbq__*`, made after f344bca), next to the
`--norm train` runs (`*_7fc7396`, own chr1 scored in reviewer_cpu_20261013) the table used before. All of these call
sets (and the multi-tumor ones) were made with the off-reference rescue; it does not change F1 after the PoN. Both are seed 0, HPRC v2.1 d46, chr1, benchmark BED, after
the PoN, PASS at the checkpoint's validation recall-0.9 threshold.

Per call set: P / R / F1 (eval_bed/report.json), and recall of UV-context C>T (T[C>T]N + C[C>T]N) and of T>C
(pyrimidine-centred, from rtg tp-baseline + fn; the E18_sbs96 code). Writes tables.md and summary.json next to itself.
    /wanglab/jshen/anaconda3/bin/python r17_normbq.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1].parent  # experiments/
sys.path.insert(0, str(EXP / "E18_sbs96"))
sys.path.insert(0, str(EXP / "E16_bootstrap"))
import callsets as C  # noqa: E402
from e18_sbs96 import Genome, snvs, spectrum  # noqa: E402

R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
MODEL_PLAT = {"HG008": {"Illumina": "Illumina", "PacBio": "PacBio", "ONT": "ONT"},
              "COLO829": {"Illumina": "Illumina", "PacBio": "fiberseq", "ONT": "ONT"}}
TESTS = [("HG008T", "Illumina"), ("HG008T", "PacBio"), ("HG008T", "ONT"), ("COLO829T", "Illumina"),
         ("COLO829T", "fiberseq"), ("COLO829T", "ONT"), ("HCC1395T", "Illumina"), ("HCC1395T", "PacBio"),
         ("HCC1395T", "ONT")]
NORMS = {"train": "7fc7396", "dataset": "normbq"}


def run_dir(model, tumor, plat, norm):
    mp = MODEL_PLAT[model]["PacBio" if plat == "fiberseq" else plat]
    run = f"{model}_{mp}_SNV_d46_{NORMS[norm]}"
    if tumor == model + "T":  # own chr1: --norm train scored under reviewer_cpu_20261013, normbq in its RIS run dir
        return R / ("reviewer_cpu_20261013" if norm == "train" else "ris_7fc7396_20261012") / run, "W"
    return R / "reviewer_cpu_20261013" / f"xs__{run}__{tumor}_{plat}_d46", "H"


def score(d, genome):
    e = d / "vcf_chr1" / "eval_bed"
    t = json.load(open(e / "report.json"))["pon"]["at_threshold"]
    w = e / "vcfeval" / "pon_calls"
    out = dict(P=t["precision"], R=t["recall"], F1=t["f1"], TP=t["tp_baseline"], FP=t["fp"])
    tpb, _ = spectrum(genome, snvs(C.read_records(w / "tp-baseline.vcf.gz")))
    fn, _ = spectrum(genome, snvs(C.read_records(w / "fn.vcf.gz")))
    for name, sel in (("uv", lambda k: k[2:5] == "C>T" and k[0] in "CT"), ("tc", lambda k: k[2:5] == "T>C")):
        a = sum(v for k, v in tpb.items() if sel(k)); b = sum(v for k, v in fn.items() if sel(k))
        out[f"recall_{name}"] = a / (a + b) if a + b else None
        out[f"n_{name}"] = a + b
    return out


def main():
    genome = Genome(C.fasta_path())
    res = {}
    for tumor, plat in TESTS:
        for model in ("COLO829", "HG008"):
            for norm in NORMS:
                d, exp = run_dir(model, tumor, plat, norm)
                res[f"{model}|{tumor}_{plat}|{norm}"] = dict(exposure=exp, dir=str(d), **score(d, genome))
    (HERE / "summary.json").write_text(json.dumps(res, indent=1))
    f = lambda x: f"{x['P']:.3f} / {x['R']:.3f} / {x['F1']:.3f}"
    L = ["# R1.7b single-tumor SNV columns: `--norm train` vs `--norm dataset` (seed 0, chr1, BED, after the PoN)", "",
         "| test set | model | exposure | --norm train P / R / F1 | --norm dataset P / R / F1 | ΔF1 | recall UV C>T (train / dataset) | recall T>C (train / dataset) |",
         "|---|---|---|---|---|---|---|---|"]
    for tumor, plat in TESTS:
        for model in ("COLO829", "HG008"):
            a, b = res[f"{model}|{tumor}_{plat}|train"], res[f"{model}|{tumor}_{plat}|dataset"]
            L.append(f"| {tumor} {plat} | {model}T | [{a['exposure']}] | {f(a)} | {f(b)} | {b['F1'] - a['F1']:+.3f} | "
                     f"{a['recall_uv']:.3f} / {b['recall_uv']:.3f} ({a['n_uv']}) | {a['recall_tc']:.3f} / {b['recall_tc']:.3f} ({a['n_tc']}) |")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

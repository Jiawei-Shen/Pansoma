#!/usr/bin/env python3
"""Supplementary Table R1.7b (SNV, chr1, benchmark BED, after the common population filter, P / R / F1 at each model's
PASS threshold) on the latest model of every cell (user 2026-11-13). Every number is read from an eval_bed/report.json
(pon.at_threshold; comparators: E8 score eval_bed.pass_pon). Per cell, the newest existing run:

  six-tumor   snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, all three platforms (current SNV pipeline; ONT since
              job 477091, 2026-11-13).
  four-line   Illumina / PacBio HiFi: multicancer4_<P>_SNV_sitef_af05_lnf25w1 with the rule-2 VCFs (vcf_*_r095; COLO829T
              on the AF 0.05 test set) = the current pipeline. ONT: multicancer4_ONT_SNV normfix (OLD recipe; the
              current-pipeline retrain snv_newpipe_20261109/multicancer4_ONT_SNV_newpipe, job 477096, is used once its
              VCFs exist).
  COLO829T    Illumina / Fiber-seq: snv_newpipe_20261109/COLO829_<P>_SNV_d46_newpipe seed 0 (current pipeline) on its own
  model       chr1, and on HG008T / HCC1395T chr1 once the cross-tumor VCFs vcf_chr1_<set> exist (else the OLD normbq
              model's cross-tumor prediction). ONT: COLO829_ONT_SNV_d46_newpipe (job 481060, submitted 2026-11-13) once
              its VCFs exist, else the OLD normbq model.
  HG008T      OLD normbq models (user 2026-11-09: HG008T sets and models are not rebuilt).
  model

A cell that is not on the current pipeline is marked with a dagger. The second table lists every cell that differs from
the table in the letter (LETTER below, pasted by the user 2026-11-13).
    /wanglab/jshen/anaconda3/bin/python r17b_latest.py > r17b_latest.md
"""
import json
from pathlib import Path

R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
NP = R / "snv_newpipe_20261109"
TESTS = [(t, p) for t in ("HG008T", "COLO829T", "HCC1395T") for p in ("Illumina", "PacBio", "ONT")]
COLS = ["COLO829T model", "HG008T model", "four-line model", "six-tumor model", "DeepSomatic-TO", "ClairS-TO"]
LETTER = """0.552 / 0.603 / 0.576    0.426 / 0.852 / 0.568    0.543 / 0.684 / 0.605    0.592 / 0.611 / 0.601    0.376 / 0.889 / 0.529    0.410 / 0.845 / 0.552
0.565 / 0.548 / 0.556    0.471 / 0.838 / 0.603    0.535 / 0.720 / 0.614    0.586 / 0.670 / 0.625    0.425 / 0.806 / 0.557    0.455 / 0.844 / 0.591
0.522 / 0.703 / 0.599    0.433 / 0.859 / 0.576    0.492 / 0.726 / 0.587    0.635 / 0.551 / 0.590    0.310 / 0.869 / 0.457    0.414 / 0.831 / 0.552
0.867 / 0.712 / 0.782    0.700 / 0.738 / 0.719    0.827 / 0.564 / 0.671    0.880 / 0.661 / 0.755    0.427 / 0.835 / 0.565    0.431 / 0.764 / 0.551
0.863 / 0.707 / 0.777    0.816 / 0.702 / 0.755    0.892 / 0.446 / 0.595    0.880 / 0.691 / 0.774    0.732 / 0.555 / 0.632    0.779 / 0.745 / 0.762
0.814 / 0.700 / 0.753    0.407 / 0.716 / 0.519    0.497 / 0.308 / 0.380    0.817 / 0.670 / 0.737    0.582 / 0.695 / 0.633    0.749 / 0.684 / 0.715
0.873 / 0.662 / 0.753    0.796 / 0.804 / 0.800    0.877 / 0.816 / 0.845    0.879 / 0.822 / 0.850    0.781 / 0.927 / 0.848    0.806 / 0.859 / 0.832
0.868 / 0.647 / 0.741    0.792 / 0.625 / 0.699    0.856 / 0.787 / 0.820    0.843 / 0.809 / 0.825    0.766 / 0.927 / 0.839    0.805 / 0.822 / 0.814
0.781 / 0.551 / 0.646    0.624 / 0.690 / 0.655    0.650 / 0.668 / 0.658    0.767 / 0.697 / 0.730    0.755 / 0.766 / 0.760    0.791 / 0.730 / 0.760"""
LETTER = {(TESTS[i], COLS[j]): c.strip() for i, row in enumerate(LETTER.splitlines())
          for j, c in enumerate(row.split("    "))}


def ds(tumor, p):
    return "COLO829T_fiberseq" if (tumor, p) == ("COLO829T", "PacBio") else f"{tumor}_{p}"


def plat(tumor, p):
    return {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if tumor == "COLO829T" else "PacBio HiFi"}[p]


def rep(d):
    return json.loads((d / "eval_bed" / "report.json").read_text())["pon"]["at_threshold"]


def single(model, tumor, p):
    """(result, on the current pipeline?) of the single-tumor model of platform p on tumor's chr1."""
    mp = "fiberseq" if (model, p) == ("COLO829", "PacBio") else p
    old = f"{model}_{mp}_SNV_d46_normbq"
    own = tumor == model + "T"
    if model == "COLO829":
        d = NP / f"COLO829_{mp}_SNV_d46_newpipe" / f"vcf_chr1_{ds(tumor, p)}_d46"
        if (d / "eval_bed" / "report.json").exists():
            return rep(d), True
        assert not own or p == "ONT", d
    d = R / "ris_7fc7396_20261012" / old if own else R / "reviewer_cpu_20261013" / f"xs__{old}__{ds(tumor, p)}_d46"
    return rep(d / "vcf_chr1"), False


def four(tumor, p):
    if p == "ONT":
        d = NP / "multicancer4_ONT_SNV_newpipe" / f"vcf_chr1_{ds(tumor, p)}_d46"
        if (d / "eval_bed" / "report.json").exists():
            return rep(d), True
        return rep(R / "multicancer4_ONT_SNV" / f"vcf_chr1_{ds(tumor, p)}_d46_normfix"), False
    sfx = "_snvaf05_r095" if tumor == "COLO829T" else "_r095"
    return rep(R / f"multicancer4_{p}_SNV_sitef_af05_lnf25w1" / f"vcf_chr1_{ds(tumor, p)}_d46{sfx}"), True


def six(tumor, p):
    return rep(NP / f"multicancer6_{p}_SNV_newpipe" / f"vcf_chr1_{ds(tumor, p)}_d46"), True


def comparator(caller, tumor, p):
    return json.loads((E8 / caller / ds(tumor, p) / "SNV" / "report.json").read_text())["eval_bed"]["pass_pon"]["at_threshold"], True


def fmt(a):
    return f"{a['precision']:.3f} / {a['recall']:.3f} / {a['f1']:.3f}"


rows, changed = [], []
for t in TESTS:
    tumor, p = t
    cells = [single("COLO829", tumor, p), single("HG008", tumor, p), four(tumor, p), six(tumor, p),
             comparator("deepsomatic", tumor, p), comparator("clairs_to", tumor, p)]
    n = {a["tp_baseline"] + a["fn"] for a, _ in cells}
    assert len(n) == 1, (t, n)  # every call set is scored against the same truth
    out = []
    for col, (a, cur) in zip(COLS, cells):
        s = fmt(a)
        if s != LETTER[t, col]:
            changed.append((f"{tumor} {plat(tumor, p)}", col, LETTER[t, col], s))
        elif cur is False:
            pass
        out.append(s + ("" if cur else " †"))
    rows.append((f"{tumor} {plat(tumor, p)}", n.pop(), out))

print("# Supplementary Table R1.7b, latest models (r17b_latest.py)\n")
print("chr1, benchmark regions, after the common population filter; P / R / F1 at each model's PASS threshold. "
      "† = not on the current SNV pipeline (see the docstring of r17b_latest.py).\n")
print("| chr1 SNV | chr1 truth SNVs | " + " | ".join(COLS) + " |")
print("|---|---|" + "---|" * len(COLS))
for name, n, out in rows:
    print(f"| {name} | {n:,} | " + " | ".join(out) + " |")
print("\n## Cells that differ from the table in the letter\n")
print("| chr1 test set | column | letter | latest | ΔF1 |")
print("|---|---|---|---|---|")
for name, col, old, new in changed:
    print(f"| {name} | {col} | {old} | {new} | {float(new.split('/')[2]) - float(old.split('/')[2]):+.3f} |")

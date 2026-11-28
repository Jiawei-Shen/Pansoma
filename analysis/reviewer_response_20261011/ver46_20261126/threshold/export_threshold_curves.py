#!/usr/bin/env python3
"""Data of the "F1 against the p_somatic threshold" figure (Reviewer 3 point 9), read from eval_bed/report.json only.

chr1, benchmark regions, rtg vcfeval, after the population filter (report "pon"), HPRC v2.1, six-tumor models:
  SNV    snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe (seed 0), m6_seeds_20261123/..._s1, _s2; PASS rule 2
  INDEL  release_s0_<P>_INDEL; PASS ppr8 (Illumina, PacBio HiFi / Fiber-seq), pp (ONT)
Curve point at threshold t (E27 v3/pr_curves/plot.py at_p): the last point of pon.curve_roc with p_somatic >= t; no
record reaches t -> no calls (TP 0, FP 0, P / R / F1 0). P / R / F1 are computed from the counts.
Comparators: DeepSomatic-TO, ClairS-TO, the tool's PASS after the same filter (E8 score, eval_bed.pass_pon).
Writes threshold_curves.tsv (long), pass_points.tsv, comparators.tsv, f1_wide.tsv, export.md.
usage: python3 export_threshold_curves.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
TUMORS = ("HG008T", "COLO829T", "HCC1395T")
PLATFORMS = ("Illumina", "PacBio", "ONT")
TRUTH = {"SNV": {"HG008T": 697, "COLO829T": 2443, "HCC1395T": 3440}, "INDEL": {"HG008T": 587, "COLO829T": 131, "HCC1395T": 133}}
T = [i / 100 for i in range(1, 100)]


def ds(t, p):
    return "COLO829T_fiberseq" if (t, p) == ("COLO829T", "PacBio") else f"{t}_{p}"


def pname(t, p):
    return {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if t == "COLO829T" else "PacBio HiFi"}[p]


def call_dir(kind, t, p, seed):
    if kind == "INDEL":
        return R / f"release_s{seed}_{p}_INDEL" / f"vcf_{'pp' if p == 'ONT' else 'ppr8'}_chr1_{ds(t, p)}_d46"
    run = (R / f"snv_newpipe_20261109/multicancer6_{p}_SNV_newpipe" if seed == 0
           else R / f"m6_seeds_20261123/multicancer6_{p}_SNV_newpipe_s{seed}")
    hits = sorted(x for x in run.glob(f"vcf_chr1_{ds(t, p)}_d46*") if (x / "eval_bed/report.json").exists())
    assert len(hits) == 1, (run, t, p, hits)
    return hits[0]


def prf(tpb, tpc, fp, fn):
    p = tpc / (tpc + fp) if tpc + fp else 0.0
    r = tpb / (tpb + fn) if tpb + fn else 0.0
    return p, r, (2 * p * r / (p + r) if p + r else 0.0)


def at(curve, t, truth):
    sel = [c for c in curve if c["p_somatic"] is not None and c["p_somatic"] >= t]
    if not sel:
        return dict(tp=0, tp_call=0, fp=0, fn=truth, precision=0.0, recall=0.0, f1=0.0)
    c = sel[-1]
    assert c["tp_baseline"] + c["fn"] == truth
    p, r, f = prf(c["tp_baseline"], c["tp_call"], c["fp"], c["fn"])
    assert abs(f - c["f1"]) < 5.1e-5, c
    return dict(tp=c["tp_baseline"], tp_call=c["tp_call"], fp=c["fp"], fn=c["fn"], precision=p, recall=r, f1=f)


long_rows, pass_rows, wide = [], [], []
for kind, seeds in (("SNV", (0, 1, 2)), ("INDEL", (0,))):
    for t in TUMORS:
        for p in PLATFORMS:
            for seed in seeds:
                d = call_dir(kind, t, p, seed)
                j = json.loads((d / "eval_bed/report.json").read_text())
                truth = TRUTH[kind][t]
                a = j["pon"]["at_threshold"]
                assert a["tp_baseline"] + a["fn"] == truth, d
                P, Rc, F = prf(a["tp_baseline"], a["tp_call"], a["fp"], a["fn"])
                curve = j["pon"]["curve_roc"]
                thr = j["threshold"]
                c_thr, c05, b = at(curve, thr, truth), at(curve, 0.5, truth), j["pon"]["best"]
                model = "six-tumor"
                pass_rows.append([kind, t, pname(t, p), model, seed, f"{thr:.4f}", a["tp_baseline"], a["fp"], a["fn"], f"{P:.4f}", f"{Rc:.4f}",
                                  f"{F:.4f}", c05["tp"], c05["fp"], c05["fn"], f"{c05['precision']:.4f}", f"{c05['recall']:.4f}", f"{c05['f1']:.4f}",
                                  f"{c_thr['f1']:.4f}", f"{b['p_somatic']:.4f}", f"{b['f1']:.4f}", truth, str(d.relative_to(R))])
                fs = []
                for x in T:
                    c = at(curve, x, truth)
                    long_rows.append([kind, t, pname(t, p), model, seed, f"{x:.2f}", c["tp"], c["fp"], c["fn"], f"{c['precision']:.4f}",
                                      f"{c['recall']:.4f}", f"{c['f1']:.4f}"])
                    fs.append(c["f1"])
                wide.append([kind, t, pname(t, p), seed, fs])


def tsv(name, header, rows):
    (HERE / name).write_text("\n".join("\t".join(map(str, r)) for r in [header] + rows) + "\n")


tsv("threshold_curves.tsv", ["kind", "tumor", "platform", "model", "seed", "threshold", "tp", "fp", "fn", "precision", "recall", "f1"], long_rows)
PH = ["kind", "tumor", "platform", "model", "seed", "pass_threshold", "pass_tp", "pass_fp", "pass_fn", "pass_precision", "pass_recall", "pass_f1",
      "t0.5_tp", "t0.5_fp", "t0.5_fn", "t0.5_precision", "t0.5_recall", "t0.5_f1", "curve_f1_at_pass_threshold", "best_f1_threshold", "best_f1",
      "truth", "call_set"]
tsv("pass_points.tsv", PH, pass_rows)
tsv("f1_wide.tsv", ["kind", "tumor", "platform", "seed"] + [f"{x:.2f}" for x in T],
    [[k, t, p, s] + [f"{x:.3f}" for x in fs] for k, t, p, s, fs in wide])

comp = []
for kind in ("SNV", "INDEL"):
    for t in TUMORS:
        for p in PLATFORMS:
            row = [kind, t, pname(t, p)]
            for caller in ("deepsomatic", "clairs_to"):
                f = E8 / caller / ds(t, p) / kind / "report.json"
                a = json.loads(f.read_text())["eval_bed"]["pass_pon"]["at_threshold"]
                assert a["tp_baseline"] + a["fn"] == TRUTH[kind][t], f
                P, Rc, F = prf(a["tp_baseline"], a["tp_call"], a["fp"], a["fn"])
                row += [a["tp_baseline"], a["fp"], a["fn"], f"{P:.4f}", f"{Rc:.4f}", f"{F:.4f}"]
            comp.append(row)
tsv("comparators.tsv", ["kind", "tumor", "platform"] + [f"{c}_{k}" for c in ("DeepSomatic-TO", "ClairS-TO") for k in ("tp", "fp", "fn", "precision", "recall", "f1")], comp)
print(len(long_rows), "curve rows;", len(pass_rows), "call sets;", len(comp), "comparator rows")

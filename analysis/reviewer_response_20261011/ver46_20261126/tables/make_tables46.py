#!/usr/bin/env python3
"""Letter ver4.6 tables that come straight from eval_bed/report.json (chr1, benchmark regions, rtg vcfeval, PASS).

  t12   Response Tables 1.2 = 2.10 = 3.4: INDELs before the population filter, HPRC v1.1 / v2.1 columns on the
        six-tumor INDEL models WITHOUT site features (m6_nositef_graphs_20261118/multicancer6_<P>_INDEL_{d9,d46}_nositef,
        seed 0; PASS ppr8, ONT pp); the GRCh38-only column (multicancer6_<P>_INDEL_linear_nositef, seed 0) recomputed
        and compared with the letter's current cells.
  t32   Response Table 3.2, SNV rows: the six-tumor SNV models on HPRC v2.1, seeds 0 / 1 / 2, after the filter
        (seed 0 snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe; seeds 1, 2 m6_seeds_20261123/..._s<SEED>).
P / R / F1 are computed from the counts (P = TP calls / (TP calls + FP), R = TP / truth), not from the report's
4-decimal values. The truth count of every cell of a dataset is asserted equal.
usage: python3 make_tables46.py   -> tables46.json, tables46.md
"""
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
TUMORS = ("HG008T", "COLO829T", "HCC1395T")
PLATFORMS = ("Illumina", "PacBio", "ONT")
TRUTH = {"SNV": {"HG008T": 697, "COLO829T": 2443, "HCC1395T": 3440}, "INDEL": {"HG008T": 587, "COLO829T": 131, "HCC1395T": 133}}


def pname(t, p):
    return {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if t == "COLO829T" else "PacBio HiFi"}[p]


def prf(a):
    tpb, tpc, fp, fn = a["tp_baseline"], a["tp_call"], a["fp"], a["fn"]
    p = tpc / (tpc + fp) if tpc + fp else 0.0
    r = tpb / (tpb + fn)
    f = 2 * p * r / (p + r) if p + r else 0.0
    for k, v in (("precision", p), ("recall", r), ("f1", f)):
        assert abs(v - a[k]) < 5.1e-5, (k, v, a)
    return {"precision": p, "recall": r, "f1": f, "tp": tpb, "tp_call": tpc, "fp": fp, "fn": fn, "truth": tpb + fn}


def rep(d, stage):
    return prf(json.loads((d / "eval_bed/report.json").read_text())[stage]["at_threshold"])


def cell3(x):
    return f"{x['precision']:.3f} / {x['recall']:.3f} / {x['f1']:.3f}"


T, L = {}, []

# ------------------------------------------------------------------ t12: INDEL before the filter
# the letter's current GRCh38-only cells (ver4.5), to be reproduced
LETTER_LINEAR = {("HG008T", "Illumina"): "0.016 / 0.445 / 0.031", ("HG008T", "PacBio"): "0.015 / 0.475 / 0.029",
                 ("HG008T", "ONT"): "0.015 / 0.097 / 0.027", ("COLO829T", "Illumina"): "0.006 / 0.298 / 0.013",
                 ("COLO829T", "PacBio"): "0.015 / 0.198 / 0.027", ("COLO829T", "ONT"): "0.006 / 0.008 / 0.006",
                 ("HCC1395T", "Illumina"): "0.015 / 0.744 / 0.029", ("HCC1395T", "PacBio"): "0.010 / 0.797 / 0.020",
                 ("HCC1395T", "ONT"): "0.019 / 0.534 / 0.036"}
# comparator cells of the letter's table (unchanged; used for the counts in the notes): ClairS-TO, DeepSomatic-TO F1
COMP_F1 = {("HG008T", "Illumina"): (0.007, 0.032), ("HG008T", "PacBio"): (0.007, 0.051), ("HG008T", "ONT"): (0.007, 0.026),
           ("COLO829T", "Illumina"): (0.003, 0.019), ("COLO829T", "PacBio"): (0.004, 0.105), ("COLO829T", "ONT"): (0.003, 0.016),
           ("HCC1395T", "Illumina"): (0.012, 0.042), ("HCC1395T", "PacBio"): (0.011, 0.382), ("HCC1395T", "ONT"): (0.012, 0.328)}


def indel_dir(arm, p, t):
    r = "pp" if p == "ONT" else "ppr8"
    if arm == "linear":
        return R / f"m4m6_graphs_20261110/multicancer6_{p}_INDEL_linear_nositef/vcf_{r}_chr1_{t}"
    return R / f"m6_nositef_graphs_20261118/multicancer6_{p}_INDEL_{arm}_nositef/vcf_{r}_chr1_{t}"


rows, cells = [], {}
for t in TUMORS:
    for p in PLATFORMS:
        c = {arm: rep(indel_dir(arm, p, t), "raw") for arm in ("linear", "d9", "d46")}
        a = {arm: rep(indel_dir(arm, p, t), "pon") for arm in ("linear", "d9", "d46")}
        assert {x["truth"] for x in list(c.values()) + list(a.values())} == {TRUTH["INDEL"][t]}, (t, p)
        assert cell3(c["linear"]) == LETTER_LINEAR[t, p], (t, p, cell3(c["linear"]), LETTER_LINEAR[t, p])
        cells[f"{t} {pname(t, p)}"] = {"before": c, "after": a}
        rows.append([f"{t} {pname(t, p)}", cell3(c["linear"]), cell3(c["d9"]), cell3(c["d46"])])
T["t12"] = {"header": ["chr1 INDELs, no population filter (P / R / F1)", "Pansoma GRCh38-only", "Pansoma HPRC v1.1", "Pansoma HPRC v2.1"],
            "rows": rows, "cells": cells}
f = {k: {arm: v["before"][arm]["f1"] for arm in ("linear", "d9", "d46")} for k, v in cells.items()}
keys = [(t, p) for t in TUMORS for p in PLATFORMS]
notes = {
    "HPRC v1.1 F1 > GRCh38-only": sum(x["d9"] > x["linear"] for x in f.values()),
    "HPRC v2.1 F1 > GRCh38-only": sum(x["d46"] > x["linear"] for x in f.values()),
    "HPRC v2.1 F1 > HPRC v1.1": sum(x["d46"] > x["d9"] for x in f.values()),
    "HPRC v2.1 F1 > both comparators (3 decimals)": sum(round(f[f"{t} {pname(t, p)}"]["d46"], 3) > max(COMP_F1[t, p]) for t, p in keys),
    "HPRC v1.1 F1 > both comparators (3 decimals)": sum(round(f[f"{t} {pname(t, p)}"]["d9"], 3) > max(COMP_F1[t, p]) for t, p in keys),
}
for tum in TUMORS:
    for arm, nm in (("linear", "GRCh38-only"), ("d9", "HPRC v1.1"), ("d46", "HPRC v2.1")):
        v = [x[arm] for k, x in f.items() if k.startswith(tum)]
        notes[f"{tum} F1 range {nm}"] = f"{min(v):.3f}-{max(v):.3f}"
T["t12"]["notes"] = notes
L += ["## t12. Response Tables 1.2 = 2.10 = 3.4: chr1 INDELs before the population filter, P / R / F1 (seed 0, no site features in any model)", "",
      "| " + " | ".join(T["t12"]["header"]) + " |", "|---|---|---|---|"]
L += ["| " + " | ".join(r) + " |" for r in rows]
L += [""] + [f"- {k}: {v}{' of 9' if isinstance(v, int) else ''}" for k, v in notes.items()] + [""]

# ------------------------------------------------------------------ t32: six-tumor SNV, seeds 0 / 1 / 2, after the filter


def snv_dir(seed, p, t):
    run = (R / f"snv_newpipe_20261109/multicancer6_{p}_SNV_newpipe" if seed == 0
           else R / f"m6_seeds_20261123/multicancer6_{p}_SNV_newpipe_s{seed}")
    stem = "COLO829T_fiberseq" if (t, p) == ("COLO829T", "PacBio") else f"{t}_{p}"
    hits = sorted(x for x in run.glob(f"vcf_chr1_{stem}_d46*") if (x / "eval_bed/report.json").exists())
    assert len(hits) == 1, (run, stem, hits)
    return hits[0]


def seeds_cell(c, k):
    v = [x[k] for x in c]
    return f"{' / '.join(f'{x:.3f}' for x in v)} ({statistics.mean(v):.3f} ± {statistics.stdev(v):.3f})"


rows, sd, per = [], {}, {}
for t in TUMORS:
    for p in PLATFORMS:
        c = [rep(snv_dir(s, p, t), "pon") for s in (0, 1, 2)]
        assert {x["truth"] for x in c} == {TRUTH["SNV"][t]}, (t, p)
        rows.append([t, pname(t, p), "6 tumors, SNV"] + [seeds_cell(c, k) for k in ("precision", "recall", "f1")])
        sd[f"{t} {pname(t, p)}"] = statistics.stdev(x["f1"] for x in c)
        per[f"{t} {pname(t, p)}"] = c
s = sorted(sd.values())
T["t32_snv"] = {"header": ["tumor", "platform", "training set", "precision", "recall", "F1"], "rows": rows, "sd": sd, "cells": per,
                "sd_range": f"{s[0]:.3f}–{s[-1]:.3f}", "sd_median": statistics.median(s),
                "sd_max3": sorted(sd.items(), key=lambda kv: -kv[1])[:3]}
L += ["## t32. Response Table 3.2, SNV rows: six-tumor SNV models on HPRC v2.1, seeds 0 / 1 / 2 (mean ± SD), chr1 after the population filter", "",
      "| " + " | ".join(T["t32_snv"]["header"]) + " |", "|---|---|---|---|---|---|"]
L += ["| " + " | ".join(r) + " |" for r in rows]
L += ["", f"- SD of F1 over the three seeds: {T['t32_snv']['sd_range']} (median {T['t32_snv']['sd_median']:.3f}); largest: "
      + "; ".join(f"{k} {v:.3f} ({' / '.join(format(x['f1'], '.3f') for x in per[k])})" for k, v in T["t32_snv"]["sd_max3"]), ""]

(HERE / "tables46.json").write_text(json.dumps(T, indent=1))
(HERE / "tables46.md").write_text("\n".join(L) + "\n")
print("\n".join(L))

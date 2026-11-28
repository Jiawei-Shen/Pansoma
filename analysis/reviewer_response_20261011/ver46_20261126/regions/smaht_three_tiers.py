#!/usr/bin/env python3
"""SMaHT easy / difficult / extreme strata for the three-reference call sets of Response Tables 1.5 = 2.5 = 3.15
(2026-11-27; the letter's table has only the extreme tier).

Same call sets, scoring, block bootstrap and metric rules as region_nositef.py (six-tumor models, seed 0, no model with
the population site features as an input; chr1; rtg vcfeval in the benchmark BED; raw_calls = before the population
filter, pon_calls = after; paired dF1 = HPRC - GRCh38-only over 2,000 resamples of 1-Mb blocks, seed 20261015, 95 %
percentile interval). Only the strata differ: the three SMaHT v2 BEDs (disjoint; they tile chr1).
Checks (asserted): totals of every rtg directory = report.json; easy + difficult + extreme = all, for every count of
every call set; the extreme tier equals cells_six_tumor.json (region_nositef.py) cell by cell, SNV counts, dF1 and
interval; truth per tier identical for the three references of a dataset.

usage: /wanglab/jshen/anaconda3/bin/python smaht_three_tiers.py
Writes smaht_three_tiers.md, smaht_three_tiers.tsv, smaht_three_tiers.json next to this script.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments")
sys.path.insert(0, str(EXP / "E16_bootstrap"))
sys.path.insert(0, str(EXP / "E10b_giab_strata"))
import callsets as C  # noqa: E402
import e16_bootstrap as E  # noqa: E402
import giab_strata as G  # noqa: E402

SMAHT = "/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_{}_v2.bed"
TIERS = ("easy", "difficult", "extreme")
RESAMPLES, SEED = 2000, 20261015
GNAME = {"linear": "GRCh38-only", "d9": "HPRC v1.1", "d46": "HPRC v2.1"}
ARMS = ("linear", "d9", "d46")
PARTS = ("tp-baseline", "fn", "tp", "fp")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
G6 = RUNS / "m4m6_graphs_20261110"
NS = RUNS / "m6_nositef_graphs_20261118"
SETS = [(t, p) for t in ("HG008T", "COLO829T", "HCC1395T") for p in ("Illumina", "PacBio", "ONT")]
WHICH = (("raw_calls", "raw", "before"), ("pon_calls", "pon", "after"))
MINUS = "−"


def ds(t, p):
    return "COLO829T_fiberseq" if (t, p) == ("COLO829T", "PacBio") else f"{t}_{p}"


def lab(t, p):
    return f"{t} " + {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if t == "COLO829T" else "PacBio HiFi"}[p]


def vdir(kind, t, p, arm):
    r = "" if kind == "SNV" else ("pp_" if p == "ONT" else "ppr8_")
    if arm == "linear":
        return G6 / f"multicancer6_{p}_{kind}_linear_nositef" / f"vcf_{r}chr1_{t}"
    run = NS / f"multicancer6_{p}_{kind}_{arm}_nositef"
    if kind == "INDEL":
        return run / f"vcf_{r}chr1_{t}"
    hits = sorted(run.glob(f"vcf_chr1_{ds(t, p)}*"))
    assert len(hits) == 1, (run, t, p, hits)
    return hits[0]


def bed_chr1(path):
    iv = sorted((int(c[1]), int(c[2])) for c in (ln.split("\t")[:3] for ln in open(path)) if c[0] == "chr1")
    m = []
    for s, e in iv:
        if m and s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    a = np.array(m, dtype=np.int64).reshape(-1, 2)
    return a[:, 0], a[:, 1]


def metrics(tot):
    tot = np.asarray(tot, dtype=float)
    tpb, fn, tpc, fp = tot[..., 0], tot[..., 1], tot[..., 2], tot[..., 3]
    with np.errstate(invalid="ignore", divide="ignore"):
        p = tpc / (tpc + fp)
        r = tpb / (tpb + fn)
        f = 2 * p * r / (p + r)
    f = np.where(np.isnan(f) & (tpb + tpc == 0) & (fn + fp > 0), 0.0, f)
    return p, r, f


ivs = {t: bed_chr1(SMAHT.format(t)) for t in TIERS}
bp = {t: int((ivs[t][1] - ivs[t][0]).sum()) for t in TIERS}
bs = E.BlockSets(C.fasta_path(), RESAMPLES, SEED)


def score(d, which):
    d = Path(d)
    rep = json.loads((d / "eval_bed/report.json").read_text())
    w = d / "eval_bed/vcfeval" / which
    assert C._rtg_complete(w), w
    recs = {p: C.read_records(w / f"{p}.vcf.gz") for p in PARTS}
    assert all(r[0] == "chr1" for p in PARTS for r in recs[p]), w
    rec = {p: np.array([r[1] for r in recs[p]], dtype=np.int64) for p in PARTS}
    stage = rep["raw" if which == "raw_calls" else "pon"]
    a = stage["at_threshold"]
    n = (len(rec["tp-baseline"]), len(rec["fn"]), len(rec["tp"]), len(rec["fp"]))
    assert n == (a["tp_baseline"], a["fn"], a["tp_call"], a["fp"]), (w, n, a)
    blocks = bs.get("bed", rep["bed"], G.WIDTH)
    out = {}
    for name, iv in list(ivs.items()) + [("all", None)]:
        sub = {p: (rec[p] if iv is None else rec[p][G.inside(rec[p], iv)]).tolist() for p in PARTS}
        m = E.block_counts(sub, blocks)
        assert m.sum() == sum(len(v) for v in sub.values()), (w, name)
        out[name] = (m.sum(axis=0), blocks.mult @ m, blocks.key)
    assert (sum(out[t][0] for t in TIERS) == out["all"][0]).all(), (w, [out[t][0] for t in TIERS], out["all"][0])  # tiers tile chr1
    return out, {"dir": str(d), "bed": rep["bed"], "truth": rep["truth"]}


def paired(a, b):
    assert a[2] == b[2]
    d = metrics(a[1])[2] - metrics(b[1])[2]
    lo, hi, _, nnan = E.ci(d)
    return {"d": float(metrics(a[0])[2] - metrics(b[0])[2]), "lo": lo, "hi": hi, "dropped": nnan}


def f3(x):
    return "-" if x != x else f"{x:.3f}"


def sgn(x, plus=True):
    return (MINUS if x < 0 else ("+" if plus else "")) + f"{abs(x):.3f}"


def ci_cell(c):
    return f"{sgn(c['d'])} [{sgn(c['lo'], False)}, {sgn(c['hi'], False)}]"


res, meta = {}, {}
for kind in ("SNV", "INDEL"):
    for t, p in SETS:
        for arm in ARMS:
            for w, _, _ in WHICH:
                res[(kind, t, p, arm, w)], meta[(kind, t, p, arm, w)] = score(vdir(kind, t, p, arm), w)
        ms = [meta[(kind, t, p, arm, w)] for arm in ARMS for w, _, _ in WHICH]
        assert len({m["bed"] for m in ms}) == 1 and len({m["truth"] for m in ms}) == 1, (kind, t, p)
        for st in TIERS + ("all",):
            assert len({int(res[(kind, t, p, arm, w)][st][0][:2].sum()) for arm in ARMS for w, _, _ in WHICH}) == 1, (kind, t, p, st)
    print(kind, "scored", file=sys.stderr)


def cell(kind, t, p, arm, w, st):
    tot = res[(kind, t, p, arm, w)][st][0]
    pr, rc, f1 = (float(x) for x in metrics(tot))
    return {"truth": int(tot[0] + tot[1]), "tp": int(tot[0]), "fn": int(tot[1]), "tp_call": int(tot[2]), "fp": int(tot[3]),
            "calls": int(tot[2] + tot[3]), "precision": pr, "recall": rc, "f1": f1}


# ---- check: the extreme tier equals region_nositef.py's cells_six_tumor.json
old = json.loads((HERE / "cells_six_tumor.json").read_text())["cells"]
nchk = 0
for w, _, when in WHICH:
    for t, p in SETS:
        for arm in ("d9", "d46"):
            o = old[f"{when}|{ds(t, p)}|{arm}"]["values"]
            for k, a in (("extreme_linear", "linear"), ("extreme_hprc", arm)):
                c = cell("SNV", t, p, a, w, "extreme")
                assert all(c[x] == o[k][x] for x in ("truth", "calls", "fp", "tp_call")), (when, t, p, arm, k)
            d = paired(res[("SNV", t, p, arm, w)]["extreme"], res[("SNV", t, p, "linear", w)]["extreme"])
            assert all(abs(d[x] - o["extreme"][x]) < 1e-12 for x in ("d", "lo", "hi")), (when, t, p, arm, d, o["extreme"])
            nchk += 1

rows, J = [], {}
for kind in ("SNV", "INDEL"):
    for w, _, when in WHICH:
        for t, p in SETS:
            for st in TIERS + ("all",):
                c = {a: cell(kind, t, p, a, w, st) for a in ARMS}
                d = {a: paired(res[(kind, t, p, a, w)][st], res[(kind, t, p, "linear", w)][st]) for a in ("d9", "d46")}
                J[f"{kind}|{when}|{ds(t, p)}|{st}"] = {"label": lab(t, p), "arms": {GNAME[a]: c[a] for a in ARMS},
                                                         "dF1": {GNAME[a]: d[a] for a in d}}
                for a in ARMS:
                    r = {"kind": kind, "filter": when, "tumor": t, "platform": lab(t, p).split(" ", 1)[1], "tier": st, "reference": GNAME[a],
                         **{k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in c[a].items()}}
                    r.update({"dF1": "", "dF1_lo": "", "dF1_hi": ""} if a == "linear" else
                             {"dF1": f"{d[a]['d']:.4f}", "dF1_lo": f"{d[a]['lo']:.4f}", "dF1_hi": f"{d[a]['hi']:.4f}"})
                    rows.append(r)
with open(HERE / "smaht_three_tiers.tsv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    wr.writeheader()
    wr.writerows(rows)
(HERE / "smaht_three_tiers.json").write_text(json.dumps({"chr1_bp": bp, "cells": J}, indent=1) + "\n")

L = ["# SMaHT easy / difficult / extreme: chr1 calling on the three references (six-tumor models without site features)", "",
     "Generated by `smaht_three_tiers.py`. chr1, benchmark regions, PASS calls, seed 0; no model takes the population site features as an input. "
     "Call sets, scoring and bootstrap as `region_nositef.py` (Response Tables 1.5 = 2.5 = 3.15). ΔF1 = HPRC − GRCh38-only with the 95 % "
     "percentile interval of 2,000 resamples of 1-Mb blocks. SMaHT v2 tiers on chr1: " +
     ", ".join(f"{t} {bp[t] / 1e6:.1f} Mb" for t in TIERS) + " (disjoint; before the benchmark BED).", ""]
for kind in ("SNV", "INDEL"):
    for w, _, when in WHICH:
        L += [f"## chr1 {kind}, {when} the population filter", "",
              "| chr1 dataset | HPRC graph | " + " | ".join(f"SMaHT {t}: chr1 truth" + f" | SMaHT {t}: F1 GRCh38-only → HPRC | SMaHT {t}: ΔF1 [95% CI]" for t in TIERS) + " |",
              "|---|---|" + "---|" * 9]
        for t, p in SETS:
            for arm in ("d9", "d46"):
                cs = []
                for st in TIERS:
                    x = J[f"{kind}|{when}|{ds(t, p)}|{st}"]
                    cs += [f"{x['arms']['GRCh38-only']['truth']:,}", f"{f3(x['arms']['GRCh38-only']['f1'])} → {f3(x['arms'][GNAME[arm]]['f1'])}", ci_cell(x["dF1"][GNAME[arm]])]
                L.append(f"| {lab(t, p)} | {GNAME[arm][5:]} | " + " | ".join(cs) + " |")
        L += ["", f"### chr1 {kind}, {when} the filter: PASS calls (precision) per tier, GRCh38-only → v1.1 → v2.1", "",
              "| chr1 dataset | " + " | ".join(f"SMaHT {t}" for t in TIERS) + " |", "|---|---|---|---|"]
        for t, p in SETS:
            L.append(f"| {lab(t, p)} | " + " | ".join(" → ".join(f"{x['arms'][GNAME[a]]['calls']:,} ({f3(x['arms'][GNAME[a]]['precision'])})" for a in ARMS)
                                                    for x in (J[f"{kind}|{when}|{ds(t, p)}|{st}"] for st in TIERS)) + " |")
        L += ["", "| tier | pairs | ΔF1 range | ΔF1 > 0 | interval above zero | interval below zero |", "|---|---|---|---|---|---|"]
        for st in TIERS + ("all",):
            v = [J[f"{kind}|{when}|{ds(t, p)}|{st}"]["dF1"][GNAME[a]] for t, p in SETS for a in ("d9", "d46")]
            L.append(f"| SMaHT {st} | {len(v)} | {sgn(min(x['d'] for x in v))} to {sgn(max(x['d'] for x in v))} | {sum(x['d'] > 0 for x in v)} | "
                     f"{sum(x['lo'] > 0 for x in v)} | {sum(x['hi'] < 0 for x in v)} |")
        L.append("")
L += ["## Checks", "", f"- Extreme tier = `cells_six_tumor.json` (SNV counts, ΔF1 and interval) in {nchk} of {nchk} cells, before and after the filter.",
      "- easy + difficult + extreme = all, for TP / FN / TP-call / FP of each of the 108 rtg directories; totals = report.json.",
      "- Truth per tier identical for the three references of a dataset."]
(HERE / "smaht_three_tiers.md").write_text("\n".join(L) + "\n")

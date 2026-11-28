#!/usr/bin/env python3
"""ClairS-TO and DeepSomatic-TO by region, before and after population filtering, next to Pansoma (2026-11-27).

Strata: the five SNV strata of Response Tables 1.5 = 2.5 = 3.15 (SMaHT extreme; GIAB v3.5 all difficult, low mappability
or segmental duplications, tandem repeats and homopolymers, outside all difficult), SMaHT easy and difficult, and for
INDELs GIAB homopolymers >= 7 bp. chr1, benchmark BED, PASS calls; rtg vcfeval records of the existing score
directories (nothing is re-called or re-scored).
Call sets (E8 = /scratch/jshen/data/reviewer_experiments_20261013/E8_comparators; R = /scratch/jshen/data/pansoma_net_v2_runs):
  before the filter (the call sets of Response Tables 1.1 / 1.2)
    ClairS-TO        E8/score_nopon_rerun/clairs_to/<DS>/<KIND>/eval_bed/pass   (rerun with --disable_nonsomatic_tagging)
    DeepSomatic-TO   E8/score_nopon/deepsomatic/<DS>/<KIND>/eval_bed/pass       (its PON filter lifted; population AF still an input)
    Pansoma          six-tumor models without site features on GRCh38-only / HPRC v1.1 / HPRC v2.1, raw_calls
  after the filter (every caller after the same filter with our resource VCFs)
    ClairS-TO, DeepSomatic-TO   E8/score/<caller>/<DS>/<KIND>/eval_bed/pass_pon  (default run, then our filter)
    Pansoma                     the same three no-site-feature models, pon_calls, and the release models on HPRC v2.1
                                (site features as an input; the models compared with the callers in the letter):
                                SNV R/snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, INDEL R/release_s0_<P>_INDEL
Paired dF1 = Pansoma HPRC v2.1 - comparator (before: the no-site-feature model; after: the release model), 95 %
percentile interval of 2,000 resamples of 1-Mb blocks (seed 20261015, same resamples for both callers).
Checks (asserted): totals of every directory = its report.json; same truth VCF and BED for all callers of a dataset;
truth per stratum identical for all callers of a dataset; SMaHT tiers sum to the whole.

usage: /wanglab/jshen/anaconda3/bin/python comparators_regions.py
Writes comparators_regions.md / .tsv / .json next to this script.
"""
import csv
import json
import os
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
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
G6, NS = RUNS / "m4m6_graphs_20261110", RUNS / "m6_nositef_graphs_20261118"
RESAMPLES, SEED = 2000, 20261015
HP7 = "AllHomopolymers_ge7bp_imperfectge11bp_slop5"
ST = [("easy", "SMaHT easy"), ("difficult", "SMaHT difficult"), ("extreme", "SMaHT extreme"),
      ("alldifficultregions", "GIAB all difficult"), ("alllowmapandsegdupregions", "GIAB low map or segdup"),
      ("AllTandemRepeatsandHomopolymers_slop5", "GIAB tandem repeats and homopolymers"),
      ("notinalldifficultregions", "GIAB outside difficult"), (HP7, "GIAB homopolymers ≥ 7 bp"), ("all", "all benchmark regions")]
SNV_COLS = [s for s, _ in ST if s not in (HP7, "all")]
SNAME = dict(ST)
PARTS = ("tp-baseline", "fn", "tp", "fp")
SETS = [(t, p) for t in ("HG008T", "COLO829T", "HCC1395T") for p in ("Illumina", "PacBio", "ONT")]
MINUS = "−"
CALLERS = {"before": ["ClairS-TO", "DeepSomatic-TO", "Pansoma GRCh38-only", "Pansoma HPRC v1.1", "Pansoma HPRC v2.1"],
           "after": ["ClairS-TO", "DeepSomatic-TO", "Pansoma GRCh38-only", "Pansoma HPRC v1.1", "Pansoma HPRC v2.1", "Pansoma release (HPRC v2.1)"]}
REF = {"before": "Pansoma HPRC v2.1", "after": "Pansoma release (HPRC v2.1)"}
ARM = {"Pansoma GRCh38-only": "linear", "Pansoma HPRC v1.1": "d9", "Pansoma HPRC v2.1": "d46"}


def ds(t, p):
    return "COLO829T_fiberseq" if (t, p) == ("COLO829T", "PacBio") else f"{t}_{p}"


def lab(t, p):
    return f"{t} " + {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if t == "COLO829T" else "PacBio HiFi"}[p]


def source(caller, kind, t, p, when):
    """(rtg directory, report.json, key path of at_threshold in the report)."""
    if caller in ("ClairS-TO", "DeepSomatic-TO"):
        c = "clairs_to" if caller == "ClairS-TO" else "deepsomatic"
        src = "score" if when == "after" else ("score_nopon_rerun" if c == "clairs_to" else "score_nopon")
        sub = "pass_pon" if when == "after" else "pass"
        d = E8 / src / c / ds(t, p) / kind
        return d / "eval_bed" / sub, d / "report.json", ("eval_bed", sub)
    r = "" if kind == "SNV" else ("pp_" if p == "ONT" else "ppr8_")
    if caller.startswith("Pansoma release"):
        d = (RUNS / f"snv_newpipe_20261109/multicancer6_{p}_SNV_newpipe" if kind == "SNV" else RUNS / f"release_s0_{p}_INDEL") / f"vcf_{r}chr1_{ds(t, p)}_d46"
    else:
        arm = ARM[caller]
        if arm == "linear":
            d = G6 / f"multicancer6_{p}_{kind}_linear_nositef" / f"vcf_{r}chr1_{t}"
        elif kind == "INDEL":
            d = NS / f"multicancer6_{p}_INDEL_{arm}_nositef" / f"vcf_{r}chr1_{t}"
        else:
            hits = sorted((NS / f"multicancer6_{p}_SNV_{arm}_nositef").glob(f"vcf_chr1_{ds(t, p)}*"))
            assert len(hits) == 1, hits
            d = hits[0]
    w = "pon_calls" if when == "after" else "raw_calls"
    return d / "eval_bed/vcfeval" / w, d / "eval_bed/report.json", ("pon" if when == "after" else "raw",)


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


ivs = {t: bed_chr1(SMAHT.format(t)) for t in ("easy", "difficult", "extreme")}
ivs.update({name: G.chr1_intervals(G.GIAB / rel) for name, rel, _ in G.STRATA if name in SNAME})
assert set(ivs) | {"all"} == set(SNAME), set(SNAME) - set(ivs)
bs = E.BlockSets(C.fasta_path(), RESAMPLES, SEED)
BED = {}


def score(caller, kind, t, p, when):
    w, repf, keys = source(caller, kind, t, p, when)
    rep = json.loads(repf.read_text())
    assert C._rtg_complete(w), w
    recs = {x: C.read_records(w / f"{x}.vcf.gz") for x in PARTS}
    assert all(r[0] == "chr1" for x in PARTS for r in recs[x]), w
    rec = {x: np.array([r[1] for r in recs[x]], dtype=np.int64) for x in PARTS}
    stage = rep
    for k in keys:
        stage = stage[k]
    a = stage["at_threshold"]
    n = (len(rec["tp-baseline"]), len(rec["fn"]), len(rec["tp"]), len(rec["fp"]))
    assert n == (a["tp_baseline"], a["fn"], a["tp_call"], a["fp"]), (w, n, a)
    key = (kind, t, p)                     # one truth VCF and BED per dataset and kind, for every caller
    ref = BED.setdefault(key, (rep["bed"], rep["truth"]))
    assert os.path.samefile(ref[0], rep["bed"]) and os.path.samefile(ref[1], rep["truth"]), (w, ref, rep["bed"], rep["truth"])
    blocks = bs.get("bed", ref[0], G.WIDTH)
    out = {}
    for name in SNAME:
        iv = ivs.get(name)
        sub = {x: (rec[x] if iv is None else rec[x][G.inside(rec[x], iv)]).tolist() for x in PARTS}
        m = E.block_counts(sub, blocks)
        assert m.sum() == sum(len(v) for v in sub.values()), (w, name)
        out[name] = (m.sum(axis=0), blocks.mult @ m, blocks.key)
    assert (sum(out[x][0] for x in ("easy", "difficult", "extreme")) == out["all"][0]).all(), w
    return out, str(w)


def paired(a, b):
    assert a[2] == b[2]
    d = metrics(a[1])[2] - metrics(b[1])[2]
    lo, hi, _, nnan = E.ci(d)
    return {"d": float(metrics(a[0])[2] - metrics(b[0])[2]), "lo": lo, "hi": hi, "dropped": nnan}


def f3(x):
    return "-" if x != x else f"{x:.3f}"


def sgn(x, plus=True):
    return (MINUS if x < 0 else ("+" if plus else "")) + f"{abs(x):.3f}"


res, dirs = {}, {}
for kind in ("SNV", "INDEL"):
    for when in ("before", "after"):
        for t, p in SETS:
            for c in CALLERS[when]:
                res[(kind, when, t, p, c)], dirs[f"{kind}|{when}|{ds(t, p)}|{c}"] = score(c, kind, t, p, when)
    for t, p in SETS:                      # same truth per stratum for every caller, before and after
        for st in SNAME:
            n = {int(res[(kind, when, t, p, c)][st][0][:2].sum()) for when in ("before", "after") for c in CALLERS[when]}
            assert len(n) == 1, (kind, t, p, st, n)
    print(kind, "scored", file=sys.stderr)


def cell(kind, when, t, p, c, st):
    tot = res[(kind, when, t, p, c)][st][0]
    pr, rc, f1 = (float(x) for x in metrics(tot))
    return {"truth": int(tot[0] + tot[1]), "tp": int(tot[0]), "fn": int(tot[1]), "tp_call": int(tot[2]), "fp": int(tot[3]),
            "calls": int(tot[2] + tot[3]), "precision": pr, "recall": rc, "f1": f1}


rows, J = [], {}
for kind in ("SNV", "INDEL"):
    for when in ("before", "after"):
        for t, p in SETS:
            for st in SNAME:
                for c in CALLERS[when]:
                    x = cell(kind, when, t, p, c, st)
                    d = None
                    if c in ("ClairS-TO", "DeepSomatic-TO"):
                        d = paired(res[(kind, when, t, p, REF[when])][st], res[(kind, when, t, p, c)][st])
                    J[f"{kind}|{when}|{ds(t, p)}|{st}|{c}"] = {**x, "pansoma_minus_caller": d, "pansoma": REF[when] if d else None}
                    rows.append({"kind": kind, "filter": when, "tumor": t, "platform": lab(t, p).split(" ", 1)[1], "stratum": SNAME[st], "caller": c,
                                 **{k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in x.items()},
                                 "pansoma_minus_caller_dF1": f"{d['d']:.4f}" if d else "", "dF1_lo": f"{d['lo']:.4f}" if d else "",
                                 "dF1_hi": f"{d['hi']:.4f}" if d else "", "pansoma_model": REF[when] if d else ""})
with open(HERE / "comparators_regions.tsv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    wr.writeheader()
    wr.writerows(rows)
(HERE / "comparators_regions.json").write_text(json.dumps({"cells": J, "dirs": dirs}, indent=1) + "\n")

WHEN = {"before": "before the population filter", "after": "after the population filter"}
L = ["# Comparators by region, before and after population filtering (chr1)", "",
     "Generated by `comparators_regions.py`; call sets and rules in its header. chr1, benchmark regions, PASS calls. Before the filter: ClairS-TO "
     "rerun without its panel-of-normals step, DeepSomatic-TO with its PON filter lifted (it still takes population allele frequencies as an "
     "input), Pansoma six-tumor models without site features, not filtered. After: every caller after the same filter with our resource VCFs; "
     "comparators in their default configuration. Pansoma release = the six-tumor HPRC v2.1 models with site features (seed 0).", ""]
for kind, cols in (("SNV", SNV_COLS), ("INDEL", [HP7, "easy", "difficult", "extreme"])):
    for when in ("before", "after"):
        L += [f"## chr1 {kind} F1, {WHEN[when]}", "", "| chr1 dataset | caller | " + " | ".join(SNAME[s] for s in cols) + " | all benchmark regions |", "|---|---|" + "---|" * (len(cols) + 1)]
        for t, p in SETS:
            L.append(f"| {lab(t, p)} | chr1 truth | " + " | ".join(f"{cell(kind, when, t, p, 'ClairS-TO', s)['truth']:,}" for s in cols + ["all"]) + " |")
            for c in CALLERS[when]:
                L.append(f"| {lab(t, p)} | {c} | " + " | ".join(f3(cell(kind, when, t, p, c, s)["f1"]) for s in cols + ["all"]) + " |")
        L += ["", f"### chr1 {kind}, {WHEN[when]}: PASS calls (precision / recall)", "", "| chr1 dataset | caller | " + " | ".join(SNAME[s] for s in cols) + " |", "|---|---|" + "---|" * len(cols)]
        for t, p in SETS:
            for c in CALLERS[when]:
                L.append(f"| {lab(t, p)} | {c} | " + " | ".join("{calls:,} ({0} / {1})".format(f3(x["precision"]), f3(x["recall"]), **x) for x in (cell(kind, when, t, p, c, s) for s in cols)) + " |")
        L += ["", f"### chr1 {kind}, {WHEN[when]}: ΔF1 = {REF[when]} − comparator [95% CI]", "", "| chr1 dataset | comparator | " + " | ".join(SNAME[s] for s in cols) + " | all benchmark regions |", "|---|---|" + "---|" * (len(cols) + 1)]
        for t, p in SETS:
            for c in ("ClairS-TO", "DeepSomatic-TO"):
                L.append(f"| {lab(t, p)} | {c} | " + " | ".join("{} [{}, {}]".format(sgn(d["d"]), sgn(d["lo"], False), sgn(d["hi"], False))
                                                                for d in (J[f"{kind}|{when}|{ds(t, p)}|{s}|{c}"]["pansoma_minus_caller"] for s in cols + ["all"])) + " |")
        L += ["", "| stratum | comparisons | ΔF1 range | Pansoma higher | interval above zero | interval below zero |", "|---|---|---|---|---|---|"]
        for s in cols + ["all"]:
            v = [J[f"{kind}|{when}|{ds(t, p)}|{s}|{c}"]["pansoma_minus_caller"] for t, p in SETS for c in ("ClairS-TO", "DeepSomatic-TO")]
            L.append(f"| {SNAME[s]} | {len(v)} | {sgn(min(x['d'] for x in v))} to {sgn(max(x['d'] for x in v))} | {sum(x['d'] > 0 for x in v)} | {sum(x['lo'] > 0 for x in v)} | {sum(x['hi'] < 0 for x in v)} |")
        L.append("")
L += ["## Checks", "", f"- {len(res)} call sets: TP / FN / TP-call / FP = report.json; one truth VCF and BED per dataset for every caller; truth per stratum identical for every caller; SMaHT tiers sum to the whole."]
(HERE / "comparators_regions.md").write_text("\n".join(L) + "\n")

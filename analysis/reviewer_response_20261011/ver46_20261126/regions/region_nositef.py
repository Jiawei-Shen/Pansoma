#!/usr/bin/env python3
"""ver4.6 item "regions": Response Tables 1.5 = 2.5 = 3.15 on the six-tumor models, three references, with NO model
taking the population site features as an input (user 2026-11-26).

Copy of ../../ver4_20261123/regions/region_six_tumor.py; only the HPRC call sets differ (function vdir), and that
file's reproduction test of older tables is not repeated here. Reads only; writes next to itself.

Call sets (R = /scratch/jshen/data/pansoma_net_v2_runs; seed 0; chr1; rtg vcfeval in the benchmark BED, eval_bed):
  SNV   GRCh38-only  R/m4m6_graphs_20261110/multicancer6_<P>_SNV_linear_nositef/vcf_chr1_<TUMOR>          (unchanged)
        HPRC v1.1    R/m6_nositef_graphs_20261118/multicancer6_<P>_SNV_d9_nositef/vcf_chr1_<SET>          (RIS, pulled 2026-11-26)
        HPRC v2.1    R/m6_nositef_graphs_20261118/multicancer6_<P>_SNV_d46_nositef/vcf_chr1_<SET>         (RIS, pulled 2026-11-26)
  INDEL GRCh38-only  R/m4m6_graphs_20261110/multicancer6_<P>_INDEL_linear_nositef/vcf_{ppr8,pp}_chr1_<TUMOR>  (unchanged)
        HPRC v1.1    R/m6_nositef_graphs_20261118/multicancer6_<P>_INDEL_d9_nositef/vcf_{ppr8,pp}_chr1_<TUMOR>
        HPRC v2.1    R/m6_nositef_graphs_20261118/multicancer6_<P>_INDEL_d46_nositef/vcf_{ppr8,pp}_chr1_<TUMOR>
        (pp for ONT, ppr8 otherwise; lab cluster)
  "before" = eval_bed/vcfeval/raw_calls (PASS calls, no population filter), "after" = pon_calls (PASS calls after the
  population filter); the truth is complete in both.
Strata, bootstrap, metric rules and the asserted checks 1-4: exactly those of region_six_tumor.py (see its header):
  SMaHT extreme v2 BED; GIAB v3.5 GRCh38 strata; INDEL column = GIAB homopolymers >= 7 bp; paired dF1 = HPRC -
  GRCh38-only over 2,000 resamples of 1-Mb blocks (seed 20261015), 95 % percentile interval; totals of every rtg
  directory = report.json; truth, BED and test-set build identical for the three references of a dataset.

usage: /wanglab/jshen/anaconda3/bin/python region_nositef.py
Writes table_before_filter.md, table_after_filter.md, whole_chr1_paired.md, whole_chr1_paired.json, summary.md,
summary.json, cells_six_tumor.json, per_callset_strata.tsv.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

T0 = time.time()
HERE = Path(__file__).resolve().parent
EXP = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments")
sys.path.insert(0, str(EXP / "E16_bootstrap"))
sys.path.insert(0, str(EXP / "E10b_giab_strata"))
import callsets as C  # noqa: E402
import e16_bootstrap as E  # noqa: E402
import giab_strata as G  # noqa: E402

SMAHT_EXTREME = "/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_extreme_v2.bed"
RESAMPLES, SEED = 2000, 20261015
GNAME = {"linear": "GRCh38-only", "d9": "HPRC v1.1", "d46": "HPRC v2.1"}
GSHORT = {"d9": "v1.1", "d46": "v2.1"}
ARMS = ("linear", "d9", "d46")
# (stratum key, letter column name, long name)
SNV_ST = [("extreme", "SMaHT extreme", "SMaHT extreme"),
          ("alldifficultregions", "GIAB all difficult", "GIAB all difficult regions"),
          ("alllowmapandsegdupregions", "GIAB low map or segdup", "GIAB low mappability or segmental duplications"),
          ("AllTandemRepeatsandHomopolymers_slop5", "GIAB tandem repeats and homopolymers", "GIAB tandem repeats and homopolymers"),
          ("notinalldifficultregions", "GIAB outside difficult", "GIAB outside all difficult regions")]
HP7 = "AllHomopolymers_ge7bp_imperfectge11bp_slop5"
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


def rule(kind, p):
    return "" if kind == "SNV" else ("pp_" if p == "ONT" else "ppr8_")


def vdir(kind, t, p, arm):
    r = rule(kind, p)
    if arm == "linear":
        return G6 / f"multicancer6_{p}_{kind}_linear_nositef" / f"vcf_{r}chr1_{t}"
    run = NS / f"multicancer6_{p}_{kind}_{arm}_nositef"
    if kind == "INDEL":
        return run / f"vcf_{r}chr1_{t}"
    # SNV: one directory per test set, named after the set (HG008T_PacBio, COLO829T_fiberseq_snvaf05,
    # HCC1395T_PacBio_d9_snvaf05, HG008T_PacBio_d46, COLO829T_fiberseq_d46_snvaf05, ...)
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
    """P, R, F1 from (..., 4) = TPb, FN, TPc, FP. F1 = 0 when there is no TP and FP + FN > 0 (rtg formula otherwise)."""
    tot = np.asarray(tot, dtype=float)
    tpb, fn, tpc, fp = tot[..., 0], tot[..., 1], tot[..., 2], tot[..., 3]
    with np.errstate(invalid="ignore", divide="ignore"):
        p = tpc / (tpc + fp)
        r = tpb / (tpb + fn)
        f = 2 * p * r / (p + r)
    f = np.where(np.isnan(f) & (tpb + tpc == 0) & (fn + fp > 0), 0.0, f)
    return p, r, f


def nan_f1(tot):
    """Number of entries where the rtg formula alone gives no F1 (E16's metrics)."""
    return int(np.isnan(E.metrics(tot)[2]).sum())


ivs = {"extreme": bed_chr1(SMAHT_EXTREME)}
ivs.update({name: G.chr1_intervals(G.GIAB / rel) for name, rel, _ in G.STRATA
            if name in [s[0] for s in SNV_ST] or name == HP7})
assert set(ivs) == {s[0] for s in SNV_ST} | {HP7}
bs = E.BlockSets(C.fasta_path(), RESAMPLES, SEED)
print("strata + block geometry %.0fs" % (time.time() - T0), file=sys.stderr)


def score(d, which, check_report=True):
    """{stratum: (totals[4], resampled totals[RESAMPLES, 4], block key)} of one rtg directory, + meta."""
    d = Path(d)
    rep = json.loads((d / "eval_bed/report.json").read_text())
    w = d / "eval_bed/vcfeval" / which
    assert C._rtg_complete(w), w
    recs = {p: C.read_records(w / f"{p}.vcf.gz") for p in PARTS}
    assert all(r[0] == "chr1" for p in PARTS for r in recs[p]), w  # chr1 call sets only
    rec = {p: np.array([r[1] for r in recs[p]], dtype=np.int64) for p in PARTS}
    stage = rep["raw" if which == "raw_calls" else "pon"]
    a = stage["at_threshold"]
    n = (len(rec["tp-baseline"]), len(rec["fn"]), len(rec["tp"]), len(rec["fp"]))
    if check_report:
        assert n == (a["tp_baseline"], a["fn"], a["tp_call"], a["fp"]), (w, n, a)
        # "baseline" = truth records inside the benchmark BED (rtg); truth_records = chr1 truth records before the BED
        assert n[0] + n[1] == stage["baseline"] <= rep["truth_records"], (w, n, stage["baseline"])
    blocks = bs.get("bed", rep["bed"], G.WIDTH)
    out = {}
    for name, iv in list(ivs.items()) + [("all", None)]:
        sub = {p: (rec[p] if iv is None else rec[p][G.inside(rec[p], iv)]).tolist() for p in PARTS}
        m = E.block_counts(sub, blocks)
        assert m.sum() == sum(len(v) for v in sub.values()), (w, name)  # no record outside the block set
        out[name] = (m.sum(axis=0), blocks.mult @ m, blocks.key)
    if check_report:
        p, r, f = (float(x) for x in metrics(out["all"][0]))
        for mine, key in ((p, "precision"), (r, "recall"), (f, "f1")):
            assert abs(mine - a[key]) < 5.1e-5, (w, key, mine, a[key])
    stats = sorted(x.name for x in d.glob("*.graph.stats.json"))
    meta = {"dir": str(d), "bed": rep["bed"], "truth": rep["truth"], "truth_records": rep["truth_records"], "baseline": stage["baseline"],
            "threshold": rep.get("threshold"), "call": rep.get("steps", {}).get("graph_vcf", {}).get("call"),
            "test_set": stats[0].split(".")[0] if len(stats) == 1 else stats,
            "report": {k: a[k] for k in ("tp_baseline", "tp_call", "fp", "fn", "precision", "recall", "f1")}}
    return out, meta


def paired(a, b):
    """dF1 a - b: point, lo, hi, share of resamples > 0, resamples dropped (F1 undefined in either arm)."""
    assert a[2] == b[2]
    d = metrics(a[1])[2] - metrics(b[1])[2]
    lo, hi, _, nnan = E.ci(d)
    ok = d[~np.isnan(d)]
    return {"d": float(metrics(a[0])[2] - metrics(b[0])[2]), "lo": lo, "hi": hi,
            "p_gt0": float((ok > 0).mean()), "dropped": nnan, "zero_tp_resamples": nan_f1(a[1]) + nan_f1(b[1]) - 0}


def f3(x):
    return "-" if x != x else f"{x:.3f}"


def sgn(x, plus=True):
    """Letter style: +0.072 / −0.268 (Unicode minus); bounds without '+'."""
    s = f"{abs(x):.3f}"
    return (MINUS if x < 0 else ("+" if plus else "")) + s


def ci_cell(c):
    return f"{sgn(c['d'])} [{sgn(c['lo'], False)}, {sgn(c['hi'], False)}]"


def prf(tot):
    return " / ".join(f3(float(v)) for v in metrics(tot))


def calls(tot):
    return int(tot[2] + tot[3])


# ---------------------------------------------------------------------------------------------------------------------
res, meta = {}, {}
for kind in ("SNV", "INDEL"):
    for t, p in SETS:
        for arm in ARMS:
            for w, _, _ in WHICH:
                res[(kind, t, p, arm, w)], meta[(kind, t, p, arm, w)] = score(vdir(kind, t, p, arm), w)
        # checks 2 and 3: same truth, BED, test-set build and truth count per stratum in the three arms, before and after
        ms = [meta[(kind, t, p, arm, w)] for arm in ARMS for w, _, _ in WHICH]
        assert len({m["bed"] for m in ms}) == 1 and len({m["truth"] for m in ms}) == 1, (kind, t, p)
        assert len({m["truth_records"] for m in ms}) == 1 and len({m["baseline"] for m in ms}) == 1, (kind, t, p)
        builds = {arm: ("snvaf05" in str(meta[(kind, t, p, arm, "raw_calls")]["test_set"])) for arm in ARMS}
        assert len(set(builds.values())) == 1, (kind, t, p, builds)
        for st in list(ivs) + ["all"]:
            truth = {int(res[(kind, t, p, arm, w)][st][0][:2].sum()) for arm in ARMS for w, _, _ in WHICH}
            assert len(truth) == 1, (kind, t, p, st, truth)
            assert len({res[(kind, t, p, arm, w)][st][2] for arm in ARMS for w, _, _ in WHICH}) == 1
    print(kind, "scored %.0fs" % (time.time() - T0), file=sys.stderr)
assert all(b.extra == 0 for b in bs.cache.values())
nblocks = {k.split("|")[1]: len(b.ids) for k, b in bs.cache.items()}


def cell(kind, t, p, arm, w, st):
    tot = res[(kind, t, p, arm, w)][st][0]
    pr, rc, f1 = (float(x) for x in metrics(tot))
    return {"tp_baseline": int(tot[0]), "fn": int(tot[1]), "tp_call": int(tot[2]), "fp": int(tot[3]),
            "truth": int(tot[0] + tot[1]), "calls": calls(tot), "precision": pr, "recall": rc, "f1": f1}


def pair(kind, t, p, arm, w, st, base="linear"):
    return paired(res[(kind, t, p, arm, w)][st], res[(kind, t, p, base, w)][st])


# ---- 1 + 2: the letter's nine-column table, before and after the filter ---------------------------------------------
HEAD = ("| chr1 dataset | HPRC graph | SMaHT extreme PASS calls (precision) | SMaHT extreme ΔF1 | GIAB all difficult ΔF1 | "
        "GIAB low map or segdup ΔF1 | GIAB tandem repeats and homopolymers ΔF1 | GIAB outside difficult ΔF1 | "
        "INDEL F1 homopolymers ≥ 7 bp |")
cells = {}
tables = {}
for w, _, when in WHICH:
    L = [HEAD, "|---" * 9 + "|"]
    for t, p in SETS:
        for arm in ("d9", "d46"):
            c0, c1 = cell("SNV", t, p, "linear", w, "extreme"), cell("SNV", t, p, arm, w, "extreme")
            row = {"extreme_calls": f"{c0['calls']:,} ({f3(c0['precision'])}) → {c1['calls']:,} ({f3(c1['precision'])})"}
            num = {"extreme_linear": c0, "extreme_hprc": c1}
            for st, _, _ in SNV_ST:
                num[st] = pair("SNV", t, p, arm, w, st)
                num[st]["truth"] = cell("SNV", t, p, arm, w, st)["truth"]
                num[st]["f1_linear"] = cell("SNV", t, p, "linear", w, st)["f1"]
                num[st]["f1_hprc"] = cell("SNV", t, p, arm, w, st)["f1"]
                row[st] = ci_cell(num[st])
            i0, i1 = cell("INDEL", t, p, "linear", w, HP7), cell("INDEL", t, p, arm, w, HP7)
            num["indel_hp7_linear"], num["indel_hp7_hprc"] = i0, i1
            row["indel_hp7"] = f"{f3(i0['f1'])} → {f3(i1['f1'])}"
            cells[f"{when}|{ds(t, p)}|{arm}"] = {"label": lab(t, p), "graph": GNAME[arm], "text": row, "values": num}
            L.append(f"| {lab(t, p)} | {GSHORT[arm]} | " + " | ".join(row.values()) + " |")
    tables[when] = L
NOTE = {
    "before": "BEFORE the population filter (rtg raw_calls: PASS calls, not filtered by the four resource VCFs)",
    "after": "AFTER the population filter (rtg pon_calls: PASS calls after the filter with the four resource VCFs)"}
COMMON = (
    "chr1, benchmark regions, six-tumor models, seed 0, none with the population site features as a model input "
    "(generated by region_nositef.py; do not edit by hand). "
    "SNV: multicancer6_<P>_SNV on GRCh38-only (m4m6_graphs_20261110 linear_nositef), HPRC v1.1 and HPRC v2.1 "
    "(m6_nositef_graphs_20261118/multicancer6_<P>_SNV_{d9,d46}_nositef); PASS = rule 2 (threshold at post-filter "
    "validation recall 0.95 of the callable positives). INDEL column: six-tumor INDEL models without the site features "
    "(GRCh38-only m4m6_graphs_20261110 linear_nositef, v1.1 / v2.1 m6_nositef_graphs_20261118/multicancer6_<P>_INDEL_"
    "{d9,d46}_nositef), PASS = ppr8 (Illumina, PacBio HiFi / Fiber-seq), pp (ONT). "
    "Arrows: GRCh38-only → HPRC. ΔF1 = HPRC − GRCh38-only with the 95 % percentile interval of 2,000 resamples of 1-Mb "
    "blocks (seed 20261015, same resamples for the two arms). F1 is 0.000 when a call set has no TP in the stratum.")
for when, fn in (("before", "table_before_filter.md"), ("after", "table_after_filter.md")):
    (HERE / fn).write_text(
        f"# Response Table 1.5 (= 2.5 = 3.15), six-tumor models without site features, {NOTE[when]}\n\n{COMMON} INDEL column: {when} the "
        f"filter as the SNV columns.\n\n" + "\n".join(tables[when]) + "\n")

# ---- 3: whole chr1 (all benchmark regions) ---------------------------------------------------------------------------
whole = {}
W = ["# Whole chr1 (benchmark regions, no stratum): six-tumor models on three references, paired ΔF1", "",
     "Generated by region_nositef.py. P / R / F1 of the PASS calls (rtg vcfeval in the benchmark BED; the point "
     "values are asserted equal to eval_bed/report.json at_threshold: counts exactly, P / R / F1 to 4 decimals). "
     "ΔF1 = HPRC − GRCh38-only, 95 % percentile interval of 2,000 resamples of 1-Mb blocks (seed 20261015), the same "
     "resamples for the three references of a dataset. Seed 0 models. SNV PASS = rule 2; INDEL PASS = ppr8 (ONT: pp). "
     "No model takes the population site features as an input (HPRC: m6_nositef_graphs_20261118).", ""]
for kind in ("SNV", "INDEL"):
    W += [f"## {kind}", "",
          "| chr1 dataset | truth | population filter | GRCh38-only P / R / F1 | HPRC v1.1 P / R / F1 | HPRC v2.1 P / R / F1 | "
          "ΔF1 v1.1 − GRCh38-only [95 % CI] | ΔF1 v2.1 − GRCh38-only [95 % CI] |", "|---" * 8 + "|"]
    for t, p in SETS:
        for w, _, when in WHICH:
            c = {arm: cell(kind, t, p, arm, w, "all") for arm in ARMS}
            d = {arm: pair(kind, t, p, arm, w, "all") for arm in ("d9", "d46")}
            d["d46_minus_d9"] = pair(kind, t, p, "d46", w, "all", base="d9")
            whole[f"{kind}|{ds(t, p)}|{when}"] = {
                "label": lab(t, p), "truth": c["linear"]["truth"], "blocks": nblocks[meta[(kind, t, p, "linear", w)]["bed"]],
                "arms": {GNAME[a]: c[a] for a in ARMS},
                "dF1_v1.1_minus_GRCh38-only": d["d9"], "dF1_v2.1_minus_GRCh38-only": d["d46"],
                "dF1_v2.1_minus_v1.1": d["d46_minus_d9"],
                "dirs": {GNAME[a]: meta[(kind, t, p, a, w)]["dir"] for a in ARMS},
                "test_sets": {GNAME[a]: meta[(kind, t, p, a, w)]["test_set"] for a in ARMS},
                "thresholds": {GNAME[a]: meta[(kind, t, p, a, w)]["threshold"] for a in ARMS}}
            W.append(f"| {lab(t, p)} | {c['linear']['truth']:,} | {when} | "
                     + " | ".join(" / ".join(f3(c[a][k]) for k in ("precision", "recall", "f1")) for a in ARMS)
                     + f" | {ci_cell(d['d9'])} | {ci_cell(d['d46'])} |")
    W.append("")

# ---- 4: summary for the text -----------------------------------------------------------------------------------------
S = ["# Summary lines for R1.3 / R2.5 (six-tumor models without site features; generated by region_nositef.py)", "",
     "18 pairs = 9 datasets (HG008T, COLO829T, HCC1395T × three platforms) × HPRC v1.1 / v2.1, each against GRCh38-only. "
     "\"interval above zero\" = lower bound of the 95 % interval > 0; \"below zero\" = upper bound < 0.", ""]
summ = {"snv_strata": {}, "whole_chr1": {}, "indel_hp7": {}, "smaht_extreme": {}, "fp_share": {}}


def pair_name(k):
    _, s, arm = k.split("|")
    t, p = s.split("_")
    return f"{lab(t, 'PacBio' if p in ('PacBio', 'fiberseq') else p)} {GSHORT[arm]}"


def tally(v):
    """v: {pair key: paired dict} -> counts and lists."""
    pos = [k for k, x in v.items() if x["d"] > 0]
    above = [k for k, x in v.items() if x["lo"] > 0]
    below = [k for k, x in v.items() if x["hi"] < 0]
    incl = [k for k, x in v.items() if x["lo"] <= 0 <= x["hi"]]
    assert len(above) + len(below) + len(incl) == len(v)
    return {"n": len(v), "positive": len(pos), "not_positive": [pair_name(k) for k in v if k not in pos],
            "interval_above_zero": len(above), "interval_below_zero": len(below),
            "below_zero_pairs": [pair_name(k) for k in below], "includes_zero_pairs": [pair_name(k) for k in incl],
            "min": min(x["d"] for x in v.values()), "max": max(x["d"] for x in v.values()),
            "min_pair": pair_name(min(v, key=lambda k: v[k]["d"])), "max_pair": pair_name(max(v, key=lambda k: v[k]["d"]))}


for when in ("before", "after"):
    S += [f"## SNV ΔF1 (HPRC − GRCh38-only) per stratum, {when} the population filter", "",
          "| stratum | pairs | ΔF1 range | ΔF1 > 0 | interval above zero | interval below zero | interval includes zero |",
          "|---|---|---|---|---|---|---|"]
    for scope, arms in (("18 pairs", ("d9", "d46")), ("v1.1 (9)", ("d9",)), ("v2.1 (9)", ("d46",))):
        for st, short, _ in SNV_ST + [("all", "whole chr1 (benchmark regions)", "")]:
            if st == "all":
                v = {f"{when}|{ds(t, p)}|{arm}": pair("SNV", t, p, arm, WHICH[0 if when == "before" else 1][0], "all")
                     for t, p in SETS for arm in arms}
            else:
                v = {k: c["values"][st] for k, c in cells.items() if k.startswith(when + "|") and k.split("|")[2] in arms}
            ty = tally(v)
            summ["snv_strata"][f"{when}|{scope}|{st}"] = ty
            S.append(f"| {short} | {scope} | {sgn(ty['min'])} to {sgn(ty['max'])} | {ty['positive']} of {ty['n']} | "
                     f"{ty['interval_above_zero']} of {ty['n']} | {ty['interval_below_zero']}"
                     + (f" ({'; '.join(ty['below_zero_pairs'])})" if ty["below_zero_pairs"] else "")
                     + f" | {len(ty['includes_zero_pairs'])}"
                     + (f" ({'; '.join(ty['includes_zero_pairs'])})" if ty["includes_zero_pairs"] else "") + " |")
    S.append("")
for when, w in (("before", "raw_calls"), ("after", "pon_calls")):
    S += [f"## SMaHT extreme regions, SNV PASS calls (precision) and FP per reference, {when} the population filter", "",
          "| chr1 dataset | truth SNVs in SMaHT extreme | GRCh38-only calls (precision) | HPRC v1.1 calls (precision) | "
          "HPRC v2.1 calls (precision) | FP GRCh38-only / v1.1 / v2.1 | F1 GRCh38-only / v1.1 / v2.1 |", "|---" * 7 + "|"]
    for t, p in SETS:
        c = {arm: cell("SNV", t, p, arm, w, "extreme") for arm in ARMS}
        summ["smaht_extreme"][f"{when}|{ds(t, p)}"] = {GNAME[a]: c[a] for a in ARMS}
        S.append(f"| {lab(t, p)} | {c['linear']['truth']} | "
                 + " | ".join(f"{c[a]['calls']:,} ({f3(c[a]['precision'])})" for a in ARMS)
                 + " | " + " / ".join(f"{c[a]['fp']:,}" for a in ARMS)
                 + " | " + " / ".join(f3(c[a]["f1"]) for a in ARMS) + " |")
    n_fp = {arm: sum(summ["smaht_extreme"][f"{when}|{ds(t, p)}"][GNAME[arm]]["fp"]
                     < summ["smaht_extreme"][f"{when}|{ds(t, p)}"]["GRCh38-only"]["fp"] for t, p in SETS) for arm in ("d9", "d46")}
    n_pr = {arm: sum(summ["smaht_extreme"][f"{when}|{ds(t, p)}"][GNAME[arm]]["precision"]
                     > summ["smaht_extreme"][f"{when}|{ds(t, p)}"]["GRCh38-only"]["precision"] for t, p in SETS) for arm in ("d9", "d46")}
    n_f1 = {arm: sum(summ["smaht_extreme"][f"{when}|{ds(t, p)}"][GNAME[arm]]["f1"]
                     > summ["smaht_extreme"][f"{when}|{ds(t, p)}"]["GRCh38-only"]["f1"] for t, p in SETS) for arm in ("d9", "d46")}
    summ["smaht_extreme"][f"{when}|counts"] = {"fewer_fp": n_fp, "higher_precision": n_pr, "higher_f1": n_f1}
    S += ["", f"Datasets (of 9) with fewer SMaHT-extreme FP than GRCh38-only: v1.1 {n_fp['d9']}, v2.1 {n_fp['d46']}; "
          f"higher precision: v1.1 {n_pr['d9']}, v2.1 {n_pr['d46']}; higher F1: v1.1 {n_f1['d9']}, v2.1 {n_f1['d46']}.", ""]
for when, w in (("before", "raw_calls"), ("after", "pon_calls")):
    S += [f"## Share of SNV false positives that lie in SMaHT extreme regions, {when} the population filter", "",
          "| chr1 dataset | GRCh38-only FP in extreme / all (share) | HPRC v1.1 | HPRC v2.1 | "
          "FP of GRCh38-only minus FP of v2.1: in extreme | outside extreme | all |", "|---" * 7 + "|"]
    for t, p in SETS:
        e = {arm: cell("SNV", t, p, arm, w, "extreme")["fp"] for arm in ARMS}
        a = {arm: cell("SNV", t, p, arm, w, "all")["fp"] for arm in ARMS}
        ex_e, ex_a = e["linear"] - e["d46"], a["linear"] - a["d46"]
        summ["fp_share"][f"{when}|{ds(t, p)}"] = {"fp_extreme": {GNAME[x]: e[x] for x in ARMS}, "fp_all": {GNAME[x]: a[x] for x in ARMS},
                                                 "excess_extreme": ex_e, "excess_all": ex_a}
        S.append(f"| {lab(t, p)} | " + " | ".join(f"{e[x]:,} / {a[x]:,} ({100 * e[x] / a[x]:.1f} %)" for x in ARMS)
                 + " | " + " | ".join(f"{x:+,}".replace("-", MINUS) for x in (ex_e, ex_a - ex_e, ex_a)) + " |")
    S.append("")
for when, w in (("before", "raw_calls"), ("after", "pon_calls")):
    S += [f"## INDEL in homopolymers ≥ 7 bp, P / R / F1 per reference, {when} the population filter", "",
          "| chr1 dataset | truth INDELs in HP ≥ 7 bp | GRCh38-only P / R / F1 (calls) | HPRC v1.1 P / R / F1 (calls) | "
          "HPRC v2.1 P / R / F1 (calls) |", "|---" * 5 + "|"]
    hi = {"d9": [], "d46": []}
    hip = {"d9": [], "d46": []}
    for t, p in SETS:
        c = {arm: cell("INDEL", t, p, arm, w, HP7) for arm in ARMS}
        summ["indel_hp7"][f"{when}|{ds(t, p)}"] = {GNAME[a]: c[a] for a in ARMS}
        for arm in ("d9", "d46"):
            if c[arm]["f1"] > c["linear"]["f1"]:
                hi[arm].append(lab(t, p))
            if c[arm]["precision"] > c["linear"]["precision"]:
                hip[arm].append(lab(t, p))
        S.append(f"| {lab(t, p)} | {c['linear']['truth']} | "
                 + " | ".join(" / ".join(f3(c[a][k]) for k in ("precision", "recall", "f1")) + f" ({c[a]['calls']:,})" for a in ARMS) + " |")
    lower = [f"{lab(t, p)} {GSHORT[arm]}" for t, p in SETS for arm in ("d9", "d46") if lab(t, p) not in hi[arm]]
    summ["indel_hp7"][f"{when}|counts"] = {"f1_higher_on_hprc": len(hi["d9"]) + len(hi["d46"]), "of": 18,
                                           "v1.1": len(hi["d9"]), "v2.1": len(hi["d46"]), "not_higher": lower,
                                           "precision_higher_on_hprc": len(hip["d9"]) + len(hip["d46"]),
                                           "precision_higher_v1.1": len(hip["d9"]), "precision_higher_v2.1": len(hip["d46"])}
    S += ["", f"INDEL F1 in HP ≥ 7 bp higher on HPRC than on GRCh38-only in {len(hi['d9']) + len(hi['d46'])} of 18 pairs "
          f"(v1.1 {len(hi['d9'])} of 9, v2.1 {len(hi['d46'])} of 9); not higher: {'; '.join(lower) or 'none'}. "
          f"Precision higher on HPRC in {len(hip['d9']) + len(hip['d46'])} of 18 pairs (v1.1 {len(hip['d9'])}, v2.1 {len(hip['d46'])}).", ""]
for kind in ("SNV", "INDEL"):
    S += [f"## Whole chr1, {kind}: ΔF1 (HPRC − GRCh38-only) over the pairs", "",
          "| population filter | pairs | ΔF1 range | ΔF1 > 0 | interval above zero | interval below zero | interval includes zero |",
          "|---|---|---|---|---|---|---|"]
    for when in ("before", "after"):
        for scope, keys in (("18 pairs", ("dF1_v1.1_minus_GRCh38-only", "dF1_v2.1_minus_GRCh38-only")),
                            ("v1.1 (9)", ("dF1_v1.1_minus_GRCh38-only",)), ("v2.1 (9)", ("dF1_v2.1_minus_GRCh38-only",))):
            v = {f"{when}|{ds(t, p)}|{'d9' if 'v1.1' in k else 'd46'}": whole[f"{kind}|{ds(t, p)}|{when}"][k] for t, p in SETS for k in keys}
            ty = tally(v)
            summ["whole_chr1"][f"{kind}|{when}|{scope}"] = ty
            S.append(f"| {when} | {scope} | {sgn(ty['min'])} to {sgn(ty['max'])} | {ty['positive']} of {ty['n']} | "
                     f"{ty['interval_above_zero']} of {ty['n']} | {ty['interval_below_zero']}"
                     + (f" ({'; '.join(ty['below_zero_pairs'])})" if ty["below_zero_pairs"] else "")
                     + f" | {len(ty['includes_zero_pairs'])}"
                     + (f" ({'; '.join(ty['includes_zero_pairs'])})" if ty["includes_zero_pairs"] else "") + " |")
    S.append("")
# truth counts and bootstrap bookkeeping
S += ["## Truth counts per dataset (identical on the three references, before and after the filter)", "",
      "| chr1 tumor | SNV truth | in SMaHT extreme | in GIAB all difficult | in low map or segdup | in tandem repeats and homopolymers | "
      "outside difficult | INDEL truth | INDEL in HP ≥ 7 bp | 1-Mb blocks |", "|---" * 10 + "|"]
truth_counts = {}
for t in ("HG008T", "COLO829T", "HCC1395T"):
    n = {p: [cell("SNV", t, p, "linear", "raw_calls", st)["truth"] for st in ["all"] + [s[0] for s in SNV_ST]]
         + [cell("INDEL", t, p, "linear", "raw_calls", "all")["truth"], cell("INDEL", t, p, "linear", "raw_calls", HP7)["truth"]]
         for p in ("Illumina", "PacBio", "ONT")}
    assert n["Illumina"] == n["PacBio"] == n["ONT"], (t, n)  # one truth per tumor
    v = n["Illumina"]
    truth_counts[t] = dict(zip(["snv", "snv_extreme"] + [s[0] for s in SNV_ST[1:]] + ["indel", "indel_hp7"], v))
    truth_counts[t]["blocks"] = nblocks[meta[("SNV", t, "Illumina", "linear", "raw_calls")]["bed"]]
    S.append(f"| {t} | {v[0]:,} | {v[1]:,} | {v[2]:,} | {v[3]:,} | {v[4]:,} | {v[5]:,} | {v[6]:,} | {v[7]:,} | {truth_counts[t]['blocks']} |")
dropped = max(c["values"][st]["dropped"] for c in cells.values() for st, _, _ in SNV_ST)
zero_tp = max(c["values"][st]["zero_tp_resamples"] for c in cells.values() for st, _, _ in SNV_ST)
wd = max(v[k]["dropped"] for v in whole.values() for k in v if k.startswith("dF1"))
wz = max(v[k]["zero_tp_resamples"] for v in whole.values() for k in v if k.startswith("dF1"))
S += ["", f"Bootstrap bookkeeping: resamples dropped because F1 was undefined in an arm: at most {dropped} of {RESAMPLES} in a "
      f"table cell and {wd} in a whole-chr1 cell; resamples in which the F1 = 0 rule (no TP) applied: at most {zero_tp} in a "
      f"table cell and {wz} in a whole-chr1 cell (0 = identical to E16's bootstrap).", ""]
zcells = [f"{k} {kk} ({v[kk]['zero_tp_resamples']})" for k, v in whole.items() for kk in v
          if kk.startswith("dF1") and v[kk]["zero_tp_resamples"] > 0]
S += ["Whole-chr1 cells in which the F1 = 0 rule applied in some resamples (resamples): " + ("; ".join(zcells) or "none")
      + ". In these cells E16's rule would drop the resamples instead.", ""]
summ["truth_counts"] = truth_counts
summ["bootstrap"] = {"resamples": RESAMPLES, "seed": SEED, "block_bp": int(G.WIDTH), "max_dropped_table": dropped,
                     "max_zero_tp_table": zero_tp, "max_dropped_whole": wd, "max_zero_tp_whole": wz, "blocks": nblocks}

# per call set and stratum (long table)
T = ["kind\tdataset\treference\tfilter\tstratum\ttruth\ttp_baseline\tfn\ttp_call\tfp\tcalls\tprecision\trecall\tf1"]
for (kind, t, p, arm, w) in res:
    for st in list(ivs) + ["all"]:
        c = cell(kind, t, p, arm, w, st)
        T.append("\t".join(str(x) for x in (
            kind, ds(t, p), GNAME[arm], "before" if w == "raw_calls" else "after", st, c["truth"], c["tp_baseline"], c["fn"],
            c["tp_call"], c["fp"], c["calls"], f3(c["precision"]), f3(c["recall"]), f3(c["f1"]))))

(HERE / "whole_chr1_paired.md").write_text("\n".join(W) + "\n")
(HERE / "whole_chr1_paired.json").write_text(json.dumps(whole, indent=1))
(HERE / "summary.md").write_text("\n".join(S) + "\n")
(HERE / "summary.json").write_text(json.dumps(summ, indent=1))
(HERE / "per_callset_strata.tsv").write_text("\n".join(T) + "\n")
(HERE / "cells_six_tumor.json").write_text(json.dumps({
    "note": "region_nositef.py; key = <before|after>|<dataset>|<d9|d46>; text = the table cells; values = numbers",
    "columns": ["extreme_calls"] + [s[0] for s in SNV_ST] + ["indel_hp7"],
    "cells": cells,
    "callsets": {"|".join(k): v for k, v in meta.items()}}, indent=1))
print("checks passed: %d rtg directories = report.json (counts, P / R / F1); truth counts equal across references in "
      "every stratum; no record outside the blocks" % len(res))
print("total %.0fs" % (time.time() - T0), file=sys.stderr)

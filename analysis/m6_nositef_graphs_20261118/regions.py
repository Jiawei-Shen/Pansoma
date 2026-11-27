#!/usr/bin/env python3
"""regions.md: the SNV columns of Response Tables 1.5 = 2.5 = 3.15 (chr1 SNV calling by region, HPRC against
GRCh38-only) with the HPRC side on the six-tumor models of this batch, so that neither side takes the population site
features as an input (user 2026-11-24). Point estimates only: see "Not here".

Same rules as ../reviewer_response_20261011/ver4_20261123/regions/region_six_tumor.py (the source of the tables):
rtg vcfeval records in the benchmark BED, "before" = eval_bed/vcfeval/raw_calls (PASS calls, no population filter),
"after" = pon_calls; a record counts in a stratum when its VCF POS - 1 lies in one of the stratum's chr1 intervals
(strata overlap); P = TP calls / (TP calls + FP), R = TP / (TP + FN), F1 = 0 when there is no TP but there are FP or FN.

  HPRC v1.1 / v2.1  $PANSOMA_RUNS/m6_nositef_graphs_20261118/multicancer6_<P>_SNV_{d9,d46}_nositef/vcf_chr1_<SET>
                    (stratified here)
  GRCh38-only       m4m6_graphs_20261110/multicancer6_<P>_SNV_linear_nositef: its call sets are on the old cluster;
                    its counts per stratum are read from ver4_20261123/regions/per_callset_strata.tsv (unchanged), as
                    are those of the HPRC models with the site features (for the comparison table)
Strata: GIAB genome stratifications v3.5 GRCh38 (Dataset/GIAB_stratifications_v3.5, downloaded 2026-11-24 from the
NIST FTP, md5 equal to genome-stratifications-md5s.txt) and the SMaHT extreme regions (Dataset/
SMaHT_Regional_Categorization/SMaHT_extreme_hg38.bed.gz, github.com/parklab/SMaHT_Regional_Categorization eef108c;
the old cluster's SMaHT_extreme_v2.bed is not on storage2; chr1 non-N bases of this file = 25,948,324 = the value
recorded for the v2 file in E11b_simulated_alignment/README.md).

Checks (asserted): totals of every raw_calls / pon_calls directory = report.json; the truth records per stratum equal
those of the GRCh38-only rows of per_callset_strata.tsv for every dataset, stratum, before and after (the truth and
the strata are the same as the old cluster's).

Not here: the paired block-bootstrap intervals of the tables (they resample the GRCh38-only and HPRC records together;
the GRCh38-only records are on the old cluster), and the INDEL column (no INDEL model without the site features).
usage: /home/jiawei.s/miniconda3/envs/convnextv2/bin/python regions.py   (login node; needs numpy)
"""
import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import collect as C  # noqa: E402

DATA = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset")
GIAB = DATA / "GIAB_stratifications_v3.5/GRCh38@all"
OLD = HERE.parent / "reviewer_response_20261011/ver4_20261123/regions/per_callset_strata.tsv"
STRATA = [("extreme", DATA / "SMaHT_Regional_Categorization/SMaHT_extreme_hg38.bed.gz", "SMaHT extreme"),
          ("alldifficultregions", GIAB / "Union/GRCh38_alldifficultregions.bed.gz", "GIAB all difficult"),
          ("alllowmapandsegdupregions", GIAB / "Union/GRCh38_alllowmapandsegdupregions.bed.gz", "GIAB low map or segdup"),
          ("AllTandemRepeatsandHomopolymers_slop5", GIAB / "LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
           "GIAB tandem repeats and homopolymers"),
          ("notinalldifficultregions", GIAB / "Union/GRCh38_notinalldifficultregions.bed.gz", "GIAB outside difficult")]
KEYS = [k for k, _, _ in STRATA]
PARTS = ("tp-baseline", "fn", "tp", "fp")
WHICH = (("raw_calls", "raw", "before"), ("pon_calls", "pon", "after"))
GSHORT = {"HPRC v1.1": "v1.1", "HPRC v2.1": "v2.1"}


def chr1_intervals(path):
    """Merged chr1 intervals (starts, ends) of a bed.gz (giab_strata.chr1_intervals)."""
    iv = []
    with gzip.open(path, "rt") as f:
        for line in f:
            c = line.split("\t", 3)
            if c[0] == "chr1":
                iv.append((int(c[1]), int(c[2])))
    iv.sort()
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    a = np.array(merged, dtype=np.int64).reshape(-1, 2)
    return a[:, 0], a[:, 1]


def inside(pos1, iv):
    starts, ends = iv
    p0 = np.asarray(pos1, dtype=np.int64) - 1
    i = np.searchsorted(starts, p0, side="right") - 1
    ok = i >= 0
    out = np.zeros(p0.shape, dtype=bool)
    out[ok] = p0[ok] < ends[i[ok]]
    return out


def positions(path):
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] != "#":
                c = line.split("\t", 2)
                assert c[0] == "chr1", path
                out.append(int(c[1]))
    return np.array(out, dtype=np.int64)


def prf(x):
    """dict of counts -> P, R, F1 (F1 0 when no TP but FP or FN)."""
    p = x["tp"] / (x["tp"] + x["fp"]) if x["tp"] + x["fp"] else None
    r = x["tp-baseline"] / (x["tp-baseline"] + x["fn"]) if x["tp-baseline"] + x["fn"] else None
    if p is None and r is None:
        return p, r, None
    f = 2 * p * r / (p + r) if p and r else 0.0
    return p, r, f


def old_counts():
    """{(dataset, reference, before|after, stratum): counts} of the SNV rows of per_callset_strata.tsv."""
    out = {}
    with open(OLD) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["kind"] == "SNV":
                out[r["dataset"], r["reference"], r["filter"], r["stratum"]] = {
                    "tp-baseline": int(r["tp_baseline"]), "fn": int(r["fn"]), "tp": int(r["tp_call"]), "fp": int(r["fp"])}
    return out


def f3(x, sign=False):
    return "–" if x is None else (f"{x:+.3f}" if sign else f"{x:.3f}").replace("-", "−")


def main():
    ivs = {k: chr1_intervals(p) for k, p, _ in STRATA}
    old = old_counts()
    new = {}
    for t in C.TUMORS:
        for p in C.PLATFORMS:
            dsn = "COLO829T_fiberseq" if (t, p) == ("COLO829T", "PacBio") else f"{t}_{p}"
            for g, gl in C.GRAPHS:
                v = C.RUNS / f"multicancer6_{p}_SNV_{g}_nositef" / f"vcf_chr1_{C.test_set(t, p, g)}" / "eval_bed"
                rep = json.loads((v / "report.json").read_text())
                for d, key, when in WHICH:
                    pos = {part: positions(v / "vcfeval" / d / f"{part}.vcf.gz") for part in PARTS}
                    a = rep[key]["at_threshold"]
                    assert [len(pos[x]) for x in PARTS] == [a["tp_baseline"], a["fn"], a["tp_call"], a["fp"]], (v, d)
                    new[dsn, gl, when, "all"] = {x: len(pos[x]) for x in PARTS}
                    for k in KEYS:
                        c = {x: int(inside(pos[x], ivs[k]).sum()) for x in PARTS}
                        o = old[dsn, "GRCh38-only", when, k]
                        assert c["tp-baseline"] + c["fn"] == o["tp-baseline"] + o["fn"], (dsn, gl, when, k)
                        new[dsn, gl, when, k] = c
    order = [("COLO829T_fiberseq" if (t, p) == ("COLO829T", "PacBio") else f"{t}_{p}", f"{t} {C.plat(t, p)}")
             for t in ("HG008T", "COLO829T", "HCC1395T") for p in C.PLATFORMS]

    def dF(dsn, gl, when, k, src):
        return prf(src[dsn, gl, when, k])[2] - prf(old[dsn, "GRCh38-only", when, k])[2]

    out = ["# Response Tables 1.5 = 2.5 = 3.15, SNV columns, without the population site features in any model", "",
           "Generated by `regions.py`. chr1, benchmark regions; six-tumor models, seed 0; GRCh38-only = "
           "`multicancer6_<P>_SNV_linear_nositef` (counts of `ver4_20261123/regions/per_callset_strata.tsv`), HPRC v1.1 / "
           "v2.1 = `m6_nositef_graphs_20261118/multicancer6_<P>_SNV_{d9,d46}_nositef` (stratified here). ΔF1 = HPRC − "
           "GRCh38-only, point estimate; the paired bootstrap intervals of the letter's table are not recomputed (the "
           "GRCh38-only records are on the old cluster). The INDEL column of the table is not here (no INDEL model "
           "without the site features).", ""]
    summ = {}
    for _, _, when in WHICH:
        out += [f"## {when.capitalize()} the population filter" + (" (the letter's table)" if when == "before" else ""), "",
                "| chr1 dataset | HPRC graph | SMaHT extreme PASS calls (precision) | " + " | ".join(f"{n} ΔF1" for _, _, n in STRATA)
                + " | whole chr1 ΔF1 |", "|---|---|---|" + "---|" * (len(STRATA) + 1)]
        for dsn, label in order:
            for _, gl in C.GRAPHS:
                lo, hi = old[dsn, "GRCh38-only", when, "extreme"], new[dsn, gl, when, "extreme"]
                cell = (f"{lo['tp'] + lo['fp']:,} ({f3(prf(lo)[0])}) → {hi['tp'] + hi['fp']:,} ({f3(prf(hi)[0])})")
                ds_ = [dF(dsn, gl, when, k, new) for k in KEYS + ["all"]]
                out.append(f"| {label} | {GSHORT[gl]} | {cell} | " + " | ".join(f3(x, True) for x in ds_) + " |")
                for k, x in zip(KEYS + ["all"], ds_):
                    summ.setdefault((when, k), []).append((label, GSHORT[gl], x, dF(dsn, gl, when, k, old)))
                summ.setdefault((when, "fp"), []).append((label, GSHORT[gl], hi["fp"] < lo["fp"], prf(hi)[0] > prf(lo)[0]))
        out.append("")
    out += ["## Counts over the 18 dataset–graph pairs", "",
            "| filter | stratum | ΔF1 > 0, no site features | range, no site features | ΔF1 > 0, HPRC with site features (the letter's table) | range, with site features | pairs where the no-site-feature ΔF1 is the smaller |",
            "|---|---|---|---|---|---|---|"]
    names = dict([(k, n) for k, _, n in STRATA] + [("all", "whole chr1")])
    for _, _, when in WHICH:
        for k in KEYS + ["all"]:
            s = summ[when, k]
            neg = [f"{a} {b} ({f3(x, True)})" for a, b, x, _ in s if x <= 0]
            out.append(f"| {when} | {names[k]} | {sum(x > 0 for _, _, x, _ in s)} of 18" + (f"; not: {', '.join(neg)}" if neg else "")
                       + f" | {f3(min(x for _, _, x, _ in s), True)} to {f3(max(x for _, _, x, _ in s), True)} | "
                       f"{sum(y > 0 for _, _, _, y in s)} of 18 | {f3(min(y for _, _, _, y in s), True)} to "
                       f"{f3(max(y for _, _, _, y in s), True)} | {sum(x < y for _, _, x, y in s)} of 18 |")
    out += ["", "## SMaHT extreme regions: false-positive calls and precision against GRCh38-only", ""]
    for _, _, when in WHICH:
        s = summ[when, "fp"]
        for gs in ("v1.1", "v2.1"):
            q = [x for x in s if x[1] == gs]
            out.append(f"- {when} the filter, HPRC {gs}: fewer false-positive calls in {sum(x[2] for x in q)} of 9 datasets, "
                       f"higher precision in {sum(x[3] for x in q)} of 9"
                       + "".join(f"; not fewer FP: {x[0]}" for x in q if not x[2]) + "".join(f"; not higher precision: {x[0]}" for x in q if not x[3]))
    out += ["", "## Counts per call set and stratum (HPRC without site features)", "",
            "| dataset | reference | filter | stratum | truth | TP | FN | TP calls | FP | P | R | F1 |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for dsn, label in order:
        for _, gl in C.GRAPHS:
            for _, _, when in WHICH:
                for k in KEYS + ["all"]:
                    c = new[dsn, gl, when, k]
                    out.append(f"| {label} | {gl} | {when} | {names[k]} | {c['tp-baseline'] + c['fn']} | {c['tp-baseline']} | {c['fn']} | "
                               f"{c['tp']} | {c['fp']} | " + " | ".join(f3(x) for x in prf(c)) + " |")
    out += ["", "## Checks", "", "- totals of every `raw_calls` / `pon_calls` directory equal its `report.json` (asserted)",
            "- truth records per stratum equal the GRCh38-only rows of `per_callset_strata.tsv` for all 9 datasets x 5 strata, "
            "before and after (asserted): same truth, same strata as the old cluster's"]
    (HERE / "regions.md").write_text("\n".join(out) + "\n")
    print("\n".join(out[:out.index("## Counts per call set and stratum (HPRC without site features)")]))


if __name__ == "__main__":
    main()

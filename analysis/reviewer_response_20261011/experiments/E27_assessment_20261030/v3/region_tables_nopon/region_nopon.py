#!/usr/bin/env python3
"""Supplementary Table R1.3a before population filtering (user 2026-11-10: "这个analysis 要before the PoN").

Same call sets, thresholds, blocks and strata as the published table (../region_tables/make_tables.py, which reads E16
e16_strata.py and E10b giab_strata.py), with the rtg directory switched from eval_bed/vcfeval/pon_calls (PASS calls
after our population filter, PoN-tagged truth removed) to eval_bed/vcfeval/raw_calls (the same PASS calls at the same
validation recall-0.9 threshold, no population filter, complete truth in the benchmark BED).

- Call sets: the A / B uids of E16_bootstrap/strata_paired.tsv (SNV) and E10b_giab_strata/strata_paired.tsv (SNV, INDEL).
- Strata: SMaHT v2 extreme (/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_extreme_v2.bed) and the GIAB v3.5
  strata of the table; a record is in a stratum when its VCF position lies in it (GIAB: as giab_strata.py; SMaHT:
  E10 uses the first base of the normalized interval, which is the VCF position for SNVs; checked below).
- Bootstrap: E16's 1-Mb blocks of the set's BED, 2,000 resamples, seed 20261015; paired HPRC - GRCh38-only.
- Check: the same code on pon_calls must reproduce the published cells (E16 strata_ci / strata_paired, E10b
  strata_paired / strata_ci); it asserts that before writing anything.
Writes tables.md, cells.json next to this script. usage: /wanglab/jshen/anaconda3/bin/python region_nopon.py
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

SMAHT_EXTREME = "/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_extreme_v2.bed"
SETS = [("HG008T_Illumina", "HG008T Illumina"), ("HG008T_PacBio", "HG008T PacBio HiFi"), ("HG008T_ONT", "HG008T ONT"),
        ("COLO829T_Illumina", "COLO829T Illumina"), ("COLO829T_fiberseq", "COLO829T Fiber-seq"), ("COLO829T_ONT", "COLO829T ONT")]
GNAME = {"d9": "HPRC v1.1", "d46": "HPRC v2.1"}
GIAB_ST = [("alldifficultregions", "all difficult"), ("alllowmapandsegdupregions", "low mappability or segdup"),
           ("AllTandemRepeatsandHomopolymers_slop5", "tandem repeats and homopolymers"), ("notinalldifficultregions", "not difficult")]
HP7 = "AllHomopolymers_ge7bp_imperfectge11bp_slop5"
PARTS = ("tp-baseline", "fn", "tp", "fp")


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def bed_chr1(path):
    iv = sorted((int(c[1]), int(c[2])) for c in (l.split("\t")[:3] for l in open(path)) if c[0] == "chr1")
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    a = np.array(merged, dtype=np.int64).reshape(-1, 2)
    return a[:, 0], a[:, 1]


# call sets of the published table: (set, kind, arm) -> uid
pick = {}
for r in tsv(EXP / "E16_bootstrap/strata_paired.tsv"):
    if r["stratum"] == "extreme":
        pick[(r["set"], "SNV", r["arm"])] = r["A"]; pick[(r["set"], "SNV", "linear")] = r["B"]
for r in tsv(EXP / "E10b_giab_strata/strata_paired.tsv"):
    for arm, uid in ((r["arm"], r["A"]), ("linear", r["B"])):
        assert pick.setdefault((r["set"], r["kind"], arm), uid) == uid, (r["set"], r["kind"], arm)
want = set(pick.values())
units = {u.uid: u for u in C.discover(kinds=("SNV", "INDEL"), comparators=False) if u.uid in want}
assert set(units) == want, want - set(units)

ivs = {"extreme": bed_chr1(SMAHT_EXTREME)}
ivs.update({name: G.chr1_intervals(G.GIAB / rel) for name, rel, _ in G.STRATA if name in dict(GIAB_ST) or name == HP7})
bs = E.BlockSets(C.fasta_path(list(units.values())), 2000, 20261015)


def score(u, which):
    """{stratum: (point totals, bootstrap totals, block key)} for the strata of the table."""
    w = u.rtg["bed"].parent / which
    assert C._rtg_complete(w), w
    rec = {p: np.array([r[1] for r in C.read_records(w / f"{p}.vcf.gz") if r[0] == "chr1"], dtype=np.int64) for p in PARTS}
    blocks = bs.get("bed", u.bed, G.WIDTH)
    out = {}
    for name, iv in ivs.items():
        pos = {p: rec[p][G.inside(rec[p], iv)].tolist() for p in PARTS}
        m = E.block_counts(pos, blocks)
        out[name] = (m.sum(axis=0), blocks.mult @ m, blocks.key)
    return out


def paired(a, b):
    """point dF1 and its 95% interval, HPRC (a) - GRCh38-only (b), same resampled blocks."""
    assert a[2] == b[2]
    d = E.metrics(a[1])[2] - E.metrics(b[1])[2]
    lo, hi, _, _ = E.ci(d)
    return float(E.metrics(a[0])[2] - E.metrics(b[0])[2]), lo, hi


res = {w: {uid: score(u, w) for uid, u in units.items()} for w in ("pon_calls", "raw_calls")}

# ---- check: pon_calls reproduces the published cells
bad = []
e16ci = {(r["uid"], r["stratum"]): r for r in tsv(EXP / "E16_bootstrap/strata_ci.tsv")}
e16p = {(r["set"], r["arm"], r["stratum"]): r for r in tsv(EXP / "E16_bootstrap/strata_paired.tsv")}
gbp = {(r["set"], r["kind"], r["arm"], r["stratum"]): r for r in tsv(EXP / "E10b_giab_strata/strata_paired.tsv")}
gbc = {(r["uid"], r["stratum"]): r for r in tsv(EXP / "E10b_giab_strata/strata_ci.tsv")}
R = res["pon_calls"]
for s, _ in SETS:
    lin = pick[(s, "SNV", "linear")]
    for uid in (lin, pick[(s, "SNV", "d9")], pick[(s, "SNV", "d46")]):
        tot = R[uid]["extreme"][0]; ref = e16ci[(uid, "extreme")]
        if [int(x) for x in tot] != [int(ref[k]) for k in ("tp_baseline", "fn", "tp_call", "fp")]:
            bad.append(("E16 counts", uid, list(tot)))
    for arm in ("d9", "d46"):
        a = pick[(s, "SNV", arm)]
        got = paired(R[a]["extreme"], R[lin]["extreme"]); ref = e16p[(s, arm, "extreme")]
        if any(abs(g - float(ref[k])) > 1e-6 for g, k in zip(got, ("dF1", "dF1_lo", "dF1_hi"))):
            bad.append(("E16 paired", s, arm, got))
        for st, _ in GIAB_ST:
            got = paired(R[a][st], R[lin][st]); ref = gbp[(s, "SNV", arm, st)]
            if any(abs(g - float(ref[k])) > 1e-6 for g, k in zip(got, ("dF1", "dF1_lo", "dF1_hi"))):
                bad.append(("E10b paired", s, arm, st, got))
    for arm in ("linear", "d9", "d46"):
        uid = pick[(s, "INDEL", arm)]
        f = float(E.metrics(R[uid][HP7][0])[2]); ref = float(gbc[(uid, HP7)]["F1"])
        if abs(f - ref) > 1e-6:
            bad.append(("E10b INDEL HP7 F1", uid, f, ref))
    for arm in ("d9", "d46"):  # INDEL paired cells of E10b (GIAB strata incl. homopolymers >= 7 bp)
        a, lin_i = pick[(s, "INDEL", arm)], pick[(s, "INDEL", "linear")]
        for st in [k for k, _ in GIAB_ST] + [HP7]:
            got = paired(R[a][st], R[lin_i][st]); ref = gbp[(s, "INDEL", arm, st)]
            if any(abs(g - float(ref[k])) > 1e-6 for g, k in zip(got, ("dF1", "dF1_lo", "dF1_hi"))):
                bad.append(("E10b INDEL paired", s, arm, st, got))
assert not bad, bad

# ---- table before population filtering
f3 = lambda x: "-" if np.isnan(x) else f"{x:.3f}"  # noqa: E731
ci = lambda t: f"{t[0]:+.3f} [{t[1]:+.3f}, {t[2]:+.3f}]"  # noqa: E731
X = res["raw_calls"]
cells = {}
L = ["# Supplementary Table R1.3a before population filtering (generated by region_nopon.py)", "",
     "chr1, benchmark regions, PASS at the validation recall-0.9 threshold, no population filter (rtg raw_calls; "
     "complete truth). Same call sets, blocks (1 Mb, 2,000 resamples, seed 20261015) and strata as the published "
     "table; the same code on pon_calls reproduces every published cell (checked).", "",
     "| chr1 | graph | SMaHT extreme: PASS calls (precision), GRCh38-only → HPRC | SMaHT extreme ΔF1 | "
     + " | ".join(f"GIAB {n} ΔF1" for _, n in GIAB_ST) + " | INDEL F1, homopolymers ≥ 7 bp: GRCh38-only → HPRC |",
     "|---" * (4 + len(GIAB_ST) + 1) + "|"]
for s, lab in SETS:
    lin = pick[(s, "SNV", "linear")]
    t0 = X[lin]["extreme"][0]
    for arm in ("d9", "d46"):
        a = pick[(s, "SNV", arm)]
        t1 = X[a]["extreme"][0]
        p0, p1 = E.metrics(t0)[0], E.metrics(t1)[0]
        row = [f"{int(t0[2] + t0[3]):,} ({f3(p0)}) → {int(t1[2] + t1[3]):,} ({f3(p1)})", ci(paired(X[a]["extreme"], X[lin]["extreme"]))]
        row += [ci(paired(X[a][st], X[lin][st])) for st, _ in GIAB_ST]
        i0 = float(E.metrics(X[pick[(s, "INDEL", "linear")]][HP7][0])[2]); i1 = float(E.metrics(X[pick[(s, "INDEL", arm)]][HP7][0])[2])
        row.append(f"{i0:.3f} → {i1:.3f}")
        cells[f"{s}|{arm}"] = row
        L.append(f"| {lab} | {GNAME[arm]} | " + " | ".join(row) + " |")

# per call set, SMaHT extreme and GIAB all difficult, before and after, for reading the table
L += ["", "## Per call set: SNV P / R / F1 before → after population filtering", "",
      "| set | graph | SMaHT extreme | GIAB all difficult | GIAB not difficult |", "|---|---|---|---|---|"]
for s, lab in SETS:
    for arm in ("linear", "d9", "d46"):
        uid = pick[(s, "SNV", arm)]
        c = []
        for st in ("extreme", "alldifficultregions", "notinalldifficultregions"):
            b, a = E.metrics(X[uid][st][0]), E.metrics(R[uid][st][0])
            c.append(" / ".join(f3(float(v)) for v in b) + " → " + " / ".join(f3(float(v)) for v in a))
        L.append(f"| {lab} | {'GRCh38-only' if arm == 'linear' else GNAME[arm]} | " + " | ".join(c) + " |")
# ---- INDEL version (the 7fc7396 INDEL models of the three-reference comparison, those of the published HP >= 7 column)
cells_indel = {}
L += ["", "## INDEL before population filtering", "",
      "The 7fc7396 INDEL models (seed 0; the published table's homopolymer column), same thresholds, strata and blocks. "
      "SMaHT extreme by the VCF position (anchor base), as the GIAB strata; the same code on pon_calls reproduces the "
      "published E10b INDEL cells (checked).", "",
      "chr1 truth INDELs in the benchmark BED per stratum (TP-baseline + FN; SMaHT extreme / GIAB all difficult / "
      "GIAB not difficult / GIAB homopolymers ≥ 7 bp): " + "; ".join(
          f"{s.split('_')[0]} " + " / ".join(str(int(X[pick[(s, 'INDEL', 'linear')]][st][0][:2].sum()))
                                             for st in ("extreme", "alldifficultregions", "notinalldifficultregions", HP7))
          for s in ("HG008T_Illumina", "COLO829T_Illumina")) + ".", "",
      "| chr1 INDEL | graph | SMaHT extreme: PASS calls (precision), GRCh38-only → HPRC | SMaHT extreme ΔF1 | "
      + " | ".join(f"GIAB {n} ΔF1" for _, n in GIAB_ST) + " | GIAB homopolymers ≥ 7 bp ΔF1 |",
      "|---" * (4 + len(GIAB_ST) + 1) + "|"]
for s, lab in SETS:
    lin = pick[(s, "INDEL", "linear")]
    t0 = X[lin]["extreme"][0]
    for arm in ("d9", "d46"):
        a = pick[(s, "INDEL", arm)]
        t1 = X[a]["extreme"][0]
        p0, p1 = E.metrics(t0)[0], E.metrics(t1)[0]
        row = [f"{int(t0[2] + t0[3]):,} ({f3(p0)}) → {int(t1[2] + t1[3]):,} ({f3(p1)})", ci(paired(X[a]["extreme"], X[lin]["extreme"]))]
        row += [ci(paired(X[a][st], X[lin][st])) for st, _ in GIAB_ST]
        row.append(ci(paired(X[a][HP7], X[lin][HP7])))
        cells_indel[f"{s}|{arm}"] = row
        L.append(f"| {lab} | {GNAME[arm]} | " + " | ".join(row) + " |")
L += ["", "## Per call set: INDEL P / R / F1 before → after population filtering", "",
      "| set | graph | SMaHT extreme | GIAB all difficult | GIAB not difficult | GIAB homopolymers ≥ 7 bp |", "|---|---|---|---|---|---|"]
for s, lab in SETS:
    for arm in ("linear", "d9", "d46"):
        uid = pick[(s, "INDEL", arm)]
        c = []
        for st in ("extreme", "alldifficultregions", "notinalldifficultregions", HP7):
            b, a = E.metrics(X[uid][st][0]), E.metrics(R[uid][st][0])
            c.append(" / ".join(f3(float(v)) for v in b) + " → " + " / ".join(f3(float(v)) for v in a))
        L.append(f"| {lab} | {'GRCh38-only' if arm == 'linear' else GNAME[arm]} | " + " | ".join(c) + " |")
(HERE / "tables.md").write_text("\n".join(L) + "\n")
(HERE / "cells.json").write_text(json.dumps({"uids": {"|".join(k): v for k, v in pick.items()}, "cells": cells,
                                             "cells_indel": cells_indel}, indent=1) + "\n")
print("\n".join(L))

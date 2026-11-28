#!/usr/bin/env python3
"""Checks of fp_decomp_nositef.py's output (writes checks.txt):
1. columns 1-5 of both before-filter tables equal the ver4 tables' (ver4_20261123/fpdecomp, the letter ver4.5's
   Table 3.9 and columns 1-5 of Table 2.6 = 3.8): they do not depend on the HPRC v2.1 model;
2. the SNV table equals the letter ver4.5's Table 2.6 = 3.8 in all six columns (its last column was typed in on
   2026-11-24 from the RIS report.json files);
3. the last column equals the raw FP of the v2.1 no-site-feature report.json files (INDEL: also the FP that
   results_indel.md's before-filter P / R imply)."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parents[1] / "ver4_20261123" / "fpdecomp"
LETTER_SNV = """HG008T Illumina | 28,178 | 81%; 0.2% | 16%; 87% | 3% | 1,246
HG008T PacBio HiFi | 84,853 | 85%; 0.3% | 13%; 87% | 2% | 1,491
HG008T ONT | 50,403 | 83%; 0.5% | 16%; 88% | 1% | 785
COLO829T Illumina | 29,530 | 79%; 0.1% | 11%; 89% | 10% | 1,281
COLO829T Fiber-seq | 38,677 | 71%; 0.6% | 12%; 83% | 17% | 1,040
COLO829T ONT | 30,128 | 80%; 0.7% | 12%; 90% | 8% | 1,696
HCC1395T Illumina | 77,905 | 86%; 0.1% | 14%; 87% | 0.7% | 1,783
HCC1395T PacBio HiFi | 86,347 | 86%; 0.2% | 13%; 86% | 1% | 1,778
HCC1395T ONT | 115,580 | 87%; 0.7% | 12%; 86% | 0.9% | 2,906"""


def rows(path):
    out = []
    for ln in Path(path).read_text().splitlines()[2:]:
        out.append([c.strip() for c in ln.strip().strip("|").split("|")])
    return out


lines, bad = [], 0
for kind in ("SNV", "INDEL"):
    new, old = rows(HERE / f"table_{kind}_before_filter.md"), rows(OLD / f"table_{kind}_before_filter.md")
    assert len(new) == len(old) == 9
    diff = [(n[0], i, o[i], n[i]) for n, o in zip(new, old) for i in range(5) if n[i] != o[i]]
    bad += len(diff)
    lines.append(f"{'OK  ' if not diff else 'FAIL'} {kind}: columns 1-5 (dataset, GRCh38-graph FP, two germline columns, non-germline) "
                 f"equal to the ver4 table in 9 of 9 rows{'' if not diff else ': ' + str(diff)}")
    last = [(n[0], o[5], n[5]) for n, o in zip(new, old)]
    lines.append(f"     {kind}: last column, with site features -> without: " + "; ".join(f"{d} {a} -> {b}" for d, a, b in last))
new = rows(HERE / "table_SNV_before_filter.md")
let = [[c.strip() for c in ln.split("|")] for ln in LETTER_SNV.splitlines()]
diff = [(n[0], i, l[i], n[i]) for n, l in zip(new, let) for i in range(6) if n[i] != l[i]]
bad += len(diff)
lines.append(f"{'OK  ' if not diff else 'FAIL'} SNV: all six columns equal to the letter ver4.5's Response Table 2.6 = 3.8 in 9 of 9 rows"
             f"{'' if not diff else ': ' + str(diff)}")
J = json.loads((HERE / "fp_decomp.json").read_text())
n = 0
for e in J:
    rep = json.loads((Path(e["d46_dir"]) / "eval_bed/report.json").read_text())
    for st in ("raw", "pon"):
        assert e["stages"][st]["d46_fp_records"] == rep[st]["at_threshold"]["fp"], e["dataset"]
        assert e["stages"][st]["d46_fp_keys"] == e["stages"][st]["d46_fp_records"], (e["dataset"], e["kind"], st)
        assert e["stages"][st]["lin_fp_keys"] == e["stages"][st]["lin_fp_records"], (e["dataset"], e["kind"], st)
        n += 1
    assert "nositef" in e["d46_dir"] and "linear_nositef" in e["lin_dir"]
lines.append(f"OK   HPRC v2.1 FP records == report.json fp of the no-site-feature call set, FP allele keys == FP records in both arms: {n} of {n} cell x stage")
# INDEL after the filter against results_indel.md (TP / FP columns)
res = (HERE.parents[2] / "m6_nositef_graphs_20261118/results_indel.md").read_text().splitlines()
want = {}
for ln in res:
    c = [x.strip() for x in ln.strip().strip("|").split("|")]
    if len(c) == 10 and c[2] == "HPRC v2.1":
        want[f"{c[0]} {c[1]}"] = (int(c[5]), int(c[6]))
m = 0
for e in J:
    if e["kind"] == "INDEL":
        r = e["report"]["d46_pon"]
        assert want[e["dataset"]] == (r["tp"], r["fp"]), (e["dataset"], want[e["dataset"]], r)
        m += 1
lines.append(f"OK   INDEL HPRC v2.1 TP / FP after the filter == analysis/m6_nositef_graphs_20261118/results_indel.md section 1: {m} of 9 rows")
lines.append("ALL CHECKS PASSED" if not bad else f"{bad} DIFFERENCES")
(HERE / "checks.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))

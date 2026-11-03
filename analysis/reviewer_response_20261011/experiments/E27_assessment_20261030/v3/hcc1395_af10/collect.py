#!/usr/bin/env python3
"""HCC1395T chr1 SNVs (SEQC2 high-confidence regions, after the common population filter): each caller's PASS calls
without and with an allele-fraction floor of 0.1 on its own calls (Pansoma: graph_vcf --min-af 0.1 on the A1 read
fraction; ClairS-TO FORMAT AF, DeepSomatic FORMAT VAF -> LowAF). Full truth in both. Writes tables.md, summary.json."""
import json
from pathlib import Path
HERE = Path(__file__).resolve().parent
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators")
def pans(run, s, suf):
    f = R / run / f"vcf_chr1_{s}_normfix{suf}/eval_bed/report.json"
    return json.load(open(f))["pon"]["at_threshold"] if f.exists() else None
def comp(c, p, d):
    f = E8 / d / c / f"HCC1395T_{p}/SNV/report.json"
    return json.load(open(f))["eval_bed"]["pass_pon"]["at_threshold"] if f.exists() else None
res = {}
for p in ("Illumina", "PacBio", "ONT"):
    s = f"HCC1395T_{p}_d46"
    for name, a, b in (("Pansoma six-tumor", pans(f"multicancer6_{p}_SNV", s, ""), pans(f"multicancer6_{p}_SNV", s, "_minaf0.1")),
                       ("Pansoma four-line", pans(f"multicancer4_{p}_SNV", s, ""), pans(f"multicancer4_{p}_SNV", s, "_minaf0.1")),
                       ("DeepSomatic-TO", comp("deepsomatic", p, "score"), comp("deepsomatic", p, "score_minaf0.1")),
                       ("ClairS-TO", comp("clairs_to", p, "score"), comp("clairs_to", p, "score_minaf0.1"))):
        res[f"{p}|{name}"] = dict(native=a, minaf=b)
(HERE / "summary.json").write_text(json.dumps(res, indent=1))
f = lambda x: "-" if not x else f"{x['precision']:.3f} / {x['recall']:.3f} / {x['f1']:.3f} ({x['tp_baseline']} TP, {x['fp']} FP)"
L = ["| HCC1395T chr1 SNV | caller | PASS P / R / F1 | PASS and AF >= 0.1: P / R / F1 | ΔF1 |", "|---|---|---|---|---|"]
for k, v in res.items():
    p, n = k.split("|")
    d = "-" if not (v["native"] and v["minaf"]) else f"{v['minaf']['f1'] - v['native']['f1']:+.3f}"
    L.append(f"| {p} | {n} | {f(v['native'])} | {f(v['minaf'])} | {d} |")
(HERE / "tables.md").write_text("\n".join(L) + "\n")
print("\n".join(L))

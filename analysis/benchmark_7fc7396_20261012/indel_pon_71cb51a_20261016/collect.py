"""Every RIS INDEL call set under the 7fc7396 INDEL PoN rule (homopolymer >= 7 bp: AF >= 0.01, else 1e-3) and the 71cb51a
rule (>= 8 bp: 0.01, else 1e-4), side by side. Standard library only (login node):

    PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 collect.py

For each <run>/<SRC>_pon71cb51a (rescore.sh) it reads <run>/<OLD>/eval_{bed,nobed}/report.json and the new one, where
OLD is the 7fc7396-rule evaluation of the same calls: <SRC> itself, except for the 2026-10-06 runs (*_INDEL_small),
whose <SRC> = vcf_chr1 has the rule before 7fc7396 and whose 7fc7396-rule scores are vcf_chr1_pon7fc7396. Runs whose
vcf.sh ran after the 71cb51a pull (2026-10-16 18:02) have the new rule in <SRC> already: then OLD = NEW and the row
says so. Writes tables.md and summary.json next to this script.
"""
import json
import os
from pathlib import Path

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
HERE = Path(__file__).resolve().parent
NEW = "_pon71cb51a"
PULL = 1792191762  # 2026-10-16 18:02:42 -0500, the reflog time of the 71cb51a pull into the RIS checkout


def group(run, src):
    s = "HG008T" if run.startswith("HG008") else "COLO829T"
    if "_joint_readnorm" in run:
        g = "purity joint, BQ per read"
    elif "_joint_" in run:
        g = "purity joint"
    elif "_nopc_d46" in run:
        g = "no ch6, d46"
    elif "_nopc_" in run:
        g = "no ch6, d9"
    elif "_normds" in run:
        g = "--norm dataset, d46"
    elif "_d46_" in run:
        g = "d46"
    elif "_linear_" in run:
        g = "linear"
    elif run.endswith("_small"):
        g = "d9, 45d56ad lr 4e-4"
    elif run.endswith("_h90"):
        g = "d9, 2026-10-03 h90"
    elif run.endswith("_repro"):
        g = "d9, repro"
    else:
        g = "d9"
    return s, g


def scores(d):
    out = {}
    for w in ("eval_bed", "eval_nobed"):
        f = d / w / "report.json"
        if not f.exists():
            return None
        r = json.load(open(f))
        if "pon" not in r:  # evaluated before INDEL had a PoN (the 2026-10-03 h90 runs)
            return None
        p = r["pon"]
        out[w] = {"pass_f1": p["at_threshold"]["f1"], "pass_p": p["at_threshold"]["precision"],
                  "pass_r": p["at_threshold"]["recall"], "best_f1": p["best"]["f1"], "best_p": p["best"]["precision"],
                  "best_r": p["best"]["recall"], "best_tp": p["best"]["tp_baseline"], "ceiling": p["ceiling"],
                  "truth_tagged": r["truth_pon_tags"]["tagged"], "calls_tagged": r["pon_tags"]["tagged"]}
    out["mtime"] = (d / "eval_bed" / "report.json").stat().st_mtime
    return out


def main():
    rows = []
    for nd in sorted(RUNS.glob(f"*/vcf_chr1*{NEW}")):
        run, src = nd.parent.name, nd.name[:-len(NEW)]
        old_dir = nd.parent / ("vcf_chr1_pon7fc7396" if run.endswith("_small") else src)
        new, old = scores(nd), scores(old_dir)
        if new is None:
            continue
        s, g = group(run, src)
        rows.append({"sample": s, "group": g, "run": run, "src": src, "old_dir": old_dir.name,
                     "old_is_new_rule": bool(old) and old["mtime"] > PULL, "old": old, "new": new})
    (HERE / "summary.json").write_text(json.dumps(rows, indent=1) + "\n")
    f3 = lambda x: "no PoN" if x is None else f"{x:.3f}"  # noqa: E731
    L = ["# INDEL call sets under the 7fc7396 and the 71cb51a INDEL PoN rule\n",
         "(`collect.py`; chr1 after the PoN. 7fc7396 rule: gnomAD / CoLoRSdb exact allele at AF >= 0.01 when the ALT "
         "changes a homopolymer of >= 7 bp, else AF >= 1e-3; 71cb51a rule: >= 8 bp, else AF >= 1e-4. Same calls and "
         "PASS threshold; best = the chr1 curve's best point. * = the run's own vcf.sh already had the 71cb51a rule.)\n"]
    for where, title in (("eval_bed", "in the benchmark BED"), ("eval_nobed", "whole chr1, no BED")):
        L.append(f"\n## {title}\n")
        L.append("| sample | group | run (set) | PASS F1 7fc7396 -> 71cb51a | best F1 7fc7396 -> 71cb51a | best P / R (71cb51a) | "
                 "truth tagged 7fc7396 -> 71cb51a |")
        L.append("|---|---|---|---|---|---|---|")
        for r in sorted(rows, key=lambda r: (r["sample"] != "HG008T", r["group"], r["run"], r["src"])):
            n = r["new"][where]
            o = r["old"][where] if r["old"] else {"pass_f1": None, "best_f1": None, "truth_tagged": "no PoN"}
            name = r["run"] + ("" if r["src"] == "vcf_chr1" else f" ({r['src'][len('vcf_chr1_'):]})")
            star = " *" if r["old_is_new_rule"] else ""
            L.append(f"| {r['sample']} | {r['group']} | {name}{star} | {f3(o['pass_f1'])} -> {f3(n['pass_f1'])} | "
                     f"{f3(o['best_f1'])} -> {f3(n['best_f1'])} | {f3(n['best_p'])} / {f3(n['best_r'])} | "
                     f"{o['truth_tagged']} -> {n['truth_tagged']} |")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    print(f"{len(rows)} call sets")


if __name__ == "__main__":
    main()

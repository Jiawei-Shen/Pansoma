"""ClairS-TO and DeepSomatic tumor-only against Pansoma's purity models on COLO829T P100 / P75 / P50 / P25 (chr1, benchmark BED, after our PoN).

    PANSOMA_RUNS=... python3 collect.py   -> tables.md, summary.json (next to this script)

ClairS-TO / DeepSomatic: $PANSOMA_RUNS/comparators/{clairs_to,deepsomatic}/<SET>/<KIND>/report.json (score.sh): the
caller's PASS records minus our PoN (pass_pon) and the curve of all its records minus our PoN by QUAL (all_pon, best F1 =
chr1 upper bound). DeepSomatic P75-P25: run_deepsomatic.sh (submit_deepsomatic.sh); P100: the user's runs.
Pansoma INDEL: the ppr8 joint models COLO829_<P>_INDEL_joint_ppr8 (PASS, post-PoN best-F1 threshold, chr1 best).
Pansoma SNV: the 7fc7396 joint models COLO829_<P>_SNV_joint_7fc7396_gpu4 (the SNV recipe is unchanged; PASS, chr1 best).
E8 (wangcluster, rtg 3.12.1) scored both callers' P100 VCFs; the P100 rows here repeat that scoring on RIS as a check.
"""
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
CMP = RUNS / "comparators"
CALLERS = {"clairs_to": "ClairS-TO", "deepsomatic": "DeepSomatic"}
E8 = HERE.parents[1] / "reviewer_response_20261011/experiments/E8_comparators/summary.json"
PLATFORMS = ("Illumina", "fiberseq", "ONT")
PURITIES = ("P100", "P75", "P50", "P25")
PANSOMA = {"INDEL": "COLO829_{p}_INDEL_joint_ppr8", "SNV": "COLO829_{p}_SNV_joint_7fc7396_gpu4"}
# the P100 single-set models (d9, seeds 0-2) on the mixtures (../p100_on_dilutions/submit.sh); COLO829T Illumina SNV at
# the epochs the seed table reports (6 / 5 / 6)
P100_MODELS = {"INDEL": lambda p: [f"COLO829_{p}_INDEL_d9_ppr8{x}" for x in ("", "_s1", "_s2")],
               "SNV": lambda p: ["COLO829_Illumina_SNV_d9_normbq_sweep_e06", "COLO829_Illumina_SNV_d9_normbq_s1_sweep_e05",
                                 "COLO829_Illumina_SNV_d9_normbq_s2_sweep_e06"] if p == "Illumina"
               else [f"COLO829_{p}_SNV_d9_normbq{x}" for x in ("", "_s1", "_s2")]}


def set_of(p, f):
    return f"COLO829T_{p}" if f == "P100" else f"COLO829T_{p}_{f}"


def prf(a):
    return dict(p=a["precision"], r=a["recall"], f1=a["f1"], tp=a["tp_baseline"], fp=a["fp"])


def comparator(caller, s, kind):
    path = CMP / caller / s / kind / "report.json"
    if not path.exists():
        return None
    r = json.loads(path.read_text())["eval_bed"]
    return dict(**prf(r["pass_pon"]["at_threshold"]), best=(r["all_pon"]["best"] or {}).get("f1"),
                ceiling=r["all_pon"]["ceiling"])


def pansoma(run, s, suffix=""):
    path = RUNS / run / (f"vcf_chr1_{s}{suffix}" if s else f"vcf_chr1{suffix}") / "eval_bed" / "report.json"
    if not path.exists():
        return None
    r = json.loads(path.read_text())["pon"]
    return dict(**prf(r["at_threshold"]), best=(r["best"] or {}).get("f1"), ceiling=r["ceiling"])


def p100_model(kind, p, f, suffix=""):
    """The P100 models' scores on purity f (their own chr1 test at P100): mean P / R / F1 of the seeds, F1 per seed."""
    runs = P100_MODELS[kind](p)
    scores = [pansoma(r, "", suffix) if f == "P100" else pansoma(r, set_of(p, f), suffix) for r in runs]
    scores = [x for x in scores if x]
    if not scores:
        return None
    mean = {k: sum(x[k] for x in scores) / len(scores) for k in ("p", "r", "f1")}
    return dict(**mean, seeds=[x["f1"] for x in scores], runs=runs)


def cell(x):
    return "-" if x is None else f"{x['p']:.3f} / {x['r']:.3f} / **{x['f1']:.3f}**"


def f3(v):
    return "-" if v is None else f"{v:.3f}"


def main():
    out, L = {}, ["# ClairS-TO and DeepSomatic vs Pansoma purity models, COLO829T chr1\n",
                  "(`collect.py`; benchmark BED, after our PoN; cells P / R / **F1**. ClairS-TO / DeepSomatic: the PASS "
                  "records minus our PoN; best = the chr1 curve of all records after the PoN, by QUAL (an upper bound); "
                  "'-' = not scored. Pansoma: one joint "
                  "model per platform on P100-P25, PASS = its call rule; INDEL also at the post-PoN best-F1 threshold. P100 model: the "
                  "P100 single-set d9 models, seeds 0-2, mean P / R / F1 and each seed's F1; at P100 their own chr1 test.)\n"]
    for kind in ("INDEL", "SNV"):
        L += [f"\n## {kind}\n"]
        head = "| platform | purity | ClairS-TO PASS | ClairS-TO best | DeepSomatic PASS | DeepSomatic best | Pansoma PASS |"
        head += " Pansoma post-PoN best | Pansoma best |" if kind == "INDEL" else " Pansoma best |"
        head += " P100 model PASS (seeds) |" + (" P100 model post-PoN best |" if kind == "INDEL" else "")
        L += [head, "|---" * (head.count("|") - 1) + "|"]
        for p in PLATFORMS:
            run = PANSOMA[kind].format(p=p)
            for f in PURITIES:
                s = set_of(p, f)
                c, d, m = comparator("clairs_to", s, kind), comparator("deepsomatic", s, kind), pansoma(run, s)
                q = pansoma(run, s, "_postpon") if kind == "INDEL" else None
                a = p100_model(kind, p, f)
                aq = p100_model(kind, p, f, "_postpon") if kind == "INDEL" else None
                out[f"{kind}/{s}"] = dict(clairs_to=c, deepsomatic=d, pansoma=m, pansoma_postpon=q, pansoma_run=run,
                                          p100_model=a, p100_model_postpon=aq)
                row = f"| {p} | {f} | {cell(c)} | {f3(c and c['best'])} | {cell(d)} | {f3(d and d['best'])} | {cell(m)} |"
                row += f" {f3(q and q['f1'])} | {f3(m and m['best'])} |" if kind == "INDEL" else f" {f3(m and m['best'])} |"
                row += (f" {a['p']:.3f} / {a['r']:.3f} / **{a['f1']:.3f}** ({' / '.join(f'{v:.3f}' for v in a['seeds'])}) |"
                        if a else " - |")
                row += f" {f3(aq and aq['f1'])} |" if kind == "INDEL" else ""
                L.append(row)
    if E8.exists():  # the check: E8's P100 scores (wangcluster) against the RIS re-scoring
        e8 = json.loads(E8.read_text())
        L += ["\n## Check: P100, E8 (wangcluster) vs RIS, PASS F1 after our PoN\n",
              "| caller | set | kind | E8 | RIS |", "|---|---|---|---|---|"]
        rows = {(r["caller"], r["dataset"], r["kind"]): r for r in e8["datasets"]}
        for caller, label in CALLERS.items():
            for p in PLATFORMS:
                for kind in ("SNV", "INDEL"):
                    r = rows.get((caller, f"COLO829T_{p}", kind))
                    e8f = r and r["eval_bed"]["pass_pon"]["at"]["f1"]
                    c = out.get(f"{kind}/COLO829T_{p}", {}).get(caller)
                    L.append(f"| {label} | COLO829T_{p} | {kind} | {f3(e8f)} | {f3(c and c['f1'])} |")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    (HERE / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

"""Without any population filtering vs with it, chr1 in the benchmark BED, PASS: Pansoma (raw PASS = before our PoN; the
model has no PoN of its own) next to ClairS-TO / DeepSomatic-TO without their own panel-of-normals step, their default,
and default + our PoN (E8 score/). Without the PoN step:
  ClairS-TO: rerun with --disable_nonsomatic_tagging (clairs_to_rerun/; E8_comparators/score_nopon_rerun/). A FILTER
             rewrite is not equivalent: the tagging runs before realignment, postfilter and haplotype filtering, which
             only process PASS records; the rewrite (score_nopon/) is listed separately for comparison.
  DeepSomatic-TO: FILTER PON -> PASS (rewrite_nopon.py; score_nopon/), exact: postprocess_variants applies the PoN last
             and only to PASS calls.
  Pansoma SNV: E27 ablation/main_runs.json (validation-selected checkpoints, three seeds, mean): GRCh38-only = setting train
               linear, d9 / d46 = setting normbq (the E8 latest models); HCC1395T: multicancer4 d46 normfix (E12b, one run).
  Pansoma INDEL: final release models (ppr8, d46), three seeds, mean.
usage: /wanglab/jshen/anaconda3/bin/python collect_with_pansoma.py   (writes tables_with_pansoma.md / summary_with_pansoma.json here)"""
import json
import statistics as st
from pathlib import Path

HERE = Path(__file__).resolve().parent
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
MAIN = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E27_assessment_20261030/ablation/main_runs.json")
DATASETS = [("HG008T", "HG008", p) for p in ("Illumina", "PacBio", "ONT")] + \
           [("COLO829T", "COLO829", p) for p in ("Illumina", "fiberseq", "ONT")] + \
           [("HCC1395T", None, p) for p in ("Illumina", "PacBio", "ONT")]
CALLERS = (("clairs_to", "ClairS-TO"), ("deepsomatic", "DeepSomatic-TO"))


def at(block):
    a = block["at_threshold"]
    return dict(p=a["precision"], r=a["recall"], f1=a["f1"])


def mean(xs):
    return {k: st.mean(x[k] for x in xs) for k in ("p", "r", "f1")}


def fmt(x):
    return f"{x['p']:.3f} / {x['r']:.3f} / **{x['f1']:.3f}**" if x else "-"


def comparators(c, ds, kind):
    n = json.loads((E8 / ("score_nopon_rerun" if c == "clairs_to" else "score_nopon") / c / ds / kind / "report.json").read_text())["eval_bed"]
    o = json.loads((E8 / "score" / c / ds / kind / "report.json").read_text())["eval_bed"]
    return dict(none=at(n["pass"]), own=at(o["pass"]), own_ours=at(o["pass_pon"]))


def pansoma_snv(main, sample, short, p):
    out = {}
    if short is None:  # HCC1395T
        r = json.loads((RUNS / f"multicancer4_{p}_SNV" / f"vcf_chr1_HCC1395T_{p}_d46_normfix" / "eval_bed" / "report.json").read_text())
        out["d46"] = dict(none=at(r["raw"]), ours=at(r["pon"]))
        return out
    for arm, setting in (("linear", "train"), ("d9", "normbq"), ("d46", "normbq")):
        rows = main[setting][f"{short}|{p}|{arm}"]
        out[arm] = dict(none=mean([dict(p=x["rawP"], r=x["rawR"], f1=x["rawF1"]) for x in rows]),
                        ours=mean([dict(p=x["P"], r=x["R"], f1=x["F1"]) for x in rows]))
    return out


def pansoma_indel(sample, p):
    plat = "PacBio" if p == "fiberseq" else p
    reps = [json.loads(f.read_text()) for s in ("s0", "s1", "s2")
            if (f := RUNS / f"release_{s}_{plat}_INDEL" / f"vcf_ppr8_chr1_{sample}_{p}_d46" / "eval_bed" / "report.json").exists()]
    return dict(d46=dict(none=mean([at(r["raw"]) for r in reps]), ours=mean([at(r["pon"]) for r in reps]), seeds=len(reps)))


def main():
    m = json.loads(MAIN.read_text())
    summary, L = {}, ["# Without any population filtering vs with it (chr1, benchmark BED, PASS): P / R / F1", "",
                      "Pansoma: raw PASS (before our PoN; no PoN of its own) and PASS after our PoN. ClairS-TO / DeepSomatic-TO: PASS",
                      "without their own panel-of-normals step (ClairS-TO rerun with --disable_nonsomatic_tagging; DeepSomatic PON ->",
                      "PASS, exact, its GERMLINE class and population-AF input channel stay), their default PASS, and default + our PoN.", ""]
    for kind in ("SNV", "INDEL"):
        L += [f"## {kind}", "", "| dataset | caller | no population filtering | own filter (default) | + our PoN |", "|---|---|---|---|---|"]
        for sample, short, p in DATASETS:
            ds = f"{sample}_{p}"
            pz = pansoma_snv(m, sample, short, p) if kind == "SNV" else pansoma_indel(sample, p)
            for arm, v in pz.items():
                name = {"linear": "Pansoma GRCh38-only", "d9": "Pansoma HPRC v1.1", "d46": "Pansoma HPRC v2.1"}[arm]
                if sample == "HCC1395T" and kind == "SNV":
                    name += " (multicancer4)"
                L.append(f"| {ds} | {name} | {fmt(v['none'])} | - | {fmt(v['ours'])} |")
                summary[f"{kind}|{ds}|pansoma_{arm}"] = v
            for c, cname in CALLERS:
                try:
                    v = comparators(c, ds, kind)
                except FileNotFoundError:
                    L.append(f"| {ds} | {cname} | (not scored) | | |")
                    continue
                L.append(f"| {ds} | {cname} | {fmt(v['none'])} | {fmt(v['own'])} | {fmt(v['own_ours'])} |")
                summary[f"{kind}|{ds}|{c}"] = v
        L.append("")
    L += ["## ClairS-TO without its PoN step: rerun vs FILTER rewrite (PASS records, P / R / F1)", "",
          "| dataset | kind | rerun --disable_nonsomatic_tagging | rewrite NonSomatic -> PASS |", "|---|---|---|---|"]
    for kind in ("SNV", "INDEL"):
        for sample, short, p in DATASETS:
            ds = f"{sample}_{p}"
            try:
                rr = json.loads((E8 / "score_nopon_rerun" / "clairs_to" / ds / kind / "report.json").read_text())
                rw = json.loads((E8 / "score_nopon" / "clairs_to" / ds / kind / "report.json").read_text())
            except FileNotFoundError:
                continue
            L.append(f"| {ds} | {kind} | {fmt(at(rr['eval_bed']['pass']))} ({rr['counts']['pass']:,}) | "
                     f"{fmt(at(rw['eval_bed']['pass']))} ({rw['counts']['pass']:,}) |")
            summary[f"{kind}|{ds}|clairs_to_rewrite"] = at(rw["eval_bed"]["pass"])
    L.append("")
    (HERE / "tables_with_pansoma.md").write_text("\n".join(L) + "\n")
    (HERE / "summary_with_pansoma.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

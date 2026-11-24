#!/usr/bin/env python
"""Response Figure 1.3 (ver4): models trained with / without the haplotype-count (path-count) channel.

Data step. Reads the rtg vcfeval reports (eval_bed/report.json, benchmark BED) of the only with / without pairs that
exist on the current SNV pipeline: the COLO829T single-tumor SNV models on Illumina and Fiber-seq, HPRC v2.1 (d46)
and HPRC v1.1 (d9), chr1 held out.

  with the channel    R/snv_newpipe_20261109/COLO829_<P>_SNV_<g>_newpipe{,_s1,_s2}   (seeds 0, 1, 2)
  without the channel R/snv_newpipe_20261109/COLO829_<P>_SNV_nopc_<g>_newpipe        (seed 0, --drop-planes path_count)

P / R / F1 are computed here from the counts of ["pon"|"raw"]["at_threshold"] (PASS calls; "pon" = after the
population filter, "raw" = before), so that the deltas are not differences of rounded values:
  precision = tp_call / (tp_call + fp), recall = tp_baseline / baseline, F1 = their harmonic mean.

Checks (all asserted; a failure stops the script):
  - the values computed from the counts equal the report's own 4-decimal precision / recall / f1;
  - tp_baseline + fn == baseline (truth SNVs in the BED) == truth_records in every report and stage, and the truth
    count is the same in all 16 reports (the eight models of a graph share the test set; d46 and d9 share the truth and the BED);
  - truth file, BED, kind (SNV), chromosomes (chr1) and the four population resources are identical in all reports;
  - the args.json of every with-seed and of the without model differ from the seed-0 with model only in
    output / seed / drop_planes, and drop_planes is [] (with) or ["path_count"] (without);
  - every model was called by the same PASS rule (the rule text of steps.graph_vcf.call without the threshold value).

Outputs (this directory): values.json (everything, full precision), values.tsv (one line per model and stage),
table.md (the plotted values and the deltas).

Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python pathcount_data.py
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
N = R / "snv_newpipe_20261109"

PLATFORMS = (("Illumina", "Illumina"), ("fiberseq", "Fiber-seq"))
GRAPHS = (("d46", "HPRC v2.1"), ("d9", "HPRC v1.1"))
MODELS = (("with_s0", "with, seed 0"), ("with_s1", "with, seed 1"), ("with_s2", "with, seed 2"),
          ("without_s0", "without, seed 0"))
METRICS = (("precision", "pon", "P"), ("recall", "pon", "R"), ("f1", "pon", "F1"), ("f1", "raw", "F1 pre"))


def metrics_from_counts(a, baseline):
    """P / R / F1 from the counts of one at_threshold block; asserts they equal the report's rounded values."""
    tpb, tpc, fp, fn = a["tp_baseline"], a["tp_call"], a["fp"], a["fn"]
    assert tpb + fn == baseline, (tpb, fn, baseline)
    p = tpc / (tpc + fp) if tpc + fp else 0.0
    r = tpb / baseline
    f = 2 * p * r / (p + r) if p + r else 0.0
    for name, v in (("precision", p), ("recall", r), ("f1", f)):
        assert abs(v - a[name]) < 5.1e-5, (name, v, a[name])
    return {"tp_baseline": tpb, "tp_call": tpc, "fp": fp, "fn": fn, "truth": baseline,
            "precision": p, "recall": r, "f1": f,
            "report_precision": a["precision"], "report_recall": a["recall"], "report_f1": a["f1"]}


def load(report_path):
    """One eval_bed/report.json -> {"pon": {...}, "raw": {...}, "meta": {...}}."""
    d = json.load(open(report_path))
    out = {}
    assert d["pon"]["baseline"] == d["raw"]["baseline"], report_path      # truth in the BED, same before / after
    for stage in ("pon", "raw"):
        out[stage] = metrics_from_counts(d[stage]["at_threshold"], d[stage]["baseline"])
    g = d.get("steps", {}).get("graph_vcf", {})                           # older reports have no "steps" / "pon_files"
    out["meta"] = {"report": str(report_path), "calls": d["calls"], "truth": d["truth"], "bed": d["bed"],
                   "kind": d["kind"], "chromosomes": d["chromosomes"], "truth_records": d["truth_records"],
                   "truth_in_bed": d["pon"]["baseline"], "pon_files": d.get("pon_files"),
                   "threshold": d["threshold"], "call_rule": g.get("call"), "checkpoint": g.get("checkpoint")}
    return out


def run_dir(p, g, model):
    if model == "without_s0":
        return N / f"COLO829_{p}_SNV_nopc_{g}_newpipe"
    s = int(model[-1])
    return N / f"COLO829_{p}_SNV_{g}_newpipe{'' if s == 0 else f'_s{s}'}"


def vcf_dir(p, g, model):
    return run_dir(p, g, model) / f"vcf_chr1_COLO829T_{p}{'_d46' if g == 'd46' else ''}"


def collect():
    rows, checks = [], []
    for p, pname in PLATFORMS:
        for g, gname in GRAPHS:
            row = {"platform": p, "platform_name": pname, "graph": g, "graph_name": gname,
                   "label": f"COLO829T {pname}, {gname}", "models": {}}
            base_args = json.load(open(run_dir(p, g, "with_s0") / "args.json"))
            for m, _ in MODELS:
                row["models"][m] = load(vcf_dir(p, g, m) / "eval_bed" / "report.json")
                args = json.load(open(run_dir(p, g, m) / "args.json"))
                diff = sorted(k for k in set(args) | set(base_args) if args.get(k) != base_args.get(k))
                want = [] if m == "with_s0" else ["drop_planes", "output"] if m == "without_s0" else ["output", "seed"]
                assert diff == want, (p, g, m, diff)
                assert args["drop_planes"] == (["path_count"] if m == "without_s0" else []), (p, g, m)
                assert args["seed"] == (0 if m == "without_s0" else int(m[-1])), (p, g, m)
                row["models"][m]["meta"]["tensors"] = args["tensors"]
                row["models"][m]["meta"]["run"] = str(run_dir(p, g, m))
                assert Path(row["models"][m]["meta"]["checkpoint"]).parent == run_dir(p, g, m), (p, g, m)
            w0, wo = row["models"]["with_s0"], row["models"]["without_s0"]
            row["delta"], row["seed_range"], row["seed_mean"], row["without_vs_range"], row["delta_mean"] = {}, {}, {}, {}, {}
            for key, stage, short in METRICS:
                seeds = [row["models"][f"with_s{s}"][stage][key] for s in range(3)]
                lo, hi, x = min(seeds), max(seeds), wo[stage][key]
                row["delta"][short] = w0[stage][key] - x
                row["seed_range"][short] = [lo, hi]
                row["seed_mean"][short] = sum(seeds) / 3
                row["delta_mean"][short] = sum(seeds) / 3 - x
                row["without_vs_range"][short] = "below" if x < lo else "above" if x > hi else "inside"
            rows.append(row)

    # checks across the 16 reports
    metas = [m["meta"] for r in rows for m in r["models"].values()]
    truths = {m["truth_records"] for m in metas} | {m["truth_in_bed"] for m in metas}
    assert len(truths) == 1, truths                    # all chr1 truth SNVs of COLO829T lie in the benchmark BED
    for field in ("truth", "bed", "kind"):
        assert len({m[field] for m in metas}) == 1, field
    assert all(m["chromosomes"] == ["chr1"] and m["kind"] == "SNV" for m in metas)
    assert len({tuple(m["pon_files"]) for m in metas}) == 1
    rules = {re.sub(r"[0-9.]+(e-?\d+)? \(", "<thr> (", m["call_rule"], count=1) for m in metas}
    assert len(rules) == 1, rules
    for r in rows:                                     # one test set (one tensor set) per row
        assert len({json.dumps(m["meta"]["tensors"]) for m in r["models"].values()}) == 1, r["label"]
        assert len({Path(m["meta"]["calls"]).name for m in r["models"].values()}) == 1, r["label"]
    checks.append(f"truth SNVs (TP + FN = baseline = truth_records) = {truths.pop()} in all 16 reports, both stages")
    checks.append("P / R / F1 from the counts equal the 4-decimal values of every report (32 blocks)")
    checks.append("truth VCF, BED, kind SNV, chr1, four population resources identical in all 16 reports")
    checks.append("args.json: with-seeds differ from seed 0 only in output / seed; the without model only in output / drop_planes = [path_count]")
    checks.append("PASS rule identical: " + rules.pop())
    return rows, checks, metas


def f3(x):
    return f"{x:.3f}"


def d3(x):
    return f"{x:+.3f}".replace("-", "−")


def write(rows, checks, metas):
    json.dump({"rows": rows, "checks": checks, "metrics": [list(m) for m in METRICS]},
              open(HERE / "values.json", "w"), indent=1)

    with open(HERE / "values.tsv", "w") as f:
        f.write("tumor\tplatform\tgraph\tmodel\tseed\tstage\ttruth\tTP_baseline\tTP_call\tFP\tFN\tprecision\trecall\tF1\treport_json\n")
        for r in rows:
            for m, _ in MODELS:
                for stage, sname in (("pon", "after_population_filter"), ("raw", "before_population_filter")):
                    a = r["models"][m][stage]
                    f.write("\t".join(map(str, ["COLO829T", r["platform_name"], r["graph_name"],
                                                "with_channel" if m.startswith("with_") else "without_channel", m[-1], sname,
                                                a["truth"], a["tp_baseline"], a["tp_call"], a["fp"], a["fn"],
                                                f"{a['precision']:.6f}", f"{a['recall']:.6f}", f"{a['f1']:.6f}",
                                                r["models"][m]["meta"]["report"]])) + "\n")
        f.write("# deltas (with seed 0 - without) and seed range of the three with-channel seeds\n")
        f.write("#tumor\tplatform\tgraph\tmetric\tdelta_seed0_minus_without\twith_min_of_3_seeds\twith_max_of_3_seeds\twith_mean_of_3_seeds\tdelta_mean3_minus_without\twithout_vs_seed_range\n")
        for r in rows:
            for _, _, short in METRICS:
                f.write("\t".join(["#COLO829T", r["platform_name"], r["graph_name"], short.replace(" ", "_"),
                                   f"{r['delta'][short]:+.6f}", f"{r['seed_range'][short][0]:.6f}",
                                   f"{r['seed_range'][short][1]:.6f}", f"{r['seed_mean'][short]:.6f}",
                                   f"{r['delta_mean'][short]:+.6f}", r["without_vs_range"][short]]) + "\n")

    truth = metas[0]["truth_records"]
    L = []
    P = L.append
    P("# Response Figure 1.3 (ver4): with / without the haplotype-count channel, plotted values")
    P("")
    P(f"COLO829T single-tumor SNV models of the current pipeline (`snv_newpipe_20261109`), chr1 (held out), benchmark "
      f"regions, {truth:,} truth SNVs in every cell. PASS calls. P / R / F1 = after the population filter (four resource "
      "VCFs); F1 before = the same PASS calls before the filter. Values are computed from the TP / FP / FN counts of "
      "`eval_bed/report.json`; the deltas are differences of the unrounded values, so a delta can differ by 0.001 from "
      "the difference of two rounded cells.")
    P("")
    P("## Table A. Every plotted value")
    P("")
    P("| COLO829T chr1 SNV | Graph | Model | TP / FP / FN (after) | Precision | Recall | F1 | F1 before the filter |")
    P("|---|---|---|---|---|---|---|---|")
    for r in rows:
        for m, mname in MODELS:
            a, b = r["models"][m]["pon"], r["models"][m]["raw"]
            name = ("With channel, seed " + m[-1]) if m.startswith("with_") else "Without channel, seed 0"
            P(f"| {r['platform_name']} | {r['graph_name']} | {name} | {a['tp_baseline']:,} / {a['fp']:,} / {a['fn']:,} | "
              f"{f3(a['precision'])} | {f3(a['recall'])} | {f3(a['f1'])} | {f3(b['f1'])} |")
    P("")
    P("## Table B. Deltas of panel b (with, seed 0 − without) and the seed spread")
    P("")
    P("| COLO829T chr1 SNV | Graph | ΔP | ΔR | ΔF1 | ΔF1 before the filter | F1 with: range of 3 seeds (mean) | F1 without | Without vs the seed range (P / R / F1 / F1 before) | ΔF1, mean of 3 seeds − without |")
    P("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        lo, hi = r["seed_range"]["F1"]
        P(f"| {r['platform_name']} | {r['graph_name']} | " + " | ".join(d3(r["delta"][s]) for _, _, s in METRICS)
          + f" | {f3(lo)}–{f3(hi)} ({f3(r['seed_mean']['F1'])}) | {f3(r['models']['without_s0']['pon']['f1'])} | "
          + " / ".join(r["without_vs_range"][s] for _, _, s in METRICS) + f" | {d3(r['delta_mean']['F1'])} |")
    P("")
    P("## Table C. Range of the three with-channel seeds per metric")
    P("")
    P("| COLO829T chr1 SNV | Graph | Precision, min–max | Recall, min–max | F1, min–max | F1 before, min–max |")
    P("|---|---|---|---|---|---|")
    for r in rows:
        P(f"| {r['platform_name']} | {r['graph_name']} | "
          + " | ".join(f"{f3(r['seed_range'][s][0])}–{f3(r['seed_range'][s][1])}" for _, _, s in METRICS) + " |")
    P("")
    absd = [abs(r["delta"][s]) for r in rows for _, _, s in METRICS]
    P(f"Largest |Δ| in panel b: {max(absd):.3f}.")
    P("")
    P("Checks done by `pathcount_data.py` (asserted):")
    for c in checks:
        P(f"- {c}")
    P("")
    P("Runs: `" + str(N) + "/COLO829_{Illumina,fiberseq}_SNV_{d46,d9}_newpipe{,_s1,_s2}` (with) and "
      "`COLO829_{Illumina,fiberseq}_SNV_nopc_{d46,d9}_newpipe` (without); VCF dirs `vcf_chr1_COLO829T_<P>_d46` (v2.1) "
      "and `vcf_chr1_COLO829T_<P>` (v1.1).")
    (HERE / "table.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    rows, checks, metas = collect()
    write(rows, checks, metas)
    print((HERE / "table.md").read_text())
    print("tensors per row:")
    for r in rows:
        print(" ", r["label"], r["models"]["with_s0"]["meta"]["tensors"], Path(r["models"]["with_s0"]["meta"]["calls"]).name)

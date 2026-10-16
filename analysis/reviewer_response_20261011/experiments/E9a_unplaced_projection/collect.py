"""E9a/E9b/E9e: collect the per-run outputs into the experiment's tables (standard library only).

    python collect.py            (writes results/ next to this script and the per-run coordination TSVs)

Reads <scratch>/runs/<RUN>/<VCFDIR>/{projection.json,metrics.json} (project_unplaced.py, evaluate_complete.py), the
published RIS numbers of analysis/benchmark_7fc7396_20261012/*/summary.json (same run name and set), and the linear-graph
CPU runs' own eval reports. Writes:
    results/projection_by_run.tsv   unplaced PASS calls per run and their classes (GRCh38 i / ii / iii, CHM13)
    results/metrics_by_run.tsv      P / R / F1 per run, bed|nobed x raw|pon x scoring variant
    results/tables.md               the markdown tables of the README
    results/summary.json
and <scratch>/<RUN>.projected.tsv (all vcf dirs of the run), the file E9c reads.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SCRATCH = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
BENCH = REPO / "analysis" / "benchmark_7fc7396_20261012"
VARIANTS = ("placed", "placed_unplaced_fp", "projected_rtg", "projected_ii", "projected_worst")


def published():
    """(run, set or None) -> {'eval_bed': ..., 'eval_nobed': ...} from every summary.json of the RIS benchmark."""
    out = {}

    def walk(x, run=None, set_name=None):
        if isinstance(x, dict):
            run = x.get("run", run) if isinstance(x.get("run"), str) else run
            set_name = x.get("set", set_name) if isinstance(x.get("set"), str) else set_name
            if "eval_bed" in x and run:
                out.setdefault((run, set_name), dict(eval_bed=x["eval_bed"], eval_nobed=x.get("eval_nobed")))
                out.setdefault((run, None), dict(eval_bed=x["eval_bed"], eval_nobed=x.get("eval_nobed")))
            for v in x.values():
                walk(v, run, set_name)
        elif isinstance(x, list):
            for v in x:
                walk(v, run, set_name)
    for f in sorted(BENCH.glob("*/summary.json")):
        walk(json.loads(f.read_text()))
    return out


def fmt(x, d=3):
    return "-" if x is None else (f"{x:.{d}f}" if isinstance(x, float) else f"{x:,}" if isinstance(x, int) else str(x))


def pct(a, b):
    return "-" if not b else f"{100 * a / b:.1f}%"


def sample_platform(set_name):
    parts = set_name.split("_")
    return parts[0], parts[1], "_".join(parts[2:])


def main():
    pub = published()
    rows_p, rows_m, summary = [], [], []
    for mj in sorted(SCRATCH.glob("runs/*/*/projection.json")):
        pdir = mj.parent
        pj = json.loads(mj.read_text())
        m = json.loads((pdir / "metrics.json").read_text()) if (pdir / "metrics.json").exists() else None
        run, vdir = pdir.parent.name, pdir.name
        sample, platform, _ = sample_platform(pj["set"])
        cp, u = pj.get("classes_pass", {}), pj["unplaced_pass"]
        ch = pj.get("chm13_pass", {})
        det = pj.get("details_pass", {})
        cbc = pj.get("chm13_by_class_pass", {})
        rec = dict(run=run, vcf_dir=vdir, set=pj["set"], sample=sample, platform=platform, kind=pj["kind"],
                   graph=pj["graph"], unplaced_records=pj["unplaced_records"], unplaced_pass=u,
                   class_i=cp.get("i", 0), class_ii=cp.get("ii", 0), class_iii=cp.get("iii", 0),
                   i_exact=det.get("i:exact", 0), i_branch_base_differs=det.get("i:branch_base_differs", 0),
                   i_composite=det.get("i:composite", 0),
                   ii_in_insertion=det.get("ii:in_insertion", 0), ii_alt_is_grch38=det.get("ii:alt_is_grch38_base", 0),
                   ii_indel_context_gap=det.get("ii:near_branch_indel", det.get("ii:indel_context_gap", 0)), ii_too_long=det.get("ii:too_long", 0),
                   iii_no_left=det.get("iii:no_left_anchor", 0), iii_no_right=det.get("iii:no_right_anchor", 0),
                   iii_inconsistent=det.get("iii:inconsistent_anchors", 0),
                   chm13_direct=ch.get("direct", 0), chm13_anchor_interval=ch.get("anchor_interval", 0),
                   chm13_none=ch.get("none", 0),
                   not_i_chm13_direct=sum(v for k, v in cbc.items() if not k.startswith("i:") and k.endswith(":direct")),
                   merged_with_placed=pj.get("merged_with_placed", 0), seconds=round(pj.get("seconds", 0)))
        if m:
            rec.update({k: m["pass_counts"].get(k) for k in ("placed_raw", "placed_pon", "complete_raw", "complete_pon",
                                                             "frac_unplaced_raw_before", "frac_unplaced_pon_before",
                                                             "frac_unplaced_raw_after", "frac_unplaced_pon_after",
                                                             "projected_class_i_pass_pon_tagged")})
        rows_p.append(rec)
        pb = pub.get((run, pj["set"])) or pub.get((run, None))
        for ev in ("bed_raw", "bed_pon", "nobed_raw", "nobed_pon"):
            if not m:
                continue
            e = m["eval"][ev]
            b, label = ev.split("_")
            pubrow = None
            if pb and pb.get(f"eval_{b}"):
                x = pb[f"eval_{b}"][label]
                pubrow = dict(tp=x["tp"], fp=x["fp"], fn=x["fn"], precision=x["p"], recall=x["r"], f1=x["f1"])
            for var in ("published",) + VARIANTS:
                v = pubrow if var == "published" else e[var]
                if v is None:
                    continue
                rows_m.append(dict(run=run, vcf_dir=vdir, set=pj["set"], sample=sample, platform=platform,
                                   kind=pj["kind"], graph=pj["graph"], eval=ev, variant=var, tp=v["tp"], fp=v["fp"],
                                   fn=v["fn"], precision=v["precision"], recall=v["recall"], f1=v["f1"],
                                   class_ii_calls=e["class_ii_calls"], class_ii_tp=e["class_ii_tp"],
                                   class_iii_calls=e["class_iii_calls"]))
        summary.append(dict(projection=pj, metrics=m, published=pb))
    # linear-graph runs: their own reports (no projection)
    linear = []
    for rep in sorted(RUNS.glob("*/vcf_chr1*/eval_bed/report.json")):
        r = json.loads(rep.read_text())
        if "linear" not in rep.parent.parent.parent.name:
            continue
        linear.append(dict(run=rep.parent.parent.parent.name, kind=r["kind"], unplaced_pass=r.get("unplaced_pass"),
                           pon=r["pon"]["at_threshold"], raw=r["raw"]["at_threshold"]))
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    for name, rows in (("projection_by_run.tsv", rows_p), ("metrics_by_run.tsv", rows_m)):
        cols = list(rows[0]) if rows else []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
        with open(out / name, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join("" if r.get(c) is None else (f"{r[c]:.4f}" if isinstance(r.get(c), float) else str(r[c]))
                                  for c in cols) + "\n")
    (out / "summary.json").write_text(json.dumps(dict(runs=summary, linear=linear), indent=1) + "\n")
    # markdown
    md = ["## chr1: unplaced PASS calls and their projection (E9a GRCh38, E9b CHM13)", "",
          "| run | set | kind | graph | unplaced PASS | i (GRCh38 allele) | ii (interval) | iii (no anchor) | "
          "CHM13 direct | CHM13 anchor interval | CHM13 none | unplaced share of PASS, raw before -> after | "
          "after PoN before -> after |",
          "|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|"]
    for r in rows_p:
        u = r["unplaced_pass"]
        md.append(f"| `{r['run']}` | {r['set']} | {r['kind']} | {r['graph']} | {u:,} | {r['class_i']:,} ({pct(r['class_i'], u)}) | "
                  f"{r['class_ii']:,} ({pct(r['class_ii'], u)}) | {r['class_iii']:,} ({pct(r['class_iii'], u)}) | "
                  f"{r['chm13_direct']:,} ({pct(r['chm13_direct'], u)}) | {r['chm13_anchor_interval']:,} | {r['chm13_none']:,} | "
                  f"{fmt(r.get('frac_unplaced_raw_before'))} -> {fmt(r.get('frac_unplaced_raw_after'))} | "
                  f"{fmt(r.get('frac_unplaced_pon_before'))} -> {fmt(r.get('frac_unplaced_pon_after'))} |")
    md += ["", "Details of the PASS classes (counts): i exact / on a branch base that differs from GRCh38 / INDEL composite; "
           "ii in an insertion / ALT = GRCh38 base / INDEL next to a branch INDEL / too long; iii no left / no right / "
           "inconsistent anchors", "",
           "| run | set | i exact | i branch base differs | i composite | ii insertion | ii ALT = GRCh38 | ii INDEL context | "
           "ii too long | iii no left | iii no right | iii inconsistent | class-i PASS tagged by the PoN | merged with a placed call |",
           "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows_p:
        md.append(f"| `{r['run']}` | {r['set']} | {r['i_exact']} | {r['i_branch_base_differs']} | {r['i_composite']} | "
                  f"{r['ii_in_insertion']} | {r['ii_alt_is_grch38']} | {r['ii_indel_context_gap']} | {r['ii_too_long']} | "
                  f"{r['iii_no_left']} | {r['iii_no_right']} | {r['iii_inconsistent']} | "
                  f"{fmt(r.get('projected_class_i_pass_pon_tagged'))} | {r['merged_with_placed']} |")
    for ev, title in (("bed_pon", "in the BED, after the PoN"), ("bed_raw", "in the BED, before the PoN"),
                      ("nobed_pon", "whole chr1 (no BED), after the PoN")):
        md += ["", f"## chr1 PASS P / R / F1 {title}: placed-only vs complete output", "",
               "| run | set | kind | graph | published (RIS) | placed (CPU) | placed + all unplaced as FP | "
               "+ class i (rtg) | + class ii by interval | + class iii as FP (worst) | class ii TP / calls | class iii |",
               "|---|---|---|---|---|---|---|---|---|---|---|---:|"]
        by = {}
        for r in rows_m:
            if r["eval"] == ev:
                by.setdefault((r["run"], r["vcf_dir"]), {})[r["variant"]] = r
        for (run, vdir), d in by.items():
            any_r = next(iter(d.values()))
            cell = lambda v: (f"{fmt(d[v]['precision'])} / {fmt(d[v]['recall'])} / **{fmt(d[v]['f1'])}**"  # noqa: E731
                              if v in d else "-")
            md.append(f"| `{run}` | {any_r['set']} | {any_r['kind']} | {any_r['graph']} | {cell('published')} | {cell('placed')} | "
                      f"{cell('placed_unplaced_fp')} | {cell('projected_rtg')} | {cell('projected_ii')} | "
                      f"{cell('projected_worst')} | {any_r['class_ii_tp']} / {any_r['class_ii_calls']} | "
                      f"{any_r['class_iii_calls']} |")
    if linear:
        md += ["", "## Linear-graph CPU runs (no off-reference nodes): their own chr1 reports, in the BED", "",
               "| run | kind | unplaced PASS | raw P / R / F1 | PoN P / R / F1 |", "|---|---|---:|---|---|"]
        for r in linear:
            md.append(f"| `{r['run']}` | {r['kind']} | {fmt(r['unplaced_pass'])} | "
                      f"{fmt(r['raw']['precision'])} / {fmt(r['raw']['recall'])} / {fmt(r['raw']['f1'])} | "
                      f"{fmt(r['pon']['precision'])} / {fmt(r['pon']['recall'])} / {fmt(r['pon']['f1'])} |")
    (out / "tables.md").write_text("\n".join(md) + "\n")
    # coordination TSVs: <scratch>/<RUN>.projected.tsv = all vcf dirs of the run
    for rd in sorted(SCRATCH.glob("runs/*")):
        tsvs = sorted(rd.glob("*/projected.tsv"))
        if not tsvs:
            continue
        tmp = SCRATCH / f"{rd.name}.projected.tsv.tmp"
        with open(tmp, "w") as f:
            for k, t in enumerate(tsvs):
                lines = t.read_text().splitlines(True)
                f.writelines(lines if k == 0 else lines[1:])
        tmp.rename(SCRATCH / f"{rd.name}.projected.tsv")
    print(f"{len(rows_p)} projected run/vcf dirs, {sum(1 for s in summary if s['metrics'])} with metrics, "
          f"{len(linear)} linear runs -> {out}")


if __name__ == "__main__":
    sys.exit(main())

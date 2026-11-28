#!/usr/bin/env python3
"""chr1 SNVs and INDELs by the number of the normal's germline variants within 50 bp (0 / 1 / >=2), HG008T, COLO829T
and HCC1395T, six-tumor models, seed 0, none with the population site features as a model input (user 2026-11-27:
the COLO829T and HCC1395T counterparts of Response Tables 2.7 and 2.8).

E19's stratification (experiments/E19_germline_density: germline_density.process, strata_latest.stats): PASS calls
after the population filter in the benchmark BED; group = germline records of the normal (labels.manifest germline
truth, any FILTER) that intersect the variant's normalized interval widened by 50 bp, the variant's own interval
excluded; recall by the truth variant's group, precision by the call's group. R = pansoma_net_v2_runs:
  GRCh38-only  R/m4m6_graphs_20261110/multicancer6_<P>_<KIND>_linear_nositef
  HPRC v1.1    R/m6_nositef_graphs_20261118/multicancer6_<P>_<KIND>_d9_nositef
  HPRC v2.1    R/m6_nositef_graphs_20261118/multicancer6_<P>_<KIND>_d46_nositef
  SNV: vcf_chr1_<SET> (PASS rule 2); INDEL: vcf_<rule>_chr1_<T>, rule = ppr8 (Illumina, PacBio HiFi / Fiber-seq), pp (ONT).
Second table per tumor and kind: the same, only variants inside the normal's germline confident BED (the manifest's
germline BED): outside it the germline truth has no records, so a variant there is in group 0 by construction.
Checks (asserted): totals of every call set = its eval_bed/report.json (pon, at_threshold); one truth split per tumor
and kind. Per-record rows are cached in ./cache.
    /wanglab/jshen/anaconda3/bin/python strata_nositef_three_tumors.py  -> strata_nositef_three_tumors.md, .json
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
E19 = HERE.parents[1] / "experiments/E19_germline_density"
sys.path.insert(0, str(E19))
import strata_latest as S  # noqa: E402

G, L = S.G, S.L
G.OUT = HERE / "cache"
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
ARMS = (("linear", "GRCh38-only"), ("d9", "HPRC v1.1"), ("d46", "HPRC v2.1"))
SL = ("0", "1", "2+")
HEAD = ["| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | "
        "ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |", "|---|---|---|---|---|---|---|"]


def call_set(kind, t, p, arm):
    run = (R / f"m4m6_graphs_20261110/multicancer6_{p}_{kind}_linear_nositef" if arm == "linear"
           else R / f"m6_nositef_graphs_20261118/multicancer6_{p}_{kind}_{arm}_nositef")
    pre = "vcf_chr1_" if kind == "SNV" else f"vcf_{'pp' if p == 'ONT' else 'ppr8'}_chr1_"
    vs = sorted(x for x in run.glob(f"{pre}{t}*") if L.complete(x))
    assert len(vs) == 1, (kind, t, p, arm, vs)
    v = vs[0]
    lin = sorted(v.glob(f"*.{kind}.linear.vcf.gz"))
    assert len(lin) == 1, v
    raw = lin[0].name[: -len(f".{kind}.linear.vcf.gz")]
    st = raw.replace(".tensors", "")  # the RIS call sets carry the tensor directory's name
    sample, platform, a, purity = L.set_parts(st[: -len("_snvaf05")] if st.endswith("_snvaf05") else st)
    assert sample == t and a == arm and purity == "P100", (st, t, arm)
    return L.CallSet(run=f"{run.parent.name}__{run.name}", vcf_dir=v, group="insample", variant="main", model=run.name,
                     model_sample="m6", model_platform=p, kind=kind, set=st, sample=t, platform=platform, arm=arm,
                     purity=purity, prefix=f"{raw}.{kind}")


def sgn(x):
    return "–" if x is None else f"{x:+.3f}".replace("-", "−")


def diff(a, b):
    return None if a is None or b is None else a - b


def table(t, res):
    """Markdown rows and the ΔF1 lists of the six HPRC rows."""
    out, ds = [], {}
    for p in S.PLATS:
        for i, (arm, gname) in enumerate(ARMS):
            st = res[p, arm]
            d = None if arm == "linear" else [diff(st[s]["f"], res[p, "linear"][s]["f"]) for s in SL]
            if d:
                ds[p, arm] = d
            out.append(f"| {S.pname(t, p) if i == 0 else ''} | {gname} | {S.tri([st[s]['fp'] for s in SL])} | {S.tri([st[s]['p'] for s in SL])} | "
                       f"{S.tri([st[s]['r'] for s in SL])} | {S.tri([st[s]['f'] for s in SL])} | {'–' if d is None else ' / '.join(sgn(x) for x in d)} |")
    return out, ds


def ranges(ds):
    v = list(ds.values())
    col = lambda j: [d[j] for d in v if d[j] is not None]  # noqa: E731
    pair = lambda j: sum(d[j] is not None and d[0] is not None and d[j] > d[0] for d in v)  # noqa: E731
    return dict(none=(min(col(0)), max(col(0))), one=(min(col(1)), max(col(1))), two=(min(col(2)), max(col(2))),
                two_gt_none=pair(2), one_gt_none=pair(1), two_positive=sum(x > 0 for x in col(2)), none_positive=sum(x > 0 for x in col(0)))


def range_lines(r):
    return [f"- ΔF1 over the six HPRC rows: none {sgn(r['none'][0])} to {sgn(r['none'][1])} (> 0 in {r['none_positive']} of 6); one "
            f"{sgn(r['one'][0])} to {sgn(r['one'][1])}; at least two {sgn(r['two'][0])} to {sgn(r['two'][1])} (> 0 in {r['two_positive']} of 6)",
            f"- rows with ΔF1 (≥2) > ΔF1 (none): {r['two_gt_none']} of 6; ΔF1 (one) > ΔF1 (none): {r['one_gt_none']} of 6"]


def main():
    G.OUT.mkdir(exist_ok=True)
    out = ["# chr1 SNVs and INDELs by nearby germline variants of the normal (0 / 1 / ≥2 within 50 bp), three tumors", "",
           "Generated by `strata_nositef_three_tumors.py`. chr1, benchmark regions, PASS calls after the population filter. "
           "Six-tumor models, seed 0, none with the population site features as a model input. ΔF1 = HPRC − GRCh38-only of "
           "the same platform, from the unrounded F1. \"inside the germline BED\" = only truth variants and calls inside the "
           "confident regions of the normal's germline truth.", ""]
    J = {}
    for kind in ("SNV", "INDEL"):
        for t in S.TUMORS:
            res, inb, share = {}, {}, {}
            for p in S.PLATS:
                for arm, _ in ARMS:
                    cs = call_set(kind, t, p, arm)
                    rows = G.process(cs)
                    rep = cs.report("eval_bed")["pon"]["at_threshold"]
                    a = S.stats(rows, "all")
                    assert (a["tpb"], a["truth"] - a["tpb"], a["tpc"], a["fp"]) == (rep["tp_baseline"], rep["fn"], rep["tp_call"], rep["fp"]), cs.id
                    res[p, arm] = {s: S.stats(rows, s) for s in SL + ("all",)}
                    gb = L.bed(L.manifest(cs.set, kind)["germline_bed"])
                    rin = [r for r in rows if gb.contains_point(int(r["pos"]) - 1)]
                    inb[p, arm] = {s: S.stats(rin, s) for s in SL + ("all",)}
                    fp0 = [r for r in rows if r["part"] == "fp" and int(r["n50"]) == 0]
                    share[p, arm] = (sum(not gb.contains_point(int(r["pos"]) - 1) for r in fp0), len(fp0))
            truth = {tuple(res[k][s]["truth"] for s in SL) for k in res}
            tin = {tuple(inb[k][s]["truth"] for s in SL) for k in inb}
            assert len(truth) == 1 and len(tin) == 1, (kind, t, truth, tin)
            truth, tin = truth.pop(), tin.pop()
            rows1, d1 = table(t, res)
            rows2, d2 = table(t, inb)
            r1, r2 = ranges(d1), ranges(d2)
            out += [f"## chr1 {kind}, {t} (truth 0 / 1 / ≥2: {' / '.join(map(str, truth))})", ""] + HEAD + rows1 + [""] + range_lines(r1)
            out += ["- group-0 false positives outside the normal's germline BED: "
                    + "; ".join(f"{S.pname(t, p)} " + " / ".join(f"{share[p, a][0]} of {share[p, a][1]}" for a, _ in ARMS) for p in S.PLATS)
                    + " (GRCh38-only / v1.1 / v2.1)", "",
                    f"### chr1 {kind}, {t}, inside the germline BED only (truth 0 / 1 / ≥2: {' / '.join(map(str, tin))})", ""] + HEAD + rows2 + [""] + range_lines(r2) + [""]
            J[f"{kind} {t}"] = dict(truth=truth, truth_in_germline_bed=tin, dF1={f"{S.pname(t, p)} {a}": d for (p, a), d in d1.items()},
                                    dF1_in_germline_bed={f"{S.pname(t, p)} {a}": d for (p, a), d in d2.items()}, ranges=r1, ranges_in_germline_bed=r2,
                                    group0_fp_outside_germline_bed={f"{S.pname(t, p)} {a}": share[p, a] for (p, a) in share},
                                    cells={f"{S.pname(t, p)} {a}": {s: res[p, a][s] for s in SL + ("all",)} for (p, a) in res})
    (HERE / "strata_nositef_three_tumors.md").write_text("\n".join(out) + "\n")
    (HERE / "strata_nositef_three_tumors.json").write_text(json.dumps(J, indent=1) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()

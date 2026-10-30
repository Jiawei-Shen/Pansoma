#!/usr/bin/env python3
"""E17a on the final INDEL models (R3.10, R1.6): chr1 INDEL results after the PoN, in the benchmark BED, by INS / DEL,
length and repeat context, for the final models of the INDEL pipeline and the two comparators.

    /wanglab/jshen/anaconda3/bin/python e17a_final_models.py

Final models (tmp/indel_gap_20261020/final_20261023; PASS = ppr8, the post-PoN validation recall 0.8 of the callable
positives, stored in the checkpoint): release = trained on chr2-22 of HCC1395T (chr20 also held out), HCC1954T, H1437T,
H2009T, HG008T and COLO829T; tumor-unseen = the four cell lines only. Seeds 0, 1, 2 per platform. Call sets:
pansoma_net_v2_runs/<model>_s<seed>_<Illumina|PacBio|ONT>_INDEL/vcf_ppr8_chr1_<set>_d46/eval_bed (rtg pon_calls,
report.json, truth.INDEL.pon.vcf.gz). Comparators: the INDEL rows of ../E16_bootstrap/comparators.tsv (E8 re-score,
PASS + our PoN). Classes and counting as in e17a_indel_strata.py (imported); seeds are pooled (counts summed), so a
recall is the seed mean and a precision is weighted by calls; the overall table (section 1) gives seed means of P and F1,
as analyze.py and latest_models.py do. "retained recall" = TPb / (truth - PoN-tagged truth), the
recall among truth INDELs the population filter keeps. FP in HP>=7 = call-side FPs whose own context is HP>=7.
Writes final_models.md, final_models.tsv (every call set x axis x class, per seed) and final_models.json next to
this script.
"""
import csv
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import e17a_indel_strata as E
from e17a_indel_strata import C, ratio

HERE = Path(__file__).resolve().parent
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
DATASETS = [("HG008T", p) for p in ("Illumina", "PacBio", "ONT")] + [("COLO829T", p) for p in ("Illumina", "fiberseq", "ONT")]
MODELS = (("release", "Pansoma release"), ("unseen", "Pansoma tumor-unseen"))
COMPARATORS = (("ClairS-TO", "ClairS-TO"), ("DeepSomatic-TO", "DeepSomatic-TO"))
SEEDS = (0, 1, 2)


def final_units():
    units = []
    for sample, p in DATASETS:
        plat = "PacBio" if p == "fiberseq" else p
        for model, _ in MODELS:
            for seed in SEEDS:
                run = f"{model}_s{seed}_{plat}_INDEL"
                d = RUNS / run / f"vcf_ppr8_chr1_{sample}_{p}_d46" / "eval_bed"
                w = d / "vcfeval" / "pon_calls"
                if not ((d / "report.json").exists() and C._rtg_complete(w)):
                    continue
                u = C.Unit(uid=f"{run}/vcf_ppr8_chr1_{sample}_{p}_d46", run=run, set=f"{sample}_{p}_d46", kind="INDEL",
                           role=model, variant=f"s{seed}", model=run)
                u.rtg["bed"], u.report["bed"] = w, json.loads((d / "report.json").read_text())
                tp = d / "truth.INDEL.pon.vcf.gz"
                u.truth_pon = tp if tp.exists() else None
                units.append(u)
    return units


def comparator_units():
    out = []
    for u in C.read_comparators(("INDEL",)):
        name = next((n for prefix, n in COMPARATORS if u.model.startswith(prefix)), None)
        if name and (u.sample, u.set.split("_", 1)[1]) in DATASETS:
            u.role, u.variant = name, u.model
            out.append(u)
    return out


def pooled(rows):
    """Summed [TPb, FN, TPc, FP, tagged] per (axis|class) over a group's call sets; tagged from the first (it is the
    truth's, the same for every seed)."""
    t = defaultdict(lambda: [0, 0, 0, 0, 0])
    for i, r in enumerate(rows):
        for k, v in r["strata"].items():
            for j in range(4):
                t[k][j] += v[j]
            if i == 0:
                t[k][4] = v[4]
    return t


def main():
    units = final_units() + comparator_units()
    classify = E.Classifier(C.fasta_path(units) if units else None)
    per, rows = [], []
    for u in units:
        t, skipped, has_pon = E.strata(u, classify)
        rep = u.report.get("bed") or {}
        at = (rep.get("pon") or {}).get("at_threshold")
        allc = t.get(("all", "all"), [0] * 5)
        check = None if not at else (allc[0], allc[1], allc[2], allc[3]) == (at["tp_baseline"], at["fn"], at["tp_call"], at["fp"])
        best = (rep.get("pon") or {}).get("best") or {}
        per.append(dict(uid=u.uid, group=u.role, seed=u.variant, sample=u.sample, dataset=f"{u.sample}_{u.set.split('_')[1]}",
                        threshold=rep.get("threshold"), best_f1=best.get("f1"), ceiling=(rep.get("pon") or {}).get("ceiling"),
                        skipped=skipped, truth_pon=has_pon, report_check=check,
                        strata={f"{a}|{c}": v for (a, c), v in t.items()}))
        for (axis, cls), (tpb, fn, tpc, fp, tag) in sorted(t.items()):
            rows.append(dict(uid=u.uid, group=u.role, seed=u.variant, dataset=per[-1]["dataset"], axis=axis, cls=cls,
                             tp_baseline=tpb, fn=fn, truth=tpb + fn, recall=ratio(tpb, tpb + fn), tp_call=tpc, fp=fp,
                             calls=tpc + fp, precision=ratio(tpc, tpc + fp), truth_pon_tagged=tag if has_pon else None))
    with open(HERE / "final_models.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    groups = {}
    for sample, p in DATASETS:
        ds = f"{sample}_{p}"
        for g in [m for m, _ in MODELS] + [n for _, n in COMPARATORS]:
            members = [o for o in per if o["dataset"] == ds and o["group"] == g]
            if members:
                prf = []
                for o in members:
                    tpb, fn, tpc, fp, _ = o["strata"].get("all|all", [0] * 5)
                    pp, rr = ratio(tpc, tpc + fp), ratio(tpb, tpb + fn)
                    prf.append((pp, rr, 2 * pp * rr / (pp + rr) if pp and rr else 0.0))
                groups[(ds, g)] = dict(n=len(members), strata=pooled(members), prf=prf,
                                       best_f1=[o["best_f1"] for o in members if o["best_f1"] is not None],
                                       ceiling=[o["ceiling"] for o in members if o["ceiling"] is not None],
                                       threshold=[o["threshold"] for o in members if o["threshold"] is not None])
    (HERE / "final_models.json").write_text(json.dumps(dict(
        generated=datetime.now().isoformat(timespec="seconds"), units=per,
        groups=[dict(dataset=ds, group=g, n=v["n"], best_f1=v["best_f1"], ceiling=v["ceiling"], threshold=v["threshold"],
                     strata={k: list(x) for k, x in v["strata"].items()}) for (ds, g), v in groups.items()]), indent=1) + "\n")
    (HERE / "final_models.md").write_text(tables(per, groups) + "\n")
    print(f"{len(units)} call sets ({sum(o['group'] in ('release', 'unseen') for o in per)} Pansoma); report check "
          f"{sum(1 for o in per if o['report_check'])} / {sum(1 for o in per if o['report_check'] is not None)}")


def tables(per, groups):
    f = C.fmt
    names = dict(MODELS)
    label = lambda g: names.get(g, g)  # noqa: E731
    lines = ["# E17a on the final INDEL models (chr1, benchmark BED, PASS after the PoN)", "",
             f"Generated by `e17a_final_models.py` ({datetime.now():%Y-%m-%d %H:%M}). Pansoma: the final models of the "
             "INDEL pipeline at their stored PASS threshold (ppr8), seeds pooled (n = seeds; recall = seed mean, "
             "precision weighted by calls). Comparators: their PASS + our PoN. Classes as in `tables.md` "
             "(`truth_labels.RepeatContext`; length = |len ALT - len REF|).", "",
             "## 1. Overall", "",
             "TP / FP / FN per call set, P and F1: seed means for Pansoma. retained recall = TP / (truth - PoN-tagged truth). "
             "best F1 = the best F1 obtained by optimizing the threshold on the test set (rtg curve after the PoN; seed "
             "mean). ceiling = the highest recall on that curve.", "",
             "| dataset | call set | n | TP / FP / FN | P | R | F1 | best F1 (test-tuned) | ceiling | retained recall | "
             "PoN-tagged truth (recall bound) |", "|---|---|---:|---|---:|---:|---:|---:|---:|---:|---|"]
    for (ds, g), v in groups.items():
        tpb, fn, tpc, fp, tag = v["strata"].get("all|all", [0] * 5)
        n = v["n"]
        mean = lambda xs: sum(xs) / len(xs) if xs else None  # noqa: E731
        p, r, f1 = (mean([x[i] for x in v["prf"] if x[i] is not None]) for i in range(3))
        truth = (tpb + fn) / n
        lines.append(f"| {ds} | {label(g)} | {n} | {tpb / n:.1f} / {fp / n:.1f} / {fn / n:.1f} | {f(p)} | {f(r)} | {f(f1)} | "
                     f"{f(mean(v['best_f1']))} | {f(mean(v['ceiling']))} | {f(ratio(tpb / n, truth - tag))} | "
                     f"{tag} of {truth:.0f} ({f(1 - tag / truth) if truth else '-'}) |")
    for axis, classes in E.AXES:
        present = [c for c in classes if any(f"{axis}|{c}" in v["strata"] for v in groups.values())]
        lines += ["", f"## {2 + [a for a, _ in E.AXES].index(axis)}. By {axis}", "",
                  "Recall cell: R (truth n; PoN-tagged). Precision cell: P (calls per seed).", "",
                  "| dataset | call set | " + " | ".join(f"{c} R" for c in present) + " | " +
                  " | ".join(f"{c} P" for c in present) + " |", "|---|---|" + "---:|" * (2 * len(present))]
        for (ds, g), v in groups.items():
            rc, pc = [], []
            for c in present:
                tpb, fn, tpc, fp, tag = v["strata"].get(f"{axis}|{c}", [0] * 5)
                rc.append(f"{f(ratio(tpb, tpb + fn), 2)} ({(tpb + fn) // v['n']}; {tag})")
                pc.append(f"{f(ratio(tpc, tpc + fp), 2)} ({(tpc + fp) / v['n']:.0f})")
            lines.append(f"| {ds} | {label(g)} | " + " | ".join(rc) + " | " + " | ".join(pc) + " |")
    key = [("DEL", "HP>=7"), ("INS", "HP>=7"), ("DEL", "HP4-6"), ("INS", "HP4-6"), ("DEL", "STR"), ("INS", "STR"),
           ("DEL", "none"), ("INS", "none")]
    lines += ["", "## 5. INS / DEL by repeat context (type x context)", "", "Cell: R (truth n) / P (calls per seed).", "",
              "| dataset | call set | " + " | ".join(f"{t} {c}" for t, c in key) + " |", "|---|---|" + "---|" * len(key)]
    for (ds, g), v in groups.items():
        cells = []
        for t, c in key:
            tpb, fn, tpc, fp, _ = v["strata"].get(f"type x context|{t} {c}", [0] * 5)
            cells.append(f"{f(ratio(tpb, tpb + fn), 2)} ({(tpb + fn) // v['n']}) / {f(ratio(tpc, tpc + fp), 2)} "
                         f"({(tpc + fp) / v['n']:.0f})")
        lines.append(f"| {ds} | {label(g)} | " + " | ".join(cells) + " |")
    lines += ["", "## 6. False positives in long homopolymers", "",
              "FP per seed whose call has context HP>=7, of all FP; precision of the HP>=7 calls and of the calls outside "
              "repeats (context none).", "",
              "| dataset | call set | FP in HP>=7 / all FP | share | P HP>=7 | P none |", "|---|---|---|---:|---:|---:|"]
    for (ds, g), v in groups.items():
        n = v["n"]
        hp = v["strata"].get("context|HP>=7", [0] * 5)
        none = v["strata"].get("context|none", [0] * 5)
        allc = v["strata"].get("all|all", [0] * 5)
        lines.append(f"| {ds} | {label(g)} | {hp[3] / n:.1f} / {allc[3] / n:.1f} | {f(ratio(hp[3], allc[3]))} | "
                     f"{f(ratio(hp[2], hp[2] + hp[3]))} | {f(ratio(none[2], none[2] + none[3]))} |")
    bad = [o["uid"] for o in per if o["report_check"] is False]
    lines += ["", "## 7. Checks", "",
              f"- Pansoma call sets: {sum(o['group'] in ('release', 'unseen') for o in per)} (release / tumor-unseen x "
              f"3 platforms x 2 tumors x seeds); comparator call sets: {sum(o['group'] not in ('release', 'unseen') for o in per)}.",
              f"- Overall counts equal report.json `pon.at_threshold`: {sum(1 for o in per if o['report_check'])} of "
              f"{sum(1 for o in per if o['report_check'] is not None)}" + (f"; mismatches: {', '.join(bad)}" if bad else "") + ".",
              "- Records not classified as INDEL (same-length or symbolic): " +
              (", ".join(f"{o['uid']} {o['skipped']}" for o in per if o["skipped"]) or "none") + "."]
    return "\n".join(lines)


if __name__ == "__main__":
    main()

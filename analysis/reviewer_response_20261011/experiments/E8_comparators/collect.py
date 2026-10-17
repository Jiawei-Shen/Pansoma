#!/usr/bin/env python3
"""E8 collect: ClairS-TO / DeepSomatic tumor-only next to Pansoma (chr1, the same truth, BED, PoN rule and rtg options).

    /wanglab/jshen/anaconda3/bin/python collect.py

Reads
  - comparators: /scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score/<caller>/<dataset>/<KIND>/report.json
    (score_comparator.py) and calls/<caller>/<dataset>/command.txt;
  - Pansoma: analysis/benchmark_7fc7396_20261012/linear_grch38_7fc7396_20261012/summary.json (HG008T / COLO829T, arms
    d9 = HPRC v1.1, d46 = HPRC v2.1, linear = GRCh38-only graph) and colo829t_purity_joint_7fc7396_20261012/summary.json
    (purity joint models); Pansoma's chr1 curves for recall at matched precision and its calls for SBS96 from the CPU
    re-predictions /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/<run>/vcf_chr1*/eval_bed/ when present.
Writes tables.md and summary.json next to this script, and sbs96.tsv (R3.5).
"""
import csv
import gzip
import json
import math
import os
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
E = Path(os.environ.get("E8_DIR", "/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators"))
BENCH = HERE.parents[2] / "benchmark_7fc7396_20261012"
CPU = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
# the user's existing tumor-only runs (re-scored 2026-10-16): ClairS-TO 0.4.2 on HG008T, 0.4.1 on COLO829T; DeepSomatic 1.9.0
CALLERS = (("clairs_to", "ClairS-TO 0.4.1 / 0.4.2"), ("deepsomatic", "DeepSomatic 1.9.0"))
ARMS = (("d9", "Pansoma HPRC v1.1 d9"), ("d46", "Pansoma HPRC v2.1 d46"), ("linear", "Pansoma GRCh38-only graph"))
KINDS = ("SNV", "INDEL")


def datasets():
    rows = []
    with open(HERE / "datasets.tsv") as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        for line in f:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def load(path):
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def pansoma_main():
    """{(sample, platform, kind): {arm: entry}} from the linear summary (it carries d9, d46 and linear)."""
    out = {}
    for e in load(BENCH / "linear_grch38_7fc7396_20261012" / "summary.json") or []:
        out[(e["sample"], e["platform"], e["kind"])] = {arm: e[arm] for arm, _ in ARMS if arm in e}
    return out


def pansoma_purity():
    """{(platform, kind, purity): (run, eval_bed dict)} of the joint models."""
    s = load(BENCH / "colo829t_purity_joint_7fc7396_20261012" / "summary.json") or {}
    out = {}
    for label, j in (s.get("joint") or {}).items():
        platform, kind = label.split()
        for purity, d in j["purity"].items():
            out[(platform, kind, purity)] = (j["run"], d["eval_bed"], d.get("set"))
    return out


def indel_new_pon_rule(main_p, purity_p):
    """INDEL rows under the current INDEL PoN rule (HP >= 8: AF >= 0.01, else 1e-4; 2026-10-16): the 7fc7396 summaries
    were scored with the old rule (HP >= 7: 0.01, else 1e-3); the CPU re-predictions' INDEL call sets were re-scored in
    place (old reports kept as eval_bed.pon_hp7_1e-3, tmp/reviewer_cpu_20261013/pon_hp8_20261016). Replace raw / pon of
    every INDEL entry that has such a re-scored report; mark the others old_rule (no local call set, e.g. COLO829T d9)."""
    def conv(block):
        at = block["at_threshold"]
        return dict(tp=at["tp_baseline"], fp=at["fp"], fn=at["fn"], p=at["precision"], r=at["recall"], f1=at["f1"],
                    ceiling=block["ceiling"], best=block["best"])

    def swap(ev, vcf_dir):
        rep = load(vcf_dir / "eval_bed" / "report.json")
        if rep is None or not (vcf_dir / "eval_bed.pon_hp7_1e-3").exists():
            ev["old_rule"] = True
            return
        ev["raw"], ev["pon"] = conv(rep["raw"]), conv(rep["pon"])
        ev["truth_pon_tagged"], ev["truth_records"] = rep["truth_pon_tags"]["tagged"], rep["truth_pon_tags"]["records"]
        ev["new_rule_report"] = str(vcf_dir / "eval_bed" / "report.json")

    for (sample, platform, kind), arms in main_p.items():
        if kind == "INDEL":
            for arm, e in arms.items():
                swap(e["eval_bed"], CPU / e["run"] / "vcf_chr1")
    for (platform, kind, purity), (run, ev, set_name) in purity_p.items():
        if kind == "INDEL":
            swap(ev, CPU / run / f"vcf_chr1_{set_name}")


def sample_platform(dataset):
    sample, platform = dataset.split("_", 1)
    return sample, platform


def fmt(x, digits=3):
    return "-" if x is None else f"{x:.{digits}f}"


def prf(d):
    """'P / R / F1 (TP/FP/FN)' of an at-threshold row (comparator: describe() at_threshold; Pansoma: summary keys)."""
    if d is None:
        return "-"
    if "tp_baseline" in d:
        return (f"{fmt(d['precision'])} / {fmt(d['recall'])} / **{fmt(d['f1'])}** "
                f"({d['tp_baseline']}/{d['fp']}/{d['fn']})")
    return f"{fmt(d['p'])} / {fmt(d['r'])} / **{fmt(d['f1'])}** ({d['tp']}/{d['fp']}/{d['fn']})"


def best(d):
    if not d:
        return "-"
    return f"{fmt(d['f1'])} ({fmt(d['precision'])}, {fmt(d['recall'])})"


def recall_at_precision(roc, precision):
    """Highest recall on a curve (rows of rtg's weighted ROC) at precision >= the given one; None if never reached."""
    hits = [r["recall"] for r in roc or () if r["precision"] >= precision - 1e-9]
    return max(hits) if hits else None


def pansoma_curve(run):
    """Pansoma's chr1 PoN curve (report.json of the CPU re-prediction's vcf_chr1/eval_bed), or None."""
    for path in sorted((CPU / run).glob("vcf_chr1*/eval_bed/report.json")) if run else ():
        r = load(path)
        if r and "pon" in r:
            return r["pon"].get("curve_roc"), path
    return None, None


# ---------------------------------------------------------------- SBS96 (R3.5)
COMP = str.maketrans("ACGT", "TGCA")


def sbs96_key(fasta, chrom, pos, ref, alt):
    ctx = fasta.fetch(chrom, pos - 2, pos + 1).upper()
    if len(ctx) != 3 or "N" in ctx or ctx[1] != ref:
        return None
    if ref in "GA":
        ctx, ref, alt = ctx.translate(COMP)[::-1], ref.translate(COMP), alt.translate(COMP)
    return f"{ctx[0]}[{ref}>{alt}]{ctx[2]}"


def sbs96_counts(vcfs, fasta):
    counts = Counter()
    for vcf in vcfs:
        if not Path(vcf).exists():
            return None
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split("\t", 6)
                ref = c[3].upper()
                for alt in c[4].upper().split(","):
                    if len(ref) == 1 and len(alt) == 1 and alt in "ACGT":
                        key = sbs96_key(fasta, c[0], int(c[1]), ref, alt)
                        if key:
                            counts[key] += 1
    return counts


def cosine(a, b):
    keys = set(a) | set(b)
    dot = sum(a.get(k, 0) * b.get(k, 0) for k in keys)
    na, nb = math.sqrt(sum(v * v for v in a.values())), math.sqrt(sum(v * v for v in b.values()))
    return dot / (na * nb) if na and nb else None



# ---------------------------------------------------------------- INDEL strata (R3.10)
LENGTH_BINS = ((1, 1, "1"), (2, 4, "2-4"), (5, 15, "5-15"), (16, 50, "16-50"), (51, 10 ** 9, ">50"))


def indel_class(repeats, chrom, pos, ref, alt):
    """(type, length bin, repeat class) of a VCF INDEL: shared prefix / suffix removed, RepeatContext of the
    label rules (indexed_gam_pipeline_v4 truth_labels: HP>=7, HP4-6, STR, VNTR, none, complex)."""
    ref, alt = ref.upper(), alt.upper()
    i = 0
    while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
        i += 1
    r, a = ref[i:], alt[i:]
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    pos0 = pos - 1 + i
    kind = "INS" if not r else "DEL" if not a else "complex"
    n = abs(len(ref) - len(alt))
    length = next(label for lo, hi, label in LENGTH_BINS if lo <= n <= hi)
    context = repeats.classify(chrom, pos0, r, a) if kind != "complex" else "complex"
    return kind, length, context


def indel_records(vcf):
    if not Path(vcf).exists():
        return None
    out = []
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                c = line.split("\t", 6)
                for alt in c[4].split(","):
                    if len(alt) != len(c[3]) and not alt.startswith("<") and alt not in ("*", "."):
                        out.append((c[0], int(c[1]), c[3], alt))
    return out


def indel_strata(work, repeats):
    """{(axis, value): [tp_baseline, fn, tp_call, fp]} of an rtg output directory (axis: type, length, context)."""
    files = {k: indel_records(work / f"{k}.vcf.gz") for k in ("tp-baseline", "fn", "tp", "fp")}
    if any(v is None for v in files.values()):
        return None
    table = {}
    for col, key in enumerate(("tp-baseline", "fn", "tp", "fp")):
        for rec in files[key]:
            for axis, value in zip(("type", "length", "context"), indel_class(repeats, *rec)):
                table.setdefault((axis, value), [0, 0, 0, 0])[col] += 1
    return table

# ---------------------------------------------------------------- main
def main():
    rows = datasets()
    main_p, purity_p = pansoma_main(), pansoma_purity()
    indel_new_pon_rule(main_p, purity_p)
    summary = dict(datasets=[], pending=[])
    lines = ["# E8: ClairS-TO and DeepSomatic tumor-only under Pansoma's protocol (chr1)", "",
             "Generated by `collect.py`. chr1 only; truth and BED of the tensor set's labels; rtg vcfeval "
             "`--squash-ploidy --sample ALT,ALT`; PoN = `scripts/filter_panel_of_normals.py` with the four PoN VCFs of "
             "`jobs/vcf.sh` (INDEL rule: HP >= 8 bp AF >= 0.01, else >= 1e-4). Pansoma rows: the 7fc7396 benchmarks (PASS = "
             "validation recall-0.9 threshold; 'raw' = before the PoN); INDEL rows from the CPU re-predictions re-scored "
             "with the current INDEL rule, '(old INDEL PoN rule)' where none is local. Comparator rows: `pass` = the tool's own PASS (its own filters and PoN), "
             "`pass_pon` = PASS + our PoN, curve = all non-RefCall records + our PoN ranked by the tool's QUAL. "
             "Cells: P / R / **F1** (TP/FP/FN).", ""]
    for kind in KINDS:
        lines += [f"## chr1 {kind}, in the benchmark BED", "",
                  "| dataset | caller | PASS before our PoN | PASS after our PoN | chr1 curve best F1 after PoN (P, R) | "
                  "curve ceiling after PoN | PoN bound on recall |", "|---|---|---|---|---|---|---:|"]
        for d in rows:
            name = d["dataset"]
            sample, platform = sample_platform(name)
            bound = None
            got = False
            for caller, label in CALLERS:
                rep = load(E / "score" / caller / name / kind / "report.json")
                if rep is None:
                    summary["pending"].append(f"{caller}/{name}/{kind}")
                    continue
                got = True
                t = rep["truth_pon_tags"]
                bound = 1 - t["tagged"] / t["records"] if t and t["records"] else None
                b = rep["eval_bed"]
                lines.append(f"| {name} | {label} | {prf(b['pass']['at_threshold'])} | "
                             f"{prf(b['pass_pon']['at_threshold'])} | {best(b['all_pon']['best'])} | "
                             f"{fmt(b['all_pon']['ceiling'])} | {fmt(bound)} |")
                summary["datasets"].append(dict(dataset=name, caller=caller, kind=kind, report=str(
                    E / "score" / caller / name / kind / "report.json"), truth_records=rep["truth_records"],
                    counts=rep["counts"], pon_bound=bound,
                    eval_bed={k: dict(at=v["at_threshold"], best=v["best"], ceiling=v["ceiling"],
                                      precision_at_recall=v["precision_at_recall"]) for k, v in b.items()},
                    eval_nobed={k: dict(at=v["at_threshold"], best=v["best"], ceiling=v["ceiling"])
                                for k, v in rep["eval_nobed"].items()}))
            arms = main_p.get((sample, platform, kind), {})
            purity = name.rsplit("_", 1)[-1] if name.rsplit("_", 1)[-1] in ("P75", "P50", "P25") else None
            if sample == "COLO829T" and purity is None:
                purity = "P100"
            if purity:
                key = (platform.split("_")[0], kind, purity)
                if key in purity_p:
                    run, ev, _ = purity_p[key]
                    lines.append(f"| {name} | Pansoma joint purity model (d9) `{run}`{' (old INDEL PoN rule)' if ev.get('old_rule') else ''} | {prf(ev['raw'])} | {prf(ev['pon'])} | "
                                 f"{best(ev['pon']['best'])} | {fmt(ev['pon']['ceiling'])} | - |")
            for arm, label in ARMS:
                if arm in arms:
                    ev = arms[arm]["eval_bed"]
                    lines.append(f"| {name} | {label} `{arms[arm]['run']}`{' (old INDEL PoN rule)' if ev.get('old_rule') else ''} | {prf(ev['raw'])} | {prf(ev['pon'])} | "
                                 f"{best(ev['pon']['best'])} | {fmt(ev['pon']['ceiling'])} | "
                                 f"{fmt(1 - ev['truth_pon_tagged'] / ev['truth_records']) if ev.get('truth_records') else '-'} |")
            if not got and not arms:
                continue
        lines.append("")

    # recall at matched precision (R1.6): comparator curve at Pansoma's PASS+PoN precision, and the reverse
    lines += ["## Recall at matched precision (chr1, BED, after our PoN)", "",
              "For each Pansoma arm at its PASS + PoN point (P, R): each comparator's highest recall on its all + PoN "
              "curve at precision >= P. Reverse: Pansoma's recall on its PoN curve at the comparator's PASS + PoN "
              "precision (needs Pansoma's curve from the CPU re-predictions; '-' when not yet there).", "",
              "| dataset | kind | Pansoma arm | Pansoma P / R | ClairS-TO R at that P | DeepSomatic R at that P | "
              "ClairS-TO PASS+PoN P / R -> Pansoma R | DeepSomatic PASS+PoN P / R -> Pansoma R |",
              "|---|---|---|---|---:|---:|---|---|"]
    matched = []
    for d in rows:
        name = d["dataset"]
        sample, platform = sample_platform(name)
        for kind in KINDS:
            reps = {c: load(E / "score" / c / name / kind / "report.json") for c, _ in CALLERS}
            if not any(reps.values()):
                continue
            for arm, label in ARMS:
                entry = main_p.get((sample, platform, kind), {}).get(arm)
                if not entry:
                    continue
                pp = entry["eval_bed"]["pon"]
                curve, curve_path = pansoma_curve(entry["run"])
                cells, rev = [], []
                for c, _ in CALLERS:
                    r = reps[c]
                    if r is None:
                        cells.append("-"), rev.append("-")
                        continue
                    cells.append(fmt(recall_at_precision(r["eval_bed"]["all_pon"]["curve_roc"], pp["p"])))
                    at = r["eval_bed"]["pass_pon"]["at_threshold"]
                    pr = recall_at_precision(curve, at["precision"]) if curve and at["precision"] is not None else None
                    rev.append(f"{fmt(at['precision'])} / {fmt(at['recall'])} -> {fmt(pr)}")
                    matched.append(dict(dataset=name, kind=kind, arm=arm, caller=c, pansoma_p=pp["p"], pansoma_r=pp["r"],
                                        caller_r_at_pansoma_p=recall_at_precision(r["eval_bed"]["all_pon"]["curve_roc"], pp["p"]),
                                        caller_p=at["precision"], caller_r=at["recall"], pansoma_r_at_caller_p=pr,
                                        pansoma_curve=str(curve_path) if curve_path else None))
                lines.append(f"| {name} | {kind} | {arm} | {fmt(pp['p'])} / {fmt(pp['r'])} | {cells[0]} | {cells[1]} | "
                             f"{rev[0]} | {rev[1]} |")
    summary["matched_precision"] = matched
    lines.append("")

    # SBS96 (R3.5): truth, comparator PASS + PoN calls (TP and FP), Pansoma PASS + PoN calls when available
    try:
        import pysam
        fasta = pysam.FastaFile(FASTA)
    except Exception as exc:  # noqa: BLE001
        fasta = None
        lines.append(f"SBS96 skipped: {exc}")
    sbs_rows = []
    if fasta:
        lines += ["## SBS96 cosine similarity to the chr1 truth spectrum (SNV, BED, after our PoN; R3.5)", "",
                  "| dataset | call set | SNVs | cosine (all calls) | TP | FP | cosine (FP only) |", "|---|---|---:|---:|---:|---:|---:|"]
        for d in rows:
            name = d["dataset"]
            sample, platform = sample_platform(name)
            truth = None
            sets = []
            for c, label in CALLERS:
                w = E / "score" / c / name / "SNV" / "eval_bed" / "pass_pon"
                if (w / "tp.vcf.gz").exists():
                    if truth is None:
                        truth = sbs96_counts([w / "tp-baseline.vcf.gz", w / "fn.vcf.gz"], fasta)
                    sets.append((label, w))
            for arm, label in ARMS:
                entry = main_p.get((sample, platform, "SNV"), {}).get(arm)
                if entry:
                    for w in sorted((CPU / entry["run"]).glob("vcf_chr1*/eval_bed/vcfeval/pon_calls")):
                        if (w / "tp.vcf.gz").exists():
                            sets.append((f"{label} `{entry['run']}`", w))
                            break
            if truth is None or not sets:
                continue
            sbs_rows.append(dict(dataset=name, call_set="truth", **truth))
            lines.append(f"| {name} | truth (chr1, BED) | {sum(truth.values())} | 1 | - | - | - |")
            for label, w in sets:
                tp, fp = sbs96_counts([w / "tp.vcf.gz"], fasta), sbs96_counts([w / "fp.vcf.gz"], fasta)
                if tp is None or fp is None:
                    continue
                allc = tp + fp
                sbs_rows.append(dict(dataset=name, call_set=label + " TP", **tp))
                sbs_rows.append(dict(dataset=name, call_set=label + " FP", **fp))
                lines.append(f"| {name} | {label} | {sum(allc.values())} | {fmt(cosine(allc, truth))} | "
                             f"{sum(tp.values())} | {sum(fp.values())} | {fmt(cosine(fp, truth))} |")
        lines.append("")
    if sbs_rows:
        keys = sorted({k for r in sbs_rows for k in r} - {"dataset", "call_set"})
        with open(HERE / "sbs96.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["dataset", "call_set", *keys])
            for r in sbs_rows:
                w.writerow([r["dataset"], r["call_set"], *[r.get(k, 0) for k in keys]])

    # INDEL strata (R3.10): recall by the truth's class, precision by the calls' class, PASS + our PoN, in the BED
    import sys
    sys.path.insert(0, str(HERE.parents[3]))
    try:
        from indexed_gam_pipeline_v4.tensor_postprocessing.truth_labels import RepeatContext
        repeats = RepeatContext(FASTA)
    except Exception as exc:  # noqa: BLE001
        repeats = None
        lines.append(f"INDEL strata skipped: {exc}")
    strata_rows = []
    if repeats:
        lines += ["## INDEL strata (chr1, BED, PASS after our PoN; R3.10)", "",
                  "Recall = TP / truth of the class (truth's class); precision = TP / calls of the class (calls' class). "
                  "Context: RepeatContext of the label rules (HP>=7, HP4-6, STR unit 2-6 >= 3 units, VNTR, none).", "",
                  "| dataset | call set | HP>=7 R (n) | HP>=7 P (calls) | HP4-6 R | HP4-6 P | STR R | STR P | none R | none P | "
                  "INS R / P | DEL R / P |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for d in rows:
            name = d["dataset"]
            sample, platform = sample_platform(name)
            sets = [(label, E / "score" / c / name / "INDEL" / "eval_bed" / "pass_pon") for c, label in CALLERS]
            for arm, label in ARMS:
                entry = main_p.get((sample, platform, "INDEL"), {}).get(arm)
                if entry:
                    found = sorted((CPU / entry["run"]).glob("vcf_chr1*/eval_bed/vcfeval/pon_calls"))
                    if found:
                        sets.append((f"{label} `{entry['run']}`", found[0]))
            for label, w in sets:
                t = indel_strata(w, repeats) if (w / "tp.vcf.gz").exists() else None
                if t is None:
                    continue
                for (axis, value), (tpb, fn, tpc, fp) in sorted(t.items()):
                    strata_rows.append([name, label, axis, value, tpb, fn, tpc, fp])

                def rp(axis, value):
                    tpb, fn, tpc, fp = t.get((axis, value), [0, 0, 0, 0])
                    r = f"{tpb / (tpb + fn):.2f} ({tpb + fn})" if tpb + fn else "- (0)"
                    p = f"{tpc / (tpc + fp):.2f} ({tpc + fp})" if tpc + fp else "- (0)"
                    return r, p
                cells = [c for v in ("HP>=7", "HP4-6", "STR", "none") for c in rp("context", v)]
                ins, dele = rp("type", "INS"), rp("type", "DEL")
                lines.append(f"| {name} | {label} | " + " | ".join(cells) + f" | {ins[0]} / {ins[1]} | {dele[0]} / {dele[1]} |")
        lines.append("")
    if strata_rows:
        with open(HERE / "indel_strata.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["dataset", "call_set", "axis", "class", "tp_baseline", "fn", "tp_call", "fp"])
            w.writerows(strata_rows)

    # commands (R3.m2)
    lines += ["## Commands (R3.m2)", "", "| caller | dataset | command |", "|---|---|---|"]
    for d in rows:
        for c, label in CALLERS:
            p = E / "calls" / c / d["dataset"] / "command.txt"
            if p.exists():
                lines.append(f"| {label} | {d['dataset']} | `{p.read_text().strip()}` |")
    lines += ["", f"Pending (no report yet): {len(summary['pending'])}: " + ", ".join(summary["pending"]), ""]
    (HERE / "tables.md").write_text("\n".join(lines) + "\n")
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""E17a (R3.10): chr1 INDEL results after the PoN, in the benchmark BED, by INS / DEL, length and repeat context.

    /wanglab/jshen/anaconda3/bin/python e17a_indel_strata.py

Input: the INDEL units of ../E16_bootstrap/callsets.py (rtg eval_bed/vcfeval/pon_calls: tp-baseline, fn = truth side;
tp, fp = call side), including the vcf_chr1_valbest variants (threshold chosen on validation) and comparators.tsv
rows of kind INDEL (slot). The PoN-tagged truth (eval_bed/truth.INDEL.pon.vcf.gz) gives the per-stratum recall bound.
Classes of a record (its first ALT):
  type     INS / DEL / complex: indexed_gam_pipeline_v4 truth_labels.split_alleles (the label pipeline's trimming,
           imported, not modified); complex = REF and ALT both keep bases after trimming
  length   |len(ALT) - len(REF)|: 1, 2-4, 5-10, >10
  context  truth_labels.RepeatContext.classify (imported, not modified) on the labels' GRCh38 FASTA: HP>=7 (unit 1,
           tract >= 7 bp), HP4-6, STR (unit 2-6, >= 3 units), VNTR (unit > 6, >= 2 units), none; complex records
           'complex'
Recall of a class = TPb / (TPb + FN) by the truth record's class; precision = TPc / (TPc + FP) by the call's class
(rtg may match a call and a truth record of different representation, so the two sides are classified separately).
Writes strata.tsv (every unit x axis x class, plus the type x context and type x length crosses), summary.json and
tables.md next to this script.
"""
import csv
import json
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "E16_bootstrap"))
import callsets as C  # noqa: E402

sys.path.insert(0, str(C.REPO))
from indexed_gam_pipeline_v4.tensor_postprocessing.truth_labels import RepeatContext, split_alleles  # noqa: E402

TYPES = ("INS", "DEL", "complex")
LENGTHS = ("1", "2-4", "5-10", ">10")
CONTEXTS = ("HP>=7", "HP4-6", "STR", "VNTR", "none", "complex")
AXES = (("type", TYPES), ("length", LENGTHS), ("context", CONTEXTS))


def length_bin(n):
    return "1" if n <= 1 else "2-4" if n <= 4 else "5-10" if n <= 10 else ">10"


class Classifier:
    def __init__(self, fasta):
        self.repeats = RepeatContext(fasta)
        self.cache = {}

    def __call__(self, chrom, pos, ref, alt):
        key = (chrom, pos, ref, alt)
        if key not in self.cache:
            self.cache[key] = self._classify(chrom, pos, ref, alt)
        return self.cache[key]

    def _classify(self, chrom, pos, ref, alt):
        n = abs(len(alt) - len(ref))
        for pos0, r, a, kind, _ in split_alleles(pos, ref, [alt]):
            if kind in ("INS", "DEL"):
                return kind, length_bin(n), self.repeats.classify(chrom, pos0, r, a)
            if kind == "COMPLEX":
                return "complex", length_bin(n), "complex"
        return None  # not an INDEL (SNP / MNP after trimming, symbolic)


def strata(unit, classify):
    """{(axis, class): [TPb, FN, TPc, FP, PoN-tagged truth]} of a unit's eval_bed pon_calls; and skipped records."""
    w = unit.rtg["bed"]
    tagged = None
    if unit.truth_pon:
        tagged = {(r[0], r[1], r[2], r[3]) for r in C.read_records(unit.truth_pon) if "PanelOfNormals" in r[4].split(";")}
    table = defaultdict(lambda: [0, 0, 0, 0, 0])
    skipped = 0
    for col, name in enumerate(("tp-baseline", "fn", "tp", "fp")):
        for r in C.read_records(w / f"{name}.vcf.gz"):
            if len(r[2]) == len(r[3]):
                skipped += 1
                continue
            c = classify(r[0], r[1], r[2], r[3])
            if c is None:
                skipped += 1
                continue
            t, ln, ctx = c
            keys = [("type", t), ("length", ln), ("context", ctx), ("type x context", f"{t} {ctx}"),
                    ("type x length", f"{t} {ln}"), ("all", "all")]
            for k in keys:
                table[k][col] += 1
                if col < 2 and tagged is not None and (r[0], r[1], r[2], r[3]) in tagged:
                    table[k][4] += 1
    return dict(table), skipped, tagged is not None


def ratio(a, b):
    return a / b if b else None


def main():
    units = C.discover(kinds=("INDEL",))
    if not units:
        raise SystemExit("no finished INDEL call sets")
    classify = Classifier(C.fasta_path(units))
    rows, out = [], []
    for u in units:
        t, skipped, has_pon = strata(u, classify)
        at = (u.report.get("bed") or {}).get("pon", {}).get("at_threshold")
        allc = t.get(("all", "all"), [0] * 5)
        check = None
        if at:
            check = (allc[0], allc[1], allc[2], allc[3]) == (at["tp_baseline"], at["fn"], at["tp_call"], at["fp"])
        out.append(dict(uid=u.uid, label=u.label, role=u.role, variant=u.variant, set=u.set, threshold=(
            u.report.get("bed") or {}).get("threshold"), skipped=skipped, truth_pon=has_pon, report_check=check,
            strata={f"{a}|{c}": v for (a, c), v in t.items()}))
        for (axis, cls), (tpb, fn, tpc, fp, tag) in sorted(t.items()):
            rows.append(dict(uid=u.uid, label=u.label, role=u.role, variant=u.variant, set=u.set, axis=axis,
                             cls=cls, tp_baseline=tpb, fn=fn, truth=tpb + fn, recall=ratio(tpb, tpb + fn),
                             tp_call=tpc, fp=fp, calls=tpc + fp, precision=ratio(tpc, tpc + fp),
                             truth_pon_tagged=tag if has_pon else None,
                             recall_bound=(1 - tag / (tpb + fn)) if has_pon and tpb + fn else None))
    with open(HERE / "strata.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    (HERE / "summary.json").write_text(json.dumps(dict(generated=datetime.now().isoformat(timespec="seconds"),
                                                       units=out), indent=1) + "\n")
    (HERE / "tables.md").write_text(tables(units, out) + "\n")
    print(f"{len(units)} INDEL call sets; report check {sum(1 for o in out if o['report_check'])} / "
          f"{sum(1 for o in out if o['report_check'] is not None)}")


def tables(units, out):
    f = C.fmt
    lines = ["# E17a: INDEL strata (chr1, benchmark BED, PASS after the PoN)", "",
             f"Generated by `e17a_indel_strata.py` ({datetime.now():%Y-%m-%d %H:%M}); {len(out)} INDEL call sets. "
             "Recall by the truth record's class (TPb / (TPb + FN)), precision by the call's class (TPc / (TPc + FP)). "
             "Repeat context: `truth_labels.RepeatContext` of indexed_gam_pipeline_v4 (HP>=7, HP4-6, STR unit 2-6 "
             ">= 3 units, VNTR unit > 6 >= 2 units, none). [valbest] = PASS threshold chosen on validation (chr2-22); "
             "the plain rows use the validation recall-0.9 threshold.", "",
             "## 1. Overall (all INDELs)", "",
             "| call set | role | threshold | TP / FP / FN | P | R | F1 | truth PoN-tagged (recall bound) |",
             "|---|---|---:|---|---:|---:|---:|---|"]
    for o in out:
        tpb, fn, tpc, fp, tag = o["strata"].get("all|all", [0] * 5)
        p, r = ratio(tpc, tpc + fp), ratio(tpb, tpb + fn)
        f1 = 2 * p * r / (p + r) if p and r else None
        bound = f"{tag} ({f(1 - tag / (tpb + fn))})" if o["truth_pon"] and tpb + fn else "-"
        lines.append(f"| {o['label']} | {o['role']} | {f(o['threshold'], 4)} | {tpb} / {fp} / {fn} | {f(p)} | {f(r)} | "
                     f"{f(f1)} | {bound} |")
    for axis, classes in AXES:
        present = [c for c in classes if any(f"{axis}|{c}" in o["strata"] for o in out)]
        lines += ["", f"## {2 + [a for a, _ in AXES].index(axis)}. By {axis}", "",
                  "Recall cell: R (truth n; PoN-tagged). Precision cell: P (calls n).", "",
                  "| call set | " + " | ".join(f"{c} R" for c in present) + " | " +
                  " | ".join(f"{c} P" for c in present) + " |", "|---|" + "---:|" * (2 * len(present))]
        for o in out:
            rc, pc = [], []
            for c in present:
                tpb, fn, tpc, fp, tag = o["strata"].get(f"{axis}|{c}", [0] * 5)
                rc.append(f"{f(ratio(tpb, tpb + fn), 2)} ({tpb + fn}; {tag if o['truth_pon'] else '-'})")
                pc.append(f"{f(ratio(tpc, tpc + fp), 2)} ({tpc + fp})")
            lines.append(f"| {o['label']} | " + " | ".join(rc) + " | " + " | ".join(pc) + " |")
    # INS / DEL in homopolymers (the R3.10 claim): type x context, recall and precision
    key = [("DEL", "HP>=7"), ("INS", "HP>=7"), ("DEL", "HP4-6"), ("INS", "HP4-6"), ("DEL", "STR"), ("INS", "STR"),
           ("DEL", "none"), ("INS", "none")]
    lines += ["", "## 5. INS / DEL by repeat context (type x context)", "",
              "Cell: R (truth n) / P (calls n).", "",
              "| call set | " + " | ".join(f"{t} {c}" for t, c in key) + " |", "|---|" + "---|" * len(key)]
    for o in out:
        cells = []
        for t, c in key:
            tpb, fn, tpc, fp, _ = o["strata"].get(f"type x context|{t} {c}", [0] * 5)
            cells.append(f"{f(ratio(tpb, tpb + fn), 2)} ({tpb + fn}) / {f(ratio(tpc, tpc + fp), 2)} ({tpc + fp})")
        lines.append(f"| {o['label']} | " + " | ".join(cells) + " |")
    comp = [o for o in out if o["role"] == "comparator"]
    bad = [o["uid"] for o in out if o["report_check"] is False]
    lines += ["", "## 6. Comparators and checks", "",
              (f"- {len(comp)} comparator call sets (role 'comparator')." if comp else
               "- Comparators: none yet (`../E16_bootstrap/comparators.tsv` has no active rows; add the re-scored "
               "ClairS-TO / DeepSomatic `eval_bed/pass_pon` directories and rerun `collect.sh`)."),
              f"- Overall counts equal report.json `pon.at_threshold`: "
              f"{sum(1 for o in out if o['report_check'])} of {sum(1 for o in out if o['report_check'] is not None)}"
              + (f"; mismatches: {', '.join(bad)}" if bad else "") + ".",
              "- Records not classified as INDEL (same-length or symbolic): " +
              (", ".join(f"{o['label']} {o['skipped']}" for o in out if o["skipped"]) or "none") + "."]
    return "\n".join(lines)


if __name__ == "__main__":
    main()

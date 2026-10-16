"""E19: does the graph help more where the patient's germline variants cluster? (R2.5)

For every finished call set (calls_lib.discover; main threshold), the rtg records of the PASS calls after the PoN in
the BED (<vcf_dir>/eval_bed/vcfeval/pon_calls: tp-baseline + fn = truth, tp + fp = calls) are stratified by the
number of germline records of the normal (labels.manifest germline truth, any FILTER) that intersect the variant's
normalized reference interval widened by W = 50 and 500 bp, not counting records that intersect the variant itself.
Strata 0 / 1 / 2+. Recall = TP-baseline / truth of the stratum, precision = TP calls / calls of the stratum.
Paired arm differences per sample / platform / kind: the truth set is the same in every arm, so recall differences
are paired per truth record (discordant counts, exact McNemar p).
Per-record rows: DATA/E19/<call set>.tsv.gz. Tables: tables.md, summary.json. usage: python germline_density.py
"""
import gzip
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "call_level_common"))
import calls_lib as L  # noqa: E402

OUT = L.DATA / "E19"
WINDOWS = (50, 500)
STRATA = ("0", "1", "2+")
PARTS = ("tp-baseline", "fn", "tp", "fp")
PAIRS = (("d9", "linear"), ("d46", "linear"), ("d46", "d9"))


def pfmt(p):
    return "-" if p is None else "<0.001" if p < 0.001 else f"{p:.3f}"


def stratum(n):
    return "0" if n == 0 else "1" if n == 1 else "2+"


def process(cs):
    out = OUT / f"{cs.id}.v2.tsv.gz"
    srcs = [cs.path("eval_bed", p) for p in PARTS]
    if out.exists() and all(out.stat().st_mtime > s.stat().st_mtime for s in srcs):
        return read_rows(out)
    m = L.manifest(cs.set, cs.kind)
    germ = L.germline(m["germline_vcf"], m["fasta"])
    seq = L.chrom_seq(m["fasta"])
    tk = L.truth_tensor_keys(cs.set, cs.kind)
    rows = []
    for part in PARTS:
        for rec in L.read_records(cs.path("eval_bed", part)):
            key = None
            for alt in rec.alts:
                key = L.normalize(seq, rec.pos, rec.ref, alt)
                if key:
                    break
            if key is None:
                continue
            a, b = L.interval(*key)
            ht = int(L.has_tensor(tk, rec)) if part in ("tp-baseline", "fn") else ""
            rows.append([part, rec.pos, rec.ref, ",".join(rec.alts), ht] + [germ.count_near(a, b, w) for w in WINDOWS])
    L.write_tsv_gz(out, ["part", "pos", "ref", "alt", "has_tensor"] + [f"n{w}" for w in WINDOWS], rows)
    return read_rows(out)


def read_rows(path):
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f]


def summarize(cs, rows):
    res = dict(id=cs.id, label=cs.label, group=cs.group, variant=cs.variant, sample=cs.sample, platform=cs.platform,
               kind=cs.kind, arm=cs.arm, purity=cs.purity, model=cs.model, windows={})
    for w in WINDOWS:
        d = {s: dict(tpb=0, fn=0, tpc=0, fp=0, truth_tensor=0, tpb_tensor=0) for s in STRATA}
        tot = dict(tpb=0, fn=0, tpc=0, fp=0, truth_tensor=0, tpb_tensor=0)
        for r in rows:
            s = stratum(int(r[f"n{w}"]))
            k = {"tp-baseline": "tpb", "fn": "fn", "tp": "tpc", "fp": "fp"}[r["part"]]
            for x in (d[s], tot):
                x[k] += 1
                if r["has_tensor"] == "1":
                    x["truth_tensor"] += 1
                    x["tpb_tensor"] += k == "tpb"
        for s, x in list(d.items()) + [("all", tot)]:
            p, rr, f = L.prf(x["tpb"], x["fn"], x["tpc"], x["fp"])
            x.update(precision=p, recall=rr, f1=f, truth=x["tpb"] + x["fn"], calls=x["tpc"] + x["fp"],
                     tensor_share=x["truth_tensor"] / (x["tpb"] + x["fn"]) if x["tpb"] + x["fn"] else None,
                     recall_with_tensor=x["tpb_tensor"] / x["truth_tensor"] if x["truth_tensor"] else None)
        d["all"] = tot
        res["windows"][str(w)] = d
    rep = cs.report("eval_bed")["pon"]["at_threshold"]
    res["rtg"] = dict(tp=rep["tp_baseline"], fn=rep["fn"], tp_call=rep["tp_call"], fp=rep["fp"])
    return res


def truth_status(rows, w):
    """{(pos, ref, alt): (stratum, found)} of the truth records."""
    out = {}
    for r in rows:
        if r["part"] in ("tp-baseline", "fn"):
            out[(r["pos"], r["ref"], r["alt"])] = (stratum(int(r[f"n{w}"])), r["part"] == "tp-baseline")
    return out


def compact_rows(summ, w):
    rows = []
    for s in summ:
        d = s["windows"][str(w)]
        rows.append([s["label"]] + [f"{d[x]['truth']} / {d[x]['fp']}" for x in STRATA]
                    + [d[x]["recall"] for x in STRATA] + [d[x]["precision"] for x in STRATA] + [d[x]["f1"] for x in STRATA]
                    + [d["all"]["f1"]] + [f"{L.fmt(d[x]['tensor_share'])} / {L.fmt(d[x]['recall_with_tensor'])}" for x in STRATA])
    return rows


def main():
    t0 = time.time()
    sets = sorted(L.discover(variants=("main", "valbest")), key=L.sort_key)
    only = os.environ.get("ONLY")
    if only:
        sets = [cs for cs in sets if only in cs.id]
    L.log(f"{len(sets)} call sets")
    summ, rows_by = [], {}
    for cs in sets:
        rows = process(cs)
        rows_by[cs.id] = rows
        summ.append(summarize(cs, rows))
        L.log(f"{cs.id}: {len(rows)} records")
    out = ["# E19: precision and recall by the density of the normal's germline variants (chr1)", "",
           f"Generated by `germline_density.py` ({time.strftime('%Y-%m-%d %H:%M')}). PASS calls at the run's threshold, "
           "after the PoN, in the BED (rtg pon_calls). Stratum = germline records of the normal within W bp of the "
           "variant (the variant's own interval excluded): 0 / 1 / 2+. Recall by the truth record's stratum, precision "
           "by the call's stratum. `tensor share / R with tensor`: share of the stratum's truth with an exact tensor in "
           "the call set's own tensors (labels), and recall among those (separates candidate generation from the "
           "model).", ""]
    head = ["call set", "truth / FP (0)", "truth / FP (1)", "truth / FP (2+)", "R 0", "R 1", "R 2+", "P 0", "P 1", "P 2+",
            "F1 0", "F1 1", "F1 2+", "F1 all", "tensor share / R with tensor (0)", "(1)", "(2+)"]
    sec = 1
    for w in WINDOWS:
        for title, pred in (("in-sample", lambda s: s["group"] == "insample"),
                            ("joint model, per purity", lambda s: s["group"] == "joint"),
                            ("cross-sample", lambda s: s["group"] == "xs")):
            rs = [s for s in summ if pred(s)]
            out += [f"## {sec}. W = {w} bp, {title}", ""]
            sec += 1
            out += (L.md_table(head, compact_rows(rs, w)) if rs else ["(none finished yet)"]) + [""]

    # paired arm differences
    groups = defaultdict(dict)
    for s in summ:
        if s["group"] == "insample" and s["variant"] == "main":
            groups[(s["kind"], s["sample"], s["platform"])][s["arm"]] = s
    paired = []
    for w in WINDOWS:
        out += [f"## {sec}. Paired arm differences, W = {w} bp (in-sample, main threshold)", "",
                "`dR`, `dP`, `dF1` = arm A minus arm B per stratum. Recall discordance: truth records found by A only / "
                "by B only, with the exact two-sided McNemar p.", ""]
        sec += 1
        prow = []
        for (kind, sample, platform), arms in sorted(groups.items(), key=lambda x: (L.KIND_ORDER[x[0][0]], x[0][1], L.PLATFORM_ORDER.get(x[0][2], 9))):
            for a, b in PAIRS:
                if a not in arms or b not in arms:
                    continue
                ta = truth_status(rows_by[arms[a]["id"]], w)
                tb = truth_status(rows_by[arms[b]["id"]], w)
                for st in STRATA + ("all",):
                    da, db = arms[a]["windows"][str(w)][st], arms[b]["windows"][str(w)][st]
                    keys = [k for k in ta if k in tb and (st == "all" or ta[k][0] == st)]
                    only_a = sum(ta[k][1] and not tb[k][1] for k in keys)
                    only_b = sum(tb[k][1] and not ta[k][1] for k in keys)
                    diff = lambda x, y: None if x is None or y is None else x - y  # noqa: E731
                    rec = dict(window=w, kind=kind, sample=sample, platform=platform, a=a, b=b, stratum=st,
                               truth=len(keys), dR=diff(da["recall"], db["recall"]), dP=diff(da["precision"], db["precision"]),
                               dF1=diff(da["f1"], db["f1"]), fp_a=da["fp"], fp_b=db["fp"], only_a=only_a, only_b=only_b,
                               mcnemar_p=L.binom_two_sided(only_a, only_b))
                    paired.append(rec)
                    prow.append([f"{sample} {platform} {kind}", f"{a} - {b}", st, rec["truth"], rec["dR"], rec["dP"],
                                 rec["dF1"], f"{rec['fp_a']} / {rec['fp_b']}", f"{only_a} / {only_b}",
                                 pfmt(rec["mcnemar_p"])])
        out += (L.md_table(["set", "A - B", "stratum", "truth", "dR", "dP", "dF1", "FP A / B", "found A only / B only",
                            "McNemar p"], prow) if prow else
                ["No sample / platform / kind has two arms finished yet."]) + [""]
    (HERE / "tables.md").write_text("\n".join(out) + "\n")
    (HERE / "summary.json").write_text(json.dumps(dict(generated=time.strftime("%Y-%m-%d %H:%M"), call_sets=summ,
                                                       paired=paired), indent=1) + "\n")
    L.log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

"""Tables of the RIS GPU jobs of submit.sh / submit_e13.sh, from whatever has finished. Standard library only (login node):

    PANSOMA_RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs python3 collect.py

chr1 after the PoN, in the benchmark BED (COLO829T: the full SMaHT union, i.e. every base). SNV: PASS F1 at the
validation recall-0.9 threshold (best F1 of the chr1 curve beside it); INDEL: the best F1 of the chr1 curve (PASS F1
beside it), with the 71cb51a INDEL PoN rule (<run>/vcf_chr1_pon71cb51a, ../indel_pon_71cb51a_20261016). Writes
tables.md and summary.json next to this script.
"""
import json
import os
import statistics
from pathlib import Path

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
HERE = Path(__file__).resolve().parent
SETS = [("HG008", "Illumina"), ("HG008", "PacBio"), ("HG008", "ONT"),
        ("COLO829", "Illumina"), ("COLO829", "fiberseq"), ("COLO829", "ONT")]


def score(run, kind):
    """{pass, best, p, r, fp} of a run, or None when it is not evaluated yet."""
    d = RUNS / run
    f = d / ("vcf_chr1_pon71cb51a" if kind == "INDEL" else "vcf_chr1") / "eval_bed" / "report.json"
    if not f.exists():
        return None
    p = json.load(open(f))["pon"]
    t, b = p["at_threshold"], p["best"]
    return {"pass": t["f1"], "best": b["f1"], "p": t["precision"], "r": t["recall"], "fp": t["fp"],
            "best_tp": b["tp_baseline"]}


def main_value(s, kind):
    return None if s is None else (s["pass"] if kind == "SNV" else s["best"])


def f3(x):
    return "-" if x is None else f"{x:.3f}"


def d3(a, b):
    return "-" if a is None or b is None else f"{b - a:+.3f}"


def d9_run(s, p, kind, matched=True):
    if kind == "SNV":
        if matched or (s, p) == ("HG008", "Illumina"):
            return f"{s}_{p}_SNV_7fc7396_ddp4_b1024"
        return f"{s}_{p}_SNV_7fc7396" if s == "HG008" else f"{s}_{p}_SNV_h90"
    return f"{s}_{p}_INDEL_7fc7396"


def nopc_d9_run(s, p, kind):
    if kind == "SNV" and (s, p) == ("HG008", "Illumina"):
        return "HG008_Illumina_SNV_nopc_7fc7396_ddp4_b1024"
    return f"{s}_{p}_{kind}_nopc_7fc7396"


def main():
    out, L = {}, ["# RIS GPU jobs of the 2026-10-15 reviewer response: results so far\n",
                  "(`collect.py`; chr1 after the PoN, in the benchmark BED; SNV: PASS F1 (best F1), INDEL: best F1 (PASS F1) "
                  "with the 71cb51a INDEL PoN rule; one seed each unless stated; '-' = not finished.)\n"]

    # matched d9
    L += ["\n## Matched d9 (4 GPUs for SNV, as d46 / linear) against the earlier d9 runs and d46\n",
          "| set | kind | d9 earlier | d9 matched | d46 | d46 - d9 matched |", "|---|---|---|---|---|---|"]
    rows = []
    for s, p in SETS:
        for k in ("SNV", "INDEL"):
            early = d9_run(s, p, k, matched=False) if k == "SNV" else f"{s}_{p}_INDEL_small" if s == "COLO829" else None
            if k == "INDEL" and s == "HG008":
                early = None
            e, m, d = (score(early, k) if early else None), score(d9_run(s, p, k), k), score(f"{s}_{p}_{k}_d46_7fc7396", k)
            if k == "SNV" and (s, p) == ("HG008", "Illumina"):
                e = None  # the matched run is the earlier one
            rows.append(dict(set=f"{s}T {p}", kind=k, earlier=early, d9_earlier=e, d9_matched=m, d46=d))
            v = lambda x: main_value(x, k)  # noqa: E731
            L.append(f"| {s}T {p} | {k} | {f3(v(e))} | {f3(v(m))} | {f3(v(d))} | {d3(v(m), v(d))} |")
    out["matched_d9"] = rows

    # no ch6
    L += ["\n## No ch6 (`--drop-planes path_count`) against with ch6\n",
          "| set | kind | d9 with ch6 | d9 no ch6 | d9 change | d46 with ch6 | d46 no ch6 | d46 change |",
          "|---|---|---|---|---|---|---|---|"]
    rows = []
    for k in ("SNV", "INDEL"):
        for s, p in SETS:
            w9 = score(d9_run(s, p, k, matched=False), k)
            n9 = score(nopc_d9_run(s, p, k), k)
            w46, n46 = score(f"{s}_{p}_{k}_d46_7fc7396", k), score(f"{s}_{p}_{k}_nopc_d46_7fc7396", k)
            rows.append(dict(set=f"{s}T {p}", kind=k, d9_with=w9, d9_without=n9, d46_with=w46, d46_without=n46))
            v = lambda x: main_value(x, k)  # noqa: E731
            L.append(f"| {s}T {p} | {k} | {f3(v(w9))} | {f3(v(n9))} | {d3(v(w9), v(n9))} | {f3(v(w46))} | "
                     f"{f3(v(n46))} | {d3(v(w46), v(n46))} |")
    out["no_ch6"] = rows
    L.append("\nd9 controls: the runs on the GPU count of the no-ch6 run (SNV: HG008T Illumina "
             "`HG008_Illumina_SNV_7fc7396_ddp4_b1024`, 4 GPUs; the others `HG008_*_SNV_7fc7396` / `COLO829_*_SNV_h90`, 1 GPU; "
             "INDEL `*_INDEL_7fc7396`, 1 GPU).")

    # norm dataset
    L += ["\n## `--norm dataset` against `--norm train` (d46, in-sample)\n",
          "| set | kind | --norm train | --norm dataset | change |", "|---|---|---|---|---|"]
    rows = []
    for s, p, k in [(s, p, "SNV") for s, p in SETS] + [("HG008", "Illumina", "INDEL"), ("COLO829", "Illumina", "INDEL")]:
        a, b = score(f"{s}_{p}_{k}_d46_7fc7396", k), score(f"{s}_{p}_{k}_d46_normds", k)
        rows.append(dict(set=f"{s}T {p}", kind=k, train=a, dataset=b))
        v = lambda x: main_value(x, k)  # noqa: E731
        L.append(f"| {s}T {p} | {k} | {f3(v(a))} | {f3(v(b))} | {d3(v(a), v(b))} |")
    out["norm_dataset"] = rows

    # seeds
    L += ["\n## Seeds (SNV PASS F1 (best F1); seed 0 = the existing run, seeds 1 and 2 new)\n",
          "| set | arm | seed 0 | seed 1 | seed 2 | PASS mean (range) | best mean (range) |", "|---|---|---|---|---|---|---|"]
    rows = []
    for s, p in SETS:
        for arm, base in (("linear", f"{s}_{p}_SNV_linear_7fc7396"), ("d9", f"{s}_{p}_SNV_7fc7396_ddp4_b1024"),
                          ("d46", f"{s}_{p}_SNV_d46_7fc7396")):
            sc = [score(r, "SNV") for r in (base, base + "_s1", base + "_s2")]
            rows.append(dict(set=f"{s}T {p}", arm=arm, scores=sc))
            cells = ["-" if x is None else f"{x['pass']:.3f} ({x['best']:.3f})" for x in sc]
            summ = []
            for key in ("pass", "best"):
                got = [x[key] for x in sc if x is not None]
                summ.append(f"{statistics.mean(got):.3f} ({min(got):.3f}-{max(got):.3f})" if len(got) == 3 else "-")
            L.append(f"| {s}T {p} | {arm} | " + " | ".join(cells) + f" | {summ[0]} | {summ[1]} |")
    out["seeds"] = rows

    # E13
    L += ["\n## E13 on HG008T PacBio d46 (SNV PASS F1 (best F1))\n",
          "| run | training blocks / model | PASS F1 | best F1 |", "|---|---|---|---|"]
    rows = []
    e13 = [("HG008_PacBio_SNV_d46_7fc7396", "100 %, 50.4 M")]
    e13 += [(f"HG008_PacBio_SNV_d46_e13_keep{f}_rep{r}", f"{f} %, rep {r}") for f in (10, 25, 50, 75) for r in (1, 2)]
    e13 += [("HG008_PacBio_SNV_d46_e13_8p7M", "100 %, 8.7 M"), ("HG008_PacBio_SNV_d46_e13_199p6M", "100 %, 199.6 M")]
    for run, what in e13:
        x = score(run, "SNV")
        rows.append(dict(run=run, what=what, score=x))
        L.append(f"| {run} | {what} | {f3(x and x['pass'])} | {f3(x and x['best'])} |")
    out["e13"] = rows

    (HERE / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

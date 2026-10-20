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


PULL = 1792191762  # 2026-10-16 18:02:42 -0500: the 71cb51a INDEL PoN rule in the RIS checkout


def score(run, kind, vcf="vcf_chr1"):
    """{pass, best, p, r, fp} of a run, or None when it is not evaluated yet. INDEL: the 71cb51a PoN rule, i.e.
    <vcf>_pon71cb51a (../indel_pon_71cb51a_20261016), or <vcf> itself when its vcf.sh ran after that rule came in."""
    d = RUNS / run
    f = d / vcf / "eval_bed" / "report.json"
    if kind == "INDEL":
        f2 = d / f"{vcf}_pon71cb51a" / "eval_bed" / "report.json"
        f = f2 if f2.exists() else (f if f.exists() and f.stat().st_mtime > PULL else f2)
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
          "normds: BQ and MAPQ per dataset, predict fitted on chr1 only (2026-10-16); normds re-predicted: the same "
          "checkpoints with the checkpoint's statistics (predict of 2026-10-19); normbq: BQ only per dataset, "
          "statistics of all chromosomes (2026-10-19). Seed range: the three `--norm train` seeds.\n",
          "| set | kind | --norm train (seed range) | normds | normds re-predicted | normbq | normbq - train |",
          "|---|---|---|---|---|---|---|"]
    rows = []
    for s, p, k in [(s, p, "SNV") for s, p in SETS] + [("HG008", "Illumina", "INDEL"), ("COLO829", "Illumina", "INDEL")]:
        base = f"{s}_{p}_{k}_d46_7fc7396"
        a, b = score(base, k), score(f"{s}_{p}_{k}_d46_normds", k)
        c = score(f"{s}_{p}_{k}_d46_normds", k, vcf=f"vcf_chr1_{s}T_{p}_d46")
        q = score(f"{s}_{p}_{k}_d46_normbq", k)
        v = lambda x: main_value(x, k)  # noqa: E731
        seeds = [v(score(base + x, k)) for x in ("", "_s1", "_s2")] if k == "SNV" else []
        rng = f" ({min(seeds):.3f}-{max(seeds):.3f})" if seeds and None not in seeds else ""
        rows.append(dict(set=f"{s}T {p}", kind=k, train=a, dataset=b, dataset_repredicted=c, normbq=q))
        L.append(f"| {s}T {p} | {k} | {f3(v(a))}{rng} | {f3(v(b))} | {f3(v(c))} | {f3(v(q))} | {d3(v(a), v(q))} |")
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

    # --norm dataset (f344bca, BQ only per dataset): the reported models, three seeds for SNV
    def three(base, kind):
        sc = [score(base + x, kind) for x in ("", "_s1", "_s2")]
        got = [x for x in sc if x is not None]
        summ = {k: ((statistics.mean([x[k] for x in got]), min(x[k] for x in got), max(x[k] for x in got))
                    if len(got) == 3 else None) for k in ("pass", "best")}
        return sc, summ

    def ms(t):
        return "-" if t is None else f"{t[0]:.3f} ({t[1]:.3f}-{t[2]:.3f})"

    L += ["\n## `--norm dataset` (BQ only per dataset, f344bca): SNV over three seeds\n",
          "PASS F1 (best F1) per seed; mean (range). d46 seed 0 = `<S>_<P>_SNV_d46_normbq` (submit_normbq.sh), the rest "
          "submit_normbq_all.sh.\n",
          "| set | arm | seed 0 | seed 1 | seed 2 | PASS mean (range) | best mean (range) |", "|---|---|---|---|---|---|---|"]
    rows, means = [], {}
    for s, p in SETS:
        for arm in ("linear", "d9", "d46"):
            sc, summ = three(f"{s}_{p}_SNV_{arm}_normbq", "SNV")
            means[(s, p, arm)] = summ
            rows.append(dict(set=f"{s}T {p}", arm=arm, scores=sc))
            cells = ["-" if x is None else f"{x['pass']:.3f} ({x['best']:.3f})" for x in sc]
            L.append(f"| {s}T {p} | {arm} | " + " | ".join(cells) + f" | {ms(summ['pass'])} | {ms(summ['best'])} |")
    out["normbq_seeds"] = rows
    L += ["\nArm differences of the three-seed means (PASS F1 / best F1):\n",
          "| set | d46 - d9 | d9 - linear | d46 - linear |", "|---|---|---|---|"]
    for s, p in SETS:
        def dd(a, b):
            x, y = means[(s, p, a)], means[(s, p, b)]
            return "-" if not (x["pass"] and y["pass"]) else \
                f"{x['pass'][0] - y['pass'][0]:+.3f} / {x['best'][0] - y['best'][0]:+.3f}"
        L.append(f"| {s}T {p} | {dd('d46', 'd9')} | {dd('d9', 'linear')} | {dd('d46', 'linear')} |")

    L += ["\n## `--norm dataset`: INDEL per arm (best F1 (PASS F1), one seed)\n",
          "| set | linear | d9 | d46 |", "|---|---|---|---|"]
    rows = []
    for s, p in SETS:
        sc = {arm: score(f"{s}_{p}_INDEL_{arm}_normbq", "INDEL") for arm in ("linear", "d9", "d46")}
        rows.append(dict(set=f"{s}T {p}", scores=sc))
        L.append(f"| {s}T {p} | " + " | ".join("-" if x is None else f"{x['best']:.3f} ({x['pass']:.3f})"
                                              for x in sc.values()) + " |")
    out["normbq_indel"] = rows

    L += ["\n## `--norm dataset`: no ch6 on d46 (SNV three seeds each side; INDEL one seed)\n",
          "| set | kind | with ch6: mean (range) | no ch6: mean (range) | change of the means | no ch6 below every with-ch6 seed |",
          "|---|---|---|---|---|---|"]
    rows = []
    for k in ("SNV", "INDEL"):
        for s, p in SETS:
            if k == "SNV":
                (wsc, w), (nsc, n) = three(f"{s}_{p}_SNV_d46_normbq", "SNV"), three(f"{s}_{p}_SNV_nopc_d46_normbq", "SNV")
                key = "pass"
                below = "-" if not (w[key] and n[key]) else ("yes" if n[key][2] < w[key][1] else "no")
                L.append(f"| {s}T {p} | SNV (PASS) | {ms(w[key])} | {ms(n[key])} | "
                         f"{'-' if not (w[key] and n[key]) else f'{n[key][0] - w[key][0]:+.3f}'} | {below} |")
                rows.append(dict(set=f"{s}T {p}", kind=k, with_ch6=wsc, no_ch6=nsc))
            else:
                w, n = score(f"{s}_{p}_INDEL_d46_normbq", k), score(f"{s}_{p}_INDEL_nopc_d46_normbq", k)
                L.append(f"| {s}T {p} | INDEL (best) | {f3(w and w['best'])} | {f3(n and n['best'])} | "
                         f"{d3(w and w['best'], n and n['best'])} | - |")
                rows.append(dict(set=f"{s}T {p}", kind=k, with_ch6=w, no_ch6=n))
    out["normbq_no_ch6"] = rows

    (HERE / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

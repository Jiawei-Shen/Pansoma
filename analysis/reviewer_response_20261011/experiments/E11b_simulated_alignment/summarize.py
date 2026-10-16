#!/usr/bin/env python3
"""E11b step 5: tables from the evaluate.py outputs -> results/tables.md, results/summary.tsv, results/summary.json"""
import argparse
import json
import os

import numpy as np

ARMS = ["linear", "d9", "d46", "chm13"]
ARM_LABEL = {"linear": "GRCh38-only graph (H 17)", "d9": "HPRC v1.1 d9", "d46": "HPRC v2.1 d46",
             "chm13": "CHM13-only graph (control; CHM13 truth)"}
PLATS = ["illumina", "hifi", "ont"]


def pct(x, nd=3):
    return "n/a" if x is None else f"{100 * x:.{nd}f}%"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evals", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chm13-truth", default="/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment/"
                    "chm13_truth")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    res = {}
    for p in PLATS:
        for arm in ARMS:
            f = f"{a.evals}/{arm}.{p}.json"
            if os.path.exists(f):
                res[(arm, p)] = json.load(open(f))
    L = ["# E11b tables (simulated reads from the HG008-N v6.2 assembly, chr1 + chr20 one-to-one segments)", ""]
    L.append("correct = surjected primary alignment on the true chromosome with start or end within 50 bp of the true "
             "GRCh38 interval; wrong = placed on GRCh38 but not correct; perfect = graph alignment (GAF) covers the whole "
             "read with no mismatch / indel; MAPQ = giraffe's graph MAPQ.")
    L.append("")
    tsv = ["platform\tarm\treads\tunmapped\tgraph_only\tcorrect50\twrong\twrong_of_placed\tperfect\tperfect_wrong\t"
           "wrong_mapq_ge30\twrong_mapq60\tmapq0\toffref_reads"]
    for p in PLATS:
        rows = [(arm, res[(arm, p)]) for arm in ARMS if (arm, p) in res]
        if not rows:
            continue
        L.append(f"## {p}: overall")
        L.append("")
        L.append("| arm | reads | unmapped | graph-mapped, not surjected | correct (50 bp) | wrong | wrong / placed | "
                 "perfect | perfect but wrong | wrong at MAPQ>=30 | wrong at MAPQ 60 | MAPQ 0 | uses off-ref node |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for arm, r in rows:
            o = r["overall"]
            L.append(f"| {ARM_LABEL[arm]} | {o['reads']:,} | {o['unmapped']:,} ({pct(o['unmapped_frac'])}) | "
                     f"{o['graph_mapped_not_placed']:,} | {pct(o['correct50_frac'])} | {o['wrong']:,} | "
                     f"{pct(o['wrong_frac_of_placed'])} | {pct(o['perfect_frac'], 2)} | {o['perfect_wrong']:,} "
                     f"({pct(o['perfect_wrong_frac_of_perfect'])} of perfect) | {o['wrong_mapq_ge30']:,} "
                     f"({pct(o['wrong_frac_mapq_ge30'])}) | {o['wrong_mapq60']:,} ({pct(o['wrong_frac_mapq60'])}) | "
                     f"{pct(o['mapq0_frac'], 2)} | {o['offref_reads']:,} ({pct(o['offref_frac'], 2)}) |")
            tsv.append("\t".join(map(str, [p, arm, o["reads"], o["unmapped"], o["graph_mapped_not_placed"],
                                           o["correct50"], o["wrong"], o["wrong_frac_of_placed"], o["perfect"],
                                           o["perfect_wrong"], o["wrong_mapq_ge30"], o["wrong_mapq60"],
                                           o["mapq0_frac"], o["offref_reads"]])))
        L.append("")
        if p == "illumina":
            L.append("Error-free reads only (no simulated sequencing error; any graph mismatch is then a real "
                     "difference between the read and the graph path it took):")
            L.append("")
            L.append("| arm | reads | perfect | perfect but wrong | wrong / placed | wrong at MAPQ>=30 |")
            L.append("|---|---|---|---|---|---|")
            for arm, r in rows:
                o = r["error_free_reads"]
                L.append(f"| {ARM_LABEL[arm]} | {o['reads']:,} | {pct(o['perfect_frac'], 2)} | {o['perfect_wrong']:,} "
                         f"({pct(o['perfect_wrong_frac_of_perfect'])}) | {pct(o['wrong_frac_of_placed'])} | "
                         f"{o['wrong_mapq_ge30']:,} |")
            L.append("")
        L.append(f"## {p}: MAPQ calibration (reads placed on GRCh38; error rate = wrong / placed in the bin; "
                 "expected = mean 10^(-MAPQ/10))")
        L.append("")
        L.append("| MAPQ | " + " | ".join(f"{ARM_LABEL[arm]}: placed / wrong / rate / expected" for arm, _ in rows) + " |")
        L.append("|---|" + "---|" * len(rows))
        for i in range(len(rows[0][1]["mapq_calibration"])):
            cells = []
            for arm, r in rows:
                c = r["mapq_calibration"][i]
                ex = c["expected_from_mapq"]
                cells.append(f"{c['placed']:,} / {c['wrong']:,} / "
                             f"{'n/a' if c['error_rate'] is None else format(c['error_rate'], '.2e')} / "
                             f"{'n/a' if ex is None else format(ex, '.2e')}")
            L.append(f"| {rows[0][1]['mapq_calibration'][i]['mapq']} | " + " | ".join(cells) + " |")
        L.append("")
        strata_names = list(rows[0][1]["strata"].keys())
        for arm, r in rows:
            for s in r["strata"]:
                if s not in strata_names:
                    strata_names.append(s)
        for s in strata_names:
            L.append(f"## {p}: by {s.replace('_', ' ')}")
            L.append("")
            L.append("| stratum | arm | reads | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong "
                     "| wrong at MAPQ>=30 | MAPQ 0 |")
            L.append("|---|---|---|---|---|---|---|---|---|---|")
            keys = []
            for arm, r in rows:
                for k in r["strata"].get(s, {}):
                    if k not in keys:
                        keys.append(k)
            for k in keys:
                for arm, r in rows:
                    if s not in r["strata"] or k not in r["strata"][s]:
                        continue
                    o = r["strata"][s][k]
                    L.append(f"| {k} | {ARM_LABEL[arm]} | {o['reads']:,} | {pct(o['unmapped_frac'])} | "
                             f"{pct(o['correct50_frac'])} | {pct(o['wrong_frac_of_placed'])} | {pct(o['perfect_frac'], 2)} "
                             f"| {o['perfect_wrong']:,} | {o['wrong_mapq_ge30']:,} ({pct(o['wrong_frac_mapq_ge30'])}) | "
                             f"{pct(o['mapq0_frac'], 2)} |")
            L.append("")
        # cross-arm discordance on the same reads
        R = {}
        for arm, _ in rows:
            f = f"{a.evals}/{arm}.{p}.reads.npz"
            if os.path.exists(f):
                z = np.load(f)
                R[arm] = (z["correct50"], z["b_placed"], z["g_mapq"], z["g_perfect"], z["g_mapped"])
        ct = f"{a.chm13_truth}/{p}.chm13_truth.npz"
        if "chm13" in R and os.path.exists(ct):
            ok = np.load(ct)["ok"].astype(bool)
            L.append(f"## {p}: all arms on the reads that have a one-to-one truth on both GRCh38 and CHM13 "
                     f"({int(ok.sum()):,} of {len(ok):,} reads; GRCh38 arms scored on GRCh38, CHM13 arm on CHM13)")
            L.append("")
            L.append("| arm | unmapped | correct (50 bp) | wrong / placed | perfect | perfect but wrong | "
                     "wrong at MAPQ>=30 | MAPQ 0 |")
            L.append("|---|---|---|---|---|---|---|---|")
            for arm in ARMS:
                if arm not in R:
                    continue
                c, pl, q, pf, gm = (x[ok] for x in R[arm])
                w = pl & ~c
                L.append(f"| {ARM_LABEL[arm]} | {pct((~gm).mean())} | {pct(c.mean())} | {pct(w.sum() / max(pl.sum(), 1))} "
                         f"| {pct(pf.mean(), 2)} | {int((pf & w).sum()):,} | {int((w & (q >= 30)).sum()):,} | "
                         f"{pct(((q == 0) & gm).mean(), 2)} |")
            L.append("")
            R.pop("chm13")
        if len(R) >= 2:
            L.append(f"## {p}: same reads, two arms (correct = within 50 bp; confident wrong = placed wrong with "
                     "MAPQ>=30)")
            L.append("")
            L.append("| A vs B | correct in both | correct only in A | correct only in B | wrong/unplaced in both | "
                     "confident wrong only in A | confident wrong only in B |")
            L.append("|---|---|---|---|---|---|---|")
            arms = [x for x in ARMS if x in R]
            for i in range(len(arms)):
                for j in range(i + 1, len(arms)):
                    A, B = arms[i], arms[j]
                    ca, pa, qa = R[A][:3]
                    cb, pb, qb = R[B][:3]
                    cwa = pa & ~ca & (qa >= 30)
                    cwb = pb & ~cb & (qb >= 30)
                    L.append(f"| {A} vs {B} | {int((ca & cb).sum()):,} | {int((ca & ~cb).sum()):,} | "
                             f"{int((~ca & cb).sum()):,} | {int((~ca & ~cb).sum()):,} | {int((cwa & ~cwb).sum()):,} | "
                             f"{int((cwb & ~cwa).sum()):,} |")
            L.append("")
    with open(f"{a.out}/tables.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    # README: every section except the per-stratum ones (those stay in results/tables.md)
    keep, out_lines = True, []
    for line in L[1:]:
        if line.startswith("## "):
            keep = ": by " not in line
        if keep:
            out_lines.append(line.replace("## ", "### ", 1) if line.startswith("## ") else line)
    seg = f"{a.out}/segment_stats.json"
    if os.path.exists(seg):
        st = json.load(open(seg))
        out_lines += ["### Where the simulation samples (one-to-one segments, GRCh38 non-N bases covered)", "",
                      "| chrom | non-N bases | hap1 covered | hap2 covered | hap1 easy / difficult / extreme covered |",
                      "|---|---|---|---|---|"]
        for c, r in st.items():
            f = r["hap1_smaht_covered_frac"]
            out_lines.append(f"| {c} | {r['non_N_bases']:,} | {pct(r['hap1_covered_frac_of_non_N'], 1)} | "
                             f"{pct(r['hap2_covered_frac_of_non_N'], 1)} | {pct(f['easy'], 1)} / "
                             f"{pct(f['difficult'], 1)} / {pct(f['extreme'], 1)} |")
        out_lines.append("")
    for extra in ("truth_check.txt", "ont_error_profile.txt"):
        f = f"{a.out}/{extra}"
        if os.path.exists(f):
            out_lines += [f"### {extra[:-4].replace('_', ' ')}", "", "```", open(f).read().rstrip(), "```", ""]
    readme = os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md")
    if os.path.exists(readme):
        txt = open(readme).read()
        b, e = "<!-- RESULTS:BEGIN -->", "<!-- RESULTS:END -->"
        if b in txt and e in txt:
            body = ("Generated by summarize.py from the evaluate.py outputs (per-stratum tables: "
                    "`results/tables.md`).\n\n" + "\n".join(out_lines))
            txt = txt[:txt.index(b) + len(b)] + "\n" + body + "\n" + txt[txt.index(e):]
            open(readme, "w").write(txt)
    with open(f"{a.out}/summary.tsv", "w") as fh:
        fh.write("\n".join(tsv) + "\n")
    with open(f"{a.out}/summary.json", "w") as fh:
        json.dump({f"{arm}.{p}": r for (arm, p), r in res.items()}, fh, indent=1)
    print(f"wrote {a.out}/tables.md ({len(res)} arm x platform sets)")


if __name__ == "__main__":
    main()

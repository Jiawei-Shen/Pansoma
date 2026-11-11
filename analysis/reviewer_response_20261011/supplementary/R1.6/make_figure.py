"""Supplementary Fig. R1.6 (R1.6: "full precision-recall curves for every caller").

chr1 precision-recall curves after our population filter, benchmark regions, HG008T and COLO829T x three platforms,
(a) SNVs and (b) INDELs, for the models whose recall R1.6 / Supplementary Table R1.6 report and for both comparators:
- Pansoma SNV: HPRC v2.1 (d46) `--norm dataset` models, three per platform (<S>_<P>_SNV_d46_normbq{,_s1,_s2}, RIS copies in
  ris_7fc7396_20261012; PASS = validation recall 0.9 before the filter). Same runs as E8 latest_models.py.
- Pansoma INDEL: the final models (six tumors, HPRC v2.1), three per platform (release_s{0,1,2}_<P>_INDEL,
  vcf_ppr8_chr1_<set>_d46; PASS = ppr8).
  Curve = report.json pon.curve_roc (rtg vcfeval --all-records ranked by p_somatic; records down to p_somatic 0.01).
- DeepSomatic-TO, ClairS-TO (E8 score report.json, eval_bed): curve = all_pon (every record except RefCall / NoCall,
  whatever its FILTER, after our filter, ranked by QUAL); point = pass_pon (the tool's PASS after the same filter).
Vertical lines (chr1 truth in the benchmark regions, E27 denominators/label_status.json, as in Supplementary Table R1.6):
  dotted   = truth that is not a population allele (the bound for any caller after this filter);
  dash-dot = truth with an exact Pansoma candidate that is not a population allele (Pansoma's bound on this platform).
Writes fig_R1.6.png / .pdf and points.md (the plotted operating points and curve ends, for checking against the table).
Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_figure.py   (login node, < 1 min)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
RIS = RUNS / "ris_7fc7396_20261012"
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
LABEL_STATUS = HERE.parents[1] / "experiments" / "E27_assessment_20261030" / "denominators" / "label_status.json"

TUMORS = ("HG008T", "COLO829T")
PLATFORMS = ("Illumina", "PacBio", "ONT")
KINDS = ("SNV", "INDEL")
# validated categorical slots 1-3 (dataviz palette, light surface); comparators dashed as a second encoding
CALLERS = (("pansoma", "Pansoma", "#2a78d6", "-"), ("deepsomatic", "DeepSomatic-TO", "#eb6834", (0, (5, 2))),
           ("clairs_to", "ClairS-TO", "#1baf7a", (0, (5, 2))))
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6, "figure.facecolor": SURFACE,
                     "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE})


def ds(tumor, p, pansoma=False):
    """Data-set name: E8 uses COLO829T_fiberseq; the RIS SNV runs use HG008 / COLO829 and fiberseq."""
    fs = (tumor, p) == ("COLO829T", "PacBio")
    if pansoma:
        return f"{tumor[:-1]}_{'fiberseq' if fs else p}"
    return f"{tumor}_{'fiberseq' if fs else p}"


def plat(tumor, p):
    return {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if tumor == "COLO829T" else "PacBio HiFi"}[p]


def curve(block):
    return np.array([[c["recall"], c["precision"]] for c in block["curve_roc"]])


def point(block):
    a = block["at_threshold"]
    return a["recall"], a["precision"]


def pansoma_reports(tumor, p, kind):
    if kind == "SNV":
        fs = [RIS / f"{ds(tumor, p, True)}_SNV_d46_normbq{x}" / "vcf_chr1" / "eval_bed" / "report.json" for x in ("", "_s1", "_s2")]
    else:
        fs = [RUNS / f"release_s{s}_{p}_INDEL" / f"vcf_ppr8_chr1_{ds(tumor, p)}_d46" / "eval_bed" / "report.json" for s in (0, 1, 2)]
    return [(f, json.loads(f.read_text())["pon"]) for f in fs]


def comparator_report(caller, tumor, p, kind):
    f = E8 / caller / ds(tumor, p) / kind / "report.json"
    e = json.loads(f.read_text())["eval_bed"]
    return f, e["all_pon"], e["pass_pon"]


def bounds(status, tumor, p, kind):
    """(not a population allele, exact candidate and not a population allele) as shares of chr1 truth in the BED."""
    c = status[f"{ds(tumor, p)}_d46|{kind}"]["chr1"]
    n = c["all"]["n"]
    return (n - c["all"]["_tagged"]) / n, c["untagged"]["exact"] / n


def main():
    status = json.loads(LABEL_STATUS.read_text())
    data, rows = {}, []
    for kind in KINDS:
        for tumor in TUMORS:
            for p in PLATFORMS:
                d = dict(pansoma=pansoma_reports(tumor, p, kind),
                         **{c: comparator_report(c, tumor, p, kind) for c, *_ in CALLERS[1:]},
                         bounds=bounds(status, tumor, p, kind))
                data[(kind, tumor, p)] = d
    xmax = {}
    for kind in KINDS:
        m = 0
        for (k, _, _), d in data.items():
            if k != kind:
                continue
            m = max([m, *d["bounds"]] + [curve(b)[:, 0].max() for _, b in d["pansoma"]] +
                    [curve(d[c][1])[:, 0].max() for c, *_ in CALLERS[1:]])
        xmax[kind] = 1.0 if kind == "SNV" else min(1.0, np.ceil((m + 0.03) * 10) / 10)

    fig, axes = plt.subplots(4, 3, figsize=(7.4, 10.6), sharey=True)
    for bi, kind in enumerate(KINDS):
        for ti, tumor in enumerate(TUMORS):
            i = 2 * bi + ti
            for j, p in enumerate(PLATFORMS):
                ax = axes[i, j]
                d = data[(kind, tumor, p)]
                X = xmax[kind]
                b_any, b_exact = d["bounds"]
                ax.axvline(b_any, color=MUTED, lw=0.9, ls=(0, (1, 1.5)), zorder=1)
                ax.axvline(b_exact, color=MUTED, lw=0.9, ls=(0, (4, 1.5, 1, 1.5)), zorder=1)
                key, label, color, ls = CALLERS[0]
                prs = []
                for f, blk in d["pansoma"]:
                    c = curve(blk)
                    ax.plot(c[:, 0], c[:, 1], color=color, lw=1.1, alpha=0.85, ls=ls, zorder=3, solid_capstyle="round")
                    r, pr = point(blk)
                    prs.append((r, pr, c[:, 0].max()))
                for r, pr, _ in prs:
                    ax.plot(r, pr, marker="o", ms=4.6, color=color, mec=SURFACE, mew=0.8, ls="", zorder=6)
                cells = [f"{min(x[0] for x in prs):.3f}–{max(x[0] for x in prs):.3f}",
                         f"{min(x[1] for x in prs):.3f}–{max(x[1] for x in prs):.3f}",
                         f"{min(x[2] for x in prs):.3f}–{max(x[2] for x in prs):.3f}"]
                for key, label, color, ls in CALLERS[1:]:
                    _, allp, passp = d[key]
                    c = curve(allp)
                    ax.plot(c[:, 0], c[:, 1], color=color, lw=1.3, ls=ls, zorder=4)
                    r, pr = point(passp)
                    ax.plot(r, pr, marker="o", ms=5.2, color=color, mec=SURFACE, mew=0.8, ls="", zorder=7)
                    cells += [f"{r:.3f}", f"{pr:.3f}", f"{c[:, 0].max():.3f}"]
                rows.append(f"| {kind} | {tumor} | {plat(tumor, p)} | " + " | ".join(cells) +
                            f" | {b_exact:.3f} | {b_any:.3f} |")
                ax.set_xlim(0, X)
                ax.set_ylim(0, 1.02)
                ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
                ax.set_title(f"{tumor} · {plat(tumor, p)}", fontsize=8, color=INK, loc="left", pad=3)
                for side in ("top", "right"):
                    ax.spines[side].set_visible(False)
                ax.set_xlabel("Recall (chr1)", fontsize=7.5)
                if j == 0:
                    ax.set_ylabel("Precision (chr1)", fontsize=7.5)
    h = [Line2D([], [], color=c, lw=1.6, ls=ls, label=l) for (_, l, c, ls) in CALLERS]
    h += [Line2D([], [], color=MUTED, marker="o", ls="", ms=5, mec=SURFACE, mew=0.8, label="PASS (after the filter)"),
          Line2D([], [], color=MUTED, lw=0.9, ls=(0, (4, 1.5, 1, 1.5)), label="exact Pansoma candidate, not a population allele"),
          Line2D([], [], color=MUTED, lw=0.9, ls=(0, (1, 1.5)), label="not a population allele (bound for any caller)")]
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, fontsize=6.8, bbox_to_anchor=(0.5, 0.997),
               handlelength=2.6, columnspacing=1.3)
    fig.tight_layout(rect=(0.02, 0, 1, 0.935), h_pad=1.5, w_pad=0.8)
    for bi, kind in enumerate(KINDS):
        top = axes[2 * bi, 0].get_position().y1
        fig.text(0.005, top + 0.016, "ab"[bi], fontsize=11, fontweight="bold", color=INK, ha="left", va="bottom")
        fig.text(0.03, top + 0.017, "SNVs" if kind == "SNV" else "INDELs", fontsize=9, color=INK, ha="left", va="bottom")
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"fig_R1.6.{ext}", dpi=220)
    plt.close(fig)

    L = ["# Supplementary Fig. R1.6: plotted points (chr1, benchmark regions, after our population filter)", "",
         "Generated by `make_figure.py`. Pansoma: range over the three models per platform (SNV `--norm dataset` HPRC v2.1; "
         "INDEL final six-tumor models at ppr8). Comparators: PASS after the same filter; curve end = highest recall of the "
         "QUAL-ranked curve over every record except RefCall / NoCall after the filter.", "",
         "| kind | tumor | platform | Pansoma PASS R | Pansoma PASS P | Pansoma curve end R | DeepSomatic-TO PASS R | P | "
         "curve end R | ClairS-TO PASS R | P | curve end R | exact candidate, not population | not population (any caller) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"] + rows
    (HERE / "points.md").write_text("\n".join(L) + "\n")
    print("fig_R1.6.png/.pdf and points.md written")


if __name__ == "__main__":
    main()

"""Supplementary material of R2.7 (alignment statistics): figure and table from the E11a full-data scans.

E11a full data (`experiments/E11a_alignment_stats/full/summary.json`, `collect_full.py`): every read of each dataset
on four references (GRCh38 linear BAM, GRCh38 graph, HPRC v1.1 graph, HPRC v2.1 graph), counted when its alignment
is placed on chr1-22.
Writes alignment_counts.tsv and fig_R2.7.png next to this script.
usage: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_figures.py   (base anaconda matplotlib is built against NumPy 1.x)
"""
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
FULL = HERE.parents[1] / "experiments" / "E11a_alignment_stats" / "full" / "summary.json"

# Fixed identity per reference (categorical slots 1-4 of the validated default palette, light mode; validated
# 2026-10-19: CVD dE 9.1, normal-vision dE 22.9; aqua / yellow below 3:1 on the surface -> every bar is value-labelled).
ARMS = {"grch38_bam": ("GRCh38 (linear)", "#2a78d6"),
        "linear": ("GRCh38 graph", "#eb6834"),
        "d9": ("HPRC v1.1 graph", "#1baf7a"),
        "d46": ("HPRC v2.1 graph", "#eda100")}
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
ROWS = [("Illumina", [("HG008T_Illumina", "HG008T"), ("COLO829T_Illumina", "COLO829T")]),
        ("PacBio", [("HG008T_PacBio", "HG008T\nHiFi"), ("COLO829T_fiberseq", "COLO829T\nFiber-seq")]),
        ("ONT", [("HG008T_ONT", "HG008T"), ("COLO829T_ONT", "COLO829T")])]
METRICS = [("perfect", "perfectly aligned"), ("sub", "with substitution"), ("indel", "with indel"),
           ("clip", "with soft clipping")]


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def main():
    S = json.loads(FULL.read_text())
    with open(HERE / "alignment_counts.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tumor", "platform", "reference", "mapped reads (chr1-22)"] +
                   [x for _, lab in METRICS for x in (f"reads {lab}", f"% {lab}")])
        for _, sets in ROWS:
            for s, _ in sets:
                for a in ARMS:
                    r = S[f"{s}.{a}"]
                    w.writerow([r["tumor"], r["platform"], ARMS[a][0], r["mapped"]] +
                               [x for k, _ in METRICS for x in (r[k], f"{r[k + '_pct']:.2f}")])

    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": SURFACE})
    fig, axs = plt.subplots(3, 4, figsize=(13, 9))
    width = 0.8 / len(ARMS)
    for i, (plat, sets) in enumerate(ROWS):
        for j, (k, lab) in enumerate(METRICS):
            ax = axs[i][j]
            style(ax)
            top = 0
            for n, a in enumerate(ARMS):
                xs = [g + (n - (len(ARMS) - 1) / 2) * width for g in range(len(sets))]
                ys = [S[f"{s}.{a}"][k + "_pct"] for s, _ in sets]
                top = max(top, *ys)
                ax.bar(xs, ys, width * 0.92, color=ARMS[a][1], label=ARMS[a][0], edgecolor=SURFACE, linewidth=1)
                for x, y in zip(xs, ys):
                    ax.text(x, y, f"{y:.2f}" if y < 5 else f"{y:.1f}", ha="center", va="bottom", fontsize=6.5, color=INK2)
            ax.set_ylim(0, top * 1.15 if top else 1)
            ax.set_xticks(range(len(sets)))
            ax.set_xticklabels([l for _, l in sets], color=INK, fontsize=8)
            if i == 0:
                ax.set_title(f"reads {lab} (%)", loc="left", fontsize=10, color=INK)
            if j == 0:
                ax.set_ylabel(plat, color=INK, fontsize=10)
    h, l = axs[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(HERE / "fig_R2.7.png", dpi=200)
    plt.close(fig)
    print("wrote fig_R2.7.png, alignment_counts.tsv")


if __name__ == "__main__":
    main()

"""Supplementary material of R2.7 (alignment statistics): figures and data tables from E11a and E11b.

E11a (`experiments/E11a_alignment_stats/results/<set>.json`): the same reads on every reference.
  Illumina: window sample W (48 GRCh38 windows, half at segmental duplications; arms GRCh38 BAM, GRCh38-only graph,
  HPRC v1.1 d9, HPRC v2.1 d46). Long reads: genome-wide 10 % name-hash sample G (graph arms).
E11b (`experiments/E11b_simulated_alignment/results/summary.tsv`): reads simulated from the HG008-N assembly, true
  origin known.
Writes alignment_stats.tsv, simulated_placement.tsv and fig_R2.7{a,b,c}.png next to this script.
usage: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_figures.py   (base anaconda matplotlib is built against NumPy 1.x)
"""
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1] / "experiments"
E11A = EXP / "E11a_alignment_stats" / "results"
E11B = EXP / "E11b_simulated_alignment" / "results" / "summary.tsv"

# Fixed identity per reference (categorical slots 1-4 of the validated default palette, light mode; validated
# 2026-10-19: CVD dE 9.1, normal-vision dE 22.9; aqua / yellow below 3:1 on the surface -> every bar is value-labelled).
ARMS = {"grch38_bam": ("GRCh38 BAM (BWA-MEM)", "#2a78d6"),
        "linear": ("GRCh38-only graph", "#eb6834"),
        "d9": ("HPRC v1.1 graph (d9)", "#1baf7a"),
        "d46": ("HPRC v2.1 graph (d46)", "#eda100")}
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SETS = [("HG008T_Illumina", "W", "HG008T Illumina"), ("COLO829T_Illumina", "W", "COLO829T Illumina"),
        ("HG008T_PacBio", "G", "HG008T PacBio HiFi"), ("COLO829T_fiberseq", "G", "COLO829T Fiber-seq"),
        ("HG008T_ONT", "G", "HG008T ONT"), ("COLO829T_ONT", "G", "COLO829T ONT")]
COLS = [("reads", "reads compared"), ("unmapped_pct", "unmapped %"), ("mapq_60_60_pct", "MAPQ 60 %"),
        ("softclip_pct", "soft-clipped reads %"), ("perfect_pct", "full-length edit-free reads %"),
        ("mm_per_kb", "mismatched bases per kb"), ("indel_ev_per_kb", "indel events per kb")]


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def bars(ax, groups, series, fmt="{:.1f}"):
    """groups: x labels; series: [(arm, [value or None per group])]; 2px gaps between bars; value labels."""
    n = len(series)
    width = 0.8 / n
    for i, (arm, vals) in enumerate(series):
        xs = [g + (i - (n - 1) / 2) * width for g in range(len(groups))]
        pts = [(x, v) for x, v in zip(xs, vals) if v is not None]
        if not pts:
            continue
        ax.bar([p[0] for p in pts], [p[1] for p in pts], width * 0.92, color=ARMS[arm][1], label=ARMS[arm][0],
               edgecolor=SURFACE, linewidth=1)
        for x, v in pts:
            ax.text(x, v, fmt.format(v), ha="center", va="bottom", fontsize=7.5, color=INK2)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, color=INK)


def main():
    rows = []
    data = {}
    for name, part, label in SETS:
        d = json.loads((E11A / f"{name}.json").read_text())[part]
        data[name] = d
        for arm in d["arms"]:
            v = d["per_arm"][arm]
            rows.append([label, part, ARMS[arm][0]] + [f"{v[k]:.4g}" if k != "reads" else str(v[k]) for k, _ in COLS])
    with open(HERE / "alignment_stats.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["dataset", "sample (W window / G genome-wide 10 %)", "reference"] + [c[1] for c in COLS])
        w.writerows(rows)

    sim = list(csv.DictReader(open(E11B), delimiter="\t"))
    with open(HERE / "simulated_placement.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["platform", "reference", "reads", "correct %", "wrong of placed %", "wrong at MAPQ>=30 %",
                    "graph-only (no GRCh38 projection) %", "unmapped %"])
        for r in sim:
            n = int(r["reads"])
            w.writerow([r["platform"], ARMS[r["arm"]][0], n, f"{100 * int(r['correct50']) / n:.3f}",
                        f"{100 * float(r['wrong_of_placed']):.3f}", f"{100 * int(r['wrong_mapq_ge30']) / n:.3f}",
                        f"{100 * int(r['graph_only']) / n:.3f}", f"{100 * int(r['unmapped']) / n:.3f}"])

    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": SURFACE})

    # a. Illumina, same reads: full-length edit-free reads
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    style(ax)
    ill = ["HG008T_Illumina", "COLO829T_Illumina"]
    bars(ax, ["HG008T", "COLO829T"],
         [(a, [data[s]["per_arm"][a]["perfect_pct"] if a in data[s]["arms"] else None for s in ill]) for a in ARMS])
    ax.set_ylabel("full-length edit-free reads (%)", color=INK, fontsize=9)
    ax.set_ylim(0, 100)
    ax.set_title("a  Illumina: the same reads on four references", loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "fig_R2.7a.png", dpi=200)
    plt.close(fig)

    # b. long reads: mismatched bases and indel events per kb (two panels, one measure each)
    lr = [("HG008T_PacBio", "HG008T\nPacBio"), ("COLO829T_fiberseq", "COLO829T\nFiber-seq"), ("HG008T_ONT", "HG008T\nONT"),
          ("COLO829T_ONT", "COLO829T\nONT")]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, key, title in ((axs[0], "mm_per_kb", "b  long reads: mismatched bases per kb"),
                           (axs[1], "indel_ev_per_kb", "indel events per kb")):
        style(ax)
        bars(ax, [l for _, l in lr], [(a, [data[s]["per_arm"][a][key] for s, _ in lr]) for a in ("linear", "d9", "d46")],
             fmt="{:.2f}")
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.set_ylim(0, None)
    axs[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(HERE / "fig_R2.7b.png", dpi=200)
    plt.close(fig)

    # c. simulated reads: correct placement (dots, narrow range) and confident misplacements (bars from zero)
    plats = [("illumina", "Illumina"), ("hifi", "HiFi-like"), ("ont", "ONT-like")]
    by = {(r["platform"], r["arm"]): r for r in sim}
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    ax = axs[0]
    style(ax)
    for i, a in enumerate(("linear", "d9", "d46")):
        xs = [g + (i - 1) * 0.22 for g in range(len(plats))]
        ys = [100 * int(by[(p, a)]["correct50"]) / int(by[(p, a)]["reads"]) for p, _ in plats]
        ax.scatter(xs, ys, s=46, color=ARMS[a][1], edgecolor=SURFACE, linewidth=1.5, label=ARMS[a][0], zorder=3)
        for x, y in zip(xs, ys):
            ax.text(x, y + 0.03, f"{y:.2f}", ha="center", va="bottom", fontsize=7.5, color=INK2)
    ax.set_xticks(range(len(plats)))
    ax.set_xticklabels([l for _, l in plats])
    ax.set_title("c  simulated reads: placed correctly (%)", loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax = axs[1]
    style(ax)
    bars(ax, [l for _, l in plats],
         [(a, [1000 * int(by[(p, a)]["wrong_mapq_ge30"]) / int(by[(p, a)]["reads"]) for p, _ in plats])
          for a in ("linear", "d9", "d46")], fmt="{:.2f}")
    ax.set_title("misplaced at MAPQ >= 30 (per 1,000 reads)", loc="left", fontsize=10, color=INK)
    ax.set_ylim(0, None)
    fig.tight_layout()
    fig.savefig(HERE / "fig_R2.7c.png", dpi=200)
    plt.close(fig)
    print("wrote", ", ".join(p.name for p in sorted(HERE.glob("*.png")) + sorted(HERE.glob("*.tsv"))))


if __name__ == "__main__":
    main()

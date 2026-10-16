"""Figure of the E1.1 results (results/summary.json from collect.py) -> results/ch6_counterfactual.png.

Three small multiples, one row per ch6 edit (mode), one dot per finished run, colour = graph (HPRC v1.1 d9,
HPRC v2.1 d46, GRCh38-only), marker = platform:
  (a) per-tensor somatic AP on the chr1 test tensors;
  (b) germline test tensors with p_somatic >= the PASS threshold, relative to the run's base (log scale);
  (c) truth-level F1 at the PASS threshold after the per-tensor PoN.
usage: python plot_modes.py"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
MODES = ("base", "all_H", "train_mean", "train_median", "shuffle", "site0_low", "site0_high", "site10_low", "site10_high")
LABEL = {"base": "base (unchanged)", "all_H": "all cells -> H", "train_mean": "all cells -> training mean",
         "train_median": "all cells -> training median", "shuffle": "shuffled between tensors",
         "site0_low": "site column -> LOW", "site0_high": "site column -> H", "site10_low": "site +-10 -> LOW",
         "site10_high": "site +-10 -> H"}
COLOR = {"d9": "#2a78d6", "d9 joint": "#2a78d6", "d46": "#eb6834", "linear": "#1baf7a"}  # reference palette slots 1-3
MARKER = {"Illumina": "o", "PacBio": "s", "fiberseq": "D", "ONT": "^"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def main():
    s = json.loads((HERE / "results" / "summary.json").read_text())
    runs = s["runs"]
    modes = [m for m in MODES if any(m in r["results"] for r in runs.values())]
    fig, axes = plt.subplots(1, 3, figsize=(13, 0.42 * len(modes) + 1.8), sharey=True)
    y = {m: len(modes) - 1 - k for k, m in enumerate(modes)}
    panels = (("(a) somatic AP per tensor", lambda r, b: r["somatic_ap"], False),
              ("(b) germline tensors >= PASS threshold, x base", lambda r, b: max(r["germline"]["ge_thr"], 1) / max(b["germline"]["ge_thr"], 1), True),
              ("(c) truth F1 at PASS after the per-tensor PoN", lambda r, b: r["truth_after_pon"]["f1"], False))
    for ax, (title, fn, log) in zip(axes, panels):
        for k, (run, r) in enumerate(sorted(runs.items())):
            base = r["results"]["base"]
            plat = next((p for p in MARKER if f"_{p}_" in run), "Illumina")
            for m in modes:
                if m in r["results"]:
                    jitter = (k % 6 - 2.5) * 0.06
                    ax.scatter(fn(r["results"][m], base), y[m] + jitter, s=42, color=COLOR[r["graph"]],
                               marker=MARKER[plat], edgecolors="#fcfcfb", linewidths=1.5, zorder=3)
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        if log:
            ax.set_xscale("log")
            ax.axvline(1, color=MUTED, linewidth=1, zorder=1)
        ax.grid(axis="x", color=GRID, linewidth=1)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
    axes[0].set_yticks([y[m] for m in modes], [LABEL[m] for m in modes])
    graphs = [g for g in ("d9", "d46", "linear") if any(r["graph"].startswith(g) for r in runs.values())]
    plats = [p for p in MARKER if any(f"_{p}_" in run for run in runs)]
    handles = [Line2D([], [], marker="o", linestyle="", color=COLOR[g], label={"d9": "HPRC v1.1 d9", "d46": "HPRC v2.1 d46",
                                                                                  "linear": "GRCh38 only"}[g]) for g in graphs]
    handles += [Line2D([], [], marker=MARKER[p], linestyle="", color=MUTED, label=p) for p in plats]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, fontsize=9)
    fig.suptitle(f"E1.1: SNV models with channel 6 (haplotype path count) edited at inference, chr1 test tensors "
                 f"({len(runs)} runs)", fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(HERE / "results" / "ch6_counterfactual.png", dpi=160)
    print(HERE / "results" / "ch6_counterfactual.png")


if __name__ == "__main__":
    main()

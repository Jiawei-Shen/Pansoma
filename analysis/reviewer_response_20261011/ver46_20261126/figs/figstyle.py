"""Shared style of the bar-figure versions of the ver4.6 letter tables.

The ver4.5 figures in Pansoma_Resubmit/fig_*/ were drawn with ggplot2 (theme of fig_R1.1_SNV_noPoN.R: theme_classic,
7 pt Helvetica, light grey horizontal grid, left axis line, bold facet titles, legend row on top, 600 dpi, 6.5 in
wide). Those scripts and R are not on this cluster, so the ver4.6 figures are drawn with matplotlib to the same
layout. Arial (metric-compatible with Helvetica) is registered from the sigProfilerPlotting package, as
ver4_20261123/pathcount/pathcount_fig.py does.

Run the figure scripts with /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python.
"""
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
V46 = HERE.parent
RESUBMIT = Path("/wanglab/jshen/scratch_jshen/Pansoma_papers/ver1/Pansoma_Resubmit")
sys.path.insert(0, str(V46.parent / "ver4_20261123/build"))

_site = Path(matplotlib.__file__).resolve().parents[1]
for _f in ("Arial.ttf", "Arial Bold.ttf"):
    if (_site / "sigProfilerPlotting" / "fonts" / _f).exists():
        fm.fontManager.addfont(str(_site / "sigProfilerPlotting" / "fonts" / _f))
plt.rcParams.update({"font.family": ["Arial", "Nimbus Sans", "DejaVu Sans"], "pdf.fonttype": 42, "axes.unicode_minus": True,
                     "hatch.linewidth": 0.3})

DPI = 600
PT = 7.0
INK, MUTED, GRID = "#1A1A1A", "#4D4D4D", "#E3E3E3"
C_LINEAR, C_V11, C_V21 = "#F4B3D3", "#C51B7D", "#7B0E4C"
C_DS, C_CLAIRS = "#2CA25F", "#1F4E79"
C_FOUR, C_SNV6 = "#DE6FA9", "#1B7F6B"
C_GERM_IN, C_GERM_OUT, C_NONGERM = "#1F4E79", "#8FB8DE", "#B3B3B3"
LW = 0.3 * 96 / 72 * 0.75           # ggplot linewidth 0.3 mm-ish -> pt (0.3 * .pt / ... ); about 0.64 pt
LW = 0.64
GRID_LW = 0.53
TUMORS = ("HG008T", "COLO829T", "HCC1395T")


def platforms(tumor):
    return ("Illumina", "Fiber-seq" if tumor == "COLO829T" else "PacBio HiFi", "ONT")


def classic(ax, ylim, yticks, left=True, fmt="{:.1f}", grid=True, bottom=True):
    """theme_classic of the R script: grey grid at the y breaks, left and bottom axis lines, no x ticks."""
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_yticklabels([fmt.format(v) for v in yticks] if left else [], fontsize=PT - 0.5, color=INK)
    if grid:
        for v in yticks:
            ax.axhline(v, color=GRID, lw=GRID_LW, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
        ax.spines[s].set_linewidth(LW)
    ax.tick_params(axis="y", length=1.5 if left else 0, width=LW, color=MUTED, pad=1.5)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_facecolor("none")


def top_legend(fig, items, x=0.045, y=0.985, gap=0.012, size=PT, key=7.0):
    """One legend row at the top left: items = [(colour, text)]. Returns the right end (figure fraction)."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W, H = fig.get_size_inches()
    kw, kh = key / 72 / W, key / 72 / H
    for col, text in items:
        fig.patches.append(Rectangle((x, y - kh), kw, kh, transform=fig.transFigure, facecolor=col, edgecolor="none"))
        t = fig.text(x + kw + 3 / 72 / W, y - kh / 2, text, fontsize=size, color=INK, va="center_baseline", ha="left")
        bb = t.get_window_extent(rend)
        x = x + kw + 3 / 72 / W + bb.width / fig.dpi / W + gap + 8 / 72 / W
    return x


def save(fig, stem, folder, name):
    """PNG (600 dpi) + PDF next to the scripts, and a copy in Pansoma_Resubmit/<folder>/<name>_ver4.6.{png,pdf}.
    Only names ending in _ver4.6 are ever written there."""
    outs = []
    for ext in ("png", "pdf"):
        p = HERE / f"{stem}.{ext}"
        fig.savefig(p, dpi=DPI, facecolor="white")
        outs.append(p)
        dst = RESUBMIT / folder / f"{name}_ver4.6.{ext}"
        assert dst.parent.is_dir(), dst.parent
        assert dst.stem.endswith("_ver4.6")
        shutil.copyfile(p, dst)
        outs.append(dst)
    plt.close(fig)
    for o in outs:
        print("written", o)
    return outs

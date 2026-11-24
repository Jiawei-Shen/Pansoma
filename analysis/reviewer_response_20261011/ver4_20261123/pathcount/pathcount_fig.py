#!/usr/bin/env python
"""Response Figure 1.3 (ver4): dumbbell plots + delta heat table, with / without the haplotype-count channel.

Plot step. Reads values.json written by pathcount_data.py and draws fig_R1.3_pathcount.png (and .pdf) in the layout,
fonts and colours of the letter's ver3 image (image3.png, 4251 x 1322 px at 300 dpi; no script of that figure exists,
so its geometry and colours were measured from the image):

  panel a  four dumbbell panels (Precision, Recall, F1 after the population filter; F1 before it). Blue circle = model
           with the channel, seed 0; orange diamond = model without it (seed 0); grey bar joins the two. New in ver4:
           above every row a thin blue rug with one tick per with-channel seed (0, 1, 2) and a line over their range.
           All four panels share one x-range, so equal distances are equal differences in every panel.
  panel b  the differences (with, seed 0 - without) as a heat table, purple-white-green, symmetric, +-VMAX.

The layout is defined in pixels of the 300-dpi image (the constants below are the measured positions of the ver3
image); the image is 4251 px wide and its height follows from the number of rows (1322 px for the six rows of the
ver3 figure, 1104 px for four rows in two groups).

draw() is generic; check_october.py calls it with the six October pairs to redraw the ver3 figure for comparison.

Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python pathcount_fig.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
DPI = 300
PT = DPI / 72.0                      # pixels per point

# colours measured in the ver3 image
INK, BLUE, ORANGE = "#1a1a1a", "#0072b2", "#d55e00"
CONNECT, GUIDE, AXIS = "#9e9e9e", "#e6e6e6", "#4d4d4d"
# PRGn restricted to its middle 70 %: reproduces the cell colours of the ver3 image (see check_october.py)
CMAP = LinearSegmentedColormap.from_list("prgn70", colormaps["PRGn"](np.linspace(0.15, 0.85, 256)))

# geometry of the ver3 image, px
W = 4251
ROW0, PITCH, GROUP_GAP = 495, 109, 77          # first row, row pitch, extra space between groups
SPINE_BELOW, BOTTOM = 60, 145                  # last row -> x axis; x axis -> lower edge
PANELS_RIGHT, PANEL_GAP = 3359, 73             # right end of the fourth panel; gap between panels
LABEL_LEFT, LABEL_PAD = 70, 23                 # left end of the longest row label; label -> first panel
TITLE_Y, LETTER_Y, LEGEND_Y = 394, 117, 231    # baselines of titles / panel letters; centre line of the legend
CELL_X = (3418, 3609, 3801, 3992)              # left edges of the heat cells
CELL_W, CELL_HALF = 177, 46.5
CBAR = (3510, 4077, 194, 231)                  # x0, x1, y0, y1
RUG_DY, RUG_TICK = 33, 11                      # rug line above the row; half tick length

FS_ROW, FS_TICK, FS_TITLE, FS_LETTER, FS_CELL, FS_LEGEND = 12, 12, 13, 16, 10, 12

# The ver3 image is set in Arial. Arial.ttf / Arial Bold.ttf ship with the sigProfilerPlotting package of this
# environment; they are registered here (read only). Without them the metric-compatible Nimbus Sans is used.
import matplotlib.font_manager as fm
for _site in (Path(matplotlib.__file__).resolve().parents[1],):
    for _f in ("Arial.ttf", "Arial Bold.ttf"):
        if (_site / "sigProfilerPlotting" / "fonts" / _f).exists():
            fm.fontManager.addfont(str(_site / "sigProfilerPlotting" / "fonts" / _f))
plt.rcParams.update({"font.family": ["Arial", "Nimbus Sans", "DejaVu Sans"], "pdf.fonttype": 42, "axes.unicode_minus": True})
FONT_USED = Path(fm.findfont("Arial", fallback_to_default=False)).name if any(f.name == "Arial" for f in fm.fontManager.ttflist) else "Nimbus Sans"


def fmt_delta(x):
    return f"{x:+.3f}".replace("-", "−")


def luminance(rgb):
    c = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb[:3]]
    return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]


def draw(rows, out, xlims, xticks, vmax, vmax_label, legend, delta_title, formats=("png", "pdf")):
    """rows: list of dicts {label, group, with: [4 values], without: [4 values], seeds: [[3 values] x 4] or None,
    delta: [4 values]} in the panel order Precision, Recall, F1, F1 (pre-filter).
    xlims / xticks: one (lo, hi) and one tick list per panel. legend: list of (kind, text), kind in circle / rug / diamond."""
    # row positions
    ys, g0, extra = [], rows[0]["group"], 0
    for i, r in enumerate(rows):
        if r["group"] != g0:
            extra, g0 = extra + GROUP_GAP, r["group"]
        ys.append(ROW0 + i * PITCH + extra)
    spine = ys[-1] + SPINE_BELOW
    H = spine + BOTTOM
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.patch.set_facecolor("white")
    ov = fig.add_axes([0, 0, 1, 1])                 # overlay in image pixels (y down)
    ov.set_xlim(0, W)
    ov.set_ylim(H, 0)
    ov.axis("off")
    ov.patch.set_alpha(0)

    # row labels: measure the longest, then place the panels
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    texts = [ov.text(0, y, r["label"], fontsize=FS_ROW, color=INK, ha="right", va="center_baseline")
             for r, y in zip(rows, ys)]
    label_right = LABEL_LEFT + max(t.get_window_extent(rend).width for t in texts)
    for t in texts:
        t.set_x(label_right)
    x0 = label_right + LABEL_PAD
    pw = (PANELS_RIGHT - x0 - 3 * PANEL_GAP) / 4.0
    top = ROW0 - 62

    titles = ("Precision", "Recall", "F1", "F1 (pre-filter)")
    for j in range(4):
        px = x0 + j * (pw + PANEL_GAP)
        ax = fig.add_axes([px / W, 1 - spine / H, pw / W, (spine - top) / H])
        ax.set_xlim(*xlims[j])
        ax.set_ylim(spine, top)
        ax.patch.set_alpha(0)
        for s in ("top", "left", "right"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(AXIS)
        ax.spines["bottom"].set_linewidth(0.85)
        ax.set_yticks([])
        ax.set_xticks(xticks[j])
        ax.set_xticklabels([f"{t:g}" for t in xticks[j]])
        ax.tick_params(axis="x", length=4.3, width=0.72, color=AXIS, labelsize=FS_TICK, labelcolor=INK, pad=2.8)
        lo, hi = xlims[j]
        for r, y in zip(rows, ys):
            ax.plot([lo, hi], [y, y], color=GUIDE, lw=0.72, zorder=1, solid_capstyle="butt")
            a, b = r["with"][j], r["without"][j]
            ax.plot([b, a], [y, y], color=CONNECT, lw=1.68, zorder=2, solid_capstyle="butt")
            if r.get("seeds"):
                s = r["seeds"][j]
                yr = y - RUG_DY
                ax.plot([min(s), max(s)], [yr, yr], color=BLUE, lw=1.0, zorder=3, solid_capstyle="butt")
                for v in s:
                    ax.plot([v, v], [yr - RUG_TICK, yr + RUG_TICK], color=BLUE, lw=1.0, zorder=3, solid_capstyle="butt")
            ax.plot([a], [y], marker="o", ms=9.5, mfc=BLUE, mec="white", mew=0.8, ls="none", zorder=4, clip_on=False)
            ax.plot([b], [y], marker="D", ms=6.65, mfc=ORANGE, mec="white", mew=0.8, ls="none", zorder=5, clip_on=False)
        ov.text(px + pw / 2, TITLE_Y, titles[j], fontsize=FS_TITLE, fontweight="bold", color=INK, ha="center", va="baseline")

    # panel letters
    ov.text(69, LETTER_Y, "a", fontsize=FS_LETTER, fontweight="bold", color=INK, ha="left", va="baseline")
    ov.text(CELL_X[0], LETTER_Y, "b", fontsize=FS_LETTER, fontweight="bold", color=INK, ha="left", va="baseline")

    # legend (one line, centred over the panels; moved left when it would reach panel b)
    MW, MGAP, EGAP = 42, 42, 90
    items = [(k, ov.text(0, LEGEND_Y, t, fontsize=FS_LEGEND, color=INK, ha="left", va="center_baseline")) for k, t in legend]
    widths = [(64 if k == "rug" else MW) + MGAP + t.get_window_extent(rend).width for k, t in items]
    total = sum(widths) + EGAP * (len(items) - 1)
    lx = (x0 + PANELS_RIGHT) / 2 - total / 2
    lx = min(lx, PANELS_RIGHT - 60 - total)
    lx = max(lx, LABEL_LEFT)
    for (k, t), wd in zip(items, widths):
        mw = 64 if k == "rug" else MW
        cx = lx + mw / 2
        if k == "circle":
            ov.plot([cx], [LEGEND_Y], marker="o", ms=11.1, mfc=BLUE, mec="white", mew=0.8, ls="none", clip_on=False)
        elif k == "diamond":
            ov.plot([cx], [LEGEND_Y], marker="D", ms=7.8, mfc=ORANGE, mec="white", mew=0.8, ls="none", clip_on=False)
        else:
            ov.plot([lx, lx + mw], [LEGEND_Y, LEGEND_Y], color=BLUE, lw=1.0, solid_capstyle="butt")
            for v in (lx, lx + 0.62 * mw, lx + mw):
                ov.plot([v, v], [LEGEND_Y - RUG_TICK, LEGEND_Y + RUG_TICK], color=BLUE, lw=1.0, solid_capstyle="butt")
        t.set_x(lx + mw + MGAP)
        lx += wd + EGAP

    # panel b: heat table
    for j, head in enumerate(("P", "R", "F1", "F1 pre")):
        ov.text(CELL_X[j] + CELL_W / 2, TITLE_Y, head, fontsize=FS_TITLE, fontweight="bold", color=INK, ha="center", va="baseline")
        for r, y in zip(rows, ys):
            d = r["delta"][j]
            col = CMAP(0.5 + 0.5 * max(-1.0, min(1.0, d / vmax)))
            ov.add_patch(Rectangle((CELL_X[j], y - CELL_HALF), CELL_W, 2 * CELL_HALF, facecolor=col, edgecolor="none"))
            ov.text(CELL_X[j] + CELL_W / 2, y + 1, fmt_delta(d), fontsize=FS_CELL, ha="center", va="center_baseline",
                    color="white" if luminance(col) < 0.42 else INK)
    grad = np.linspace(0, 1, 512)[None, :]
    ov.imshow(grad, cmap=CMAP, extent=(CBAR[0], CBAR[1], CBAR[3], CBAR[2]), aspect="auto", interpolation="bilinear")
    ov.set_xlim(0, W)
    ov.set_ylim(H, 0)
    for x, s in ((CBAR[0] + 28, "−" + vmax_label), ((CBAR[0] + CBAR[1]) / 2, "0"), (CBAR[1] - 28, "+" + vmax_label)):
        ov.text(x, 292, s, fontsize=FS_TICK, color=INK, ha="center", va="baseline")
    for k, line in enumerate(delta_title):
        ov.text(CBAR[0], 109 + 55 * k, line, fontsize=FS_LEGEND, color=INK, ha="left", va="baseline")

    outs = []
    for ext in formats:
        p = Path(str(out) + "." + ext)
        fig.savefig(p, dpi=DPI, facecolor="white")
        outs.append(p)
    plt.close(fig)
    return outs, (W, H), {"label_right": label_right, "panel_x0": x0, "panel_width": pw, "rows_y": ys}


# ver4 settings -------------------------------------------------------------------------------------------------------
XLIM = (0.575, 0.925)            # one range for the four panels
XTICKS = [0.6, 0.7, 0.8, 0.9]
VMAX = 0.15                      # panel b colour scale, symmetric (ver3: 0.5)


def main():
    v = json.load(open(HERE / "values.json"))
    metrics = v["metrics"]
    rows = []
    allv = []
    for r in v["rows"]:
        m = r["models"]
        w = [m["with_s0"][st][k] for k, st, _ in metrics]
        wo = [m["without_s0"][st][k] for k, st, _ in metrics]
        seeds = [[m[f"with_s{s}"][st][k] for s in range(3)] for k, st, _ in metrics]
        delta = [r["delta"][short] for _, _, short in metrics]
        assert all(abs((a - b) - d) < 1e-12 for a, b, d in zip(w, wo, delta))
        rows.append({"label": r["label"], "group": r["platform"], "with": w, "without": wo, "seeds": seeds, "delta": delta})
        allv += w + wo + [x for s in seeds for x in s]
    assert XLIM[0] + 0.02 < min(allv) and max(allv) < XLIM[1] - 0.02, (min(allv), max(allv))
    assert max(abs(d) for r in rows for d in r["delta"]) <= VMAX
    outs, size, geo = draw(
        rows, HERE / "fig_R1.3_pathcount", [XLIM] * 4, [XTICKS] * 4, VMAX, f"{VMAX:g}",
        legend=[("circle", "With haplotype-count channel (seed 0)"),
                ("rug", "Seeds 0, 1, 2 with the channel (range)"),
                ("diamond", "Without haplotype-count channel")],
        delta_title=("Δ (with, seed 0 − without", "haplotype-count channel)"))
    print("written", [str(o) for o in outs], "size px", size, geo)
    json.dump({"png": str(outs[0]), "width_px": size[0], "height_px": size[1], "dpi": DPI, "xlim_all_panels": XLIM,
               "colour_scale": [-VMAX, VMAX], "geometry_px": geo}, open(HERE / "fig_R1.3_pathcount.meta.json", "w"), indent=1)


if __name__ == "__main__":
    main()

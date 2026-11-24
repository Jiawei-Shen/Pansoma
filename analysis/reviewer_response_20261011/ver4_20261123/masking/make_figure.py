"""Response Figure 1.2 (ver4): masking graph context outside the candidate node at inference.

Redraw of the ver3 letter figure (v3/media/image2.png, 4251 x 1322 px, six rows) for the nine comparisons of
results.json (collect.py): the six-tumor current SNV models on chr1 of HG008T, COLO829T and HCC1395T.
  a  dumbbells for Precision / Recall / F1 (chr1 SNV, benchmark regions, PASS calls after the population filter):
     blue circle = full tensor, orange diamond = own node only (every cell outside the candidate's node visit set to 0).
  b  heat table of the differences full tensor - own node; diverging purple-green scale, symmetric around 0, its limit
     the largest |difference| rounded up to the next 0.05 (printed when the script runs).
Layout in pixels of the 600-dpi image, as the ver3 figure: width 4251 px (7.085 in), row pitch 117 px, 83 px more
between the tumor groups; three more rows and one more group gap make the height 1757 px (2.928 in).
Colours as the ver3 figure: Okabe-Ito blue #0072B2 / vermillion #D55E00 (two series, shape is the second encoding),
PRGn (its 10 % - 90 % range) for the differences. Font: Nimbus Sans (Helvetica metrics, the ver3 figure's sans serif).
Writes fig_1.2_masking.png (+ .pdf). usage: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_figure.py [--out PNG] [--json FILE]
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter, MultipleLocator  # noqa: E402

HERE = Path(__file__).resolve().parent
DPI = 600
W = 4251
ROW0, PITCH, GROUP_GAP, BELOW = 445, 117, 83, 210  # first row, row pitch, extra gap between tumors, below the last row
LABEL_RIGHT = 580
PANEL_L, PANEL_R, PANEL_GAP = 605, 3402, 72
CELL_X, CELL_W, CELL_H = (3462, 3703, 3944), 226, 99
BLUE, ORANGE, INK, GRID, LINK = "#0072B2", "#D55E00", "#1a1a1a", "#e4e4e4", "#9a9a9a"
PRGN = LinearSegmentedColormap.from_list("PRGn_10_90", plt.get_cmap("PRGn")(np.linspace(0.1, 0.9, 256)))  # the ver3 figure's ends
FS = 6.0

plt.rcParams.update({"font.family": ["Nimbus Sans", "DejaVu Sans"], "font.size": FS, "axes.unicode_minus": True,
                     "pdf.fonttype": 42, "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK})


def signed(x, nd=3):
    s = f"{x:+.{nd}f}".replace("-", "−")
    return s


def nice_axis(values):
    """(lo, hi, tick step) for one metric: the data with a margin, ticks at multiples of 0.05 / 0.1 / 0.2."""
    lo, hi = min(values), max(values)
    span = max(hi - lo, 0.05)
    lo, hi = lo - 0.09 * span, hi + 0.09 * span
    for step in (0.02, 0.05, 0.1, 0.2):
        if (np.floor(hi / step + 1e-9) - np.ceil(lo / step - 1e-9) + 1) <= 5:
            break
    return lo, hi, step


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=str(HERE / "results.json"))
    ap.add_argument("--out", default=str(HERE / "fig_1.2_masking.png"))
    args = ap.parse_args()
    rows = json.load(open(args.json))["rows"]
    pending = [r["set"] for r in rows if not r.get("own_node")]
    if pending:
        raise SystemExit(f"own-node arm pending for {pending}")
    ys, y, last = [], ROW0, None
    for r in rows:
        if last is not None:
            y += PITCH + (GROUP_GAP if r["tumor"] != last else 0)
        ys.append(y)
        last = r["tumor"]
    H = ys[-1] + BELOW
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.patch.set_facecolor("white")
    fx = lambda px: px / W  # noqa: E731
    fy = lambda px: 1 - px / H  # noqa: E731
    top, axis_y = ROW0 - 62, ys[-1] + 65

    # row labels
    for r, yy in zip(rows, ys):
        fig.text(fx(LABEL_RIGHT), fy(yy), f"{r['tumor']} {r['platform']}", ha="right", va="center_baseline", fontsize=FS)

    # panel a
    fig.text(fx(69), fy(96), "a", fontsize=8.5, fontweight="bold", ha="left", va="center_baseline")
    pw = (PANEL_R - PANEL_L - 2 * PANEL_GAP) / 3
    for k, (key, name) in enumerate((("P", "Precision"), ("R", "Recall"), ("F1", "F1"))):
        x0 = PANEL_L + k * (pw + PANEL_GAP)
        ax = fig.add_axes([fx(x0), fy(axis_y), pw / W, (axis_y - top) / H])
        full = [r["full"][key] for r in rows]
        own = [r["own_node"][key] for r in rows]
        lo, hi, step = nice_axis(full + own)
        ax.set_xlim(lo, hi)
        ax.set_ylim(axis_y, top)  # pixel rows, top down
        for yy, a, b in zip(ys, full, own):
            ax.plot([lo, hi], [yy, yy], color=GRID, lw=0.5, zorder=1, solid_capstyle="butt")
            ax.plot([a, b], [yy, yy], color=LINK, lw=1.15, zorder=2, solid_capstyle="butt")
            ax.plot([a], [yy], marker="o", ms=4.7, mfc=BLUE, mec="white", mew=0.55, zorder=3, ls="none", clip_on=False)
            ax.plot([b], [yy], marker="D", ms=3.35, mfc=ORANGE, mec="white", mew=0.55, zorder=4, ls="none", clip_on=False)
        for s in ("left", "right", "top"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.6)
        ax.spines["bottom"].set_color("#333333")
        ax.set_yticks([])
        ax.xaxis.set_major_locator(MultipleLocator(step))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _, s=step: f"{v:.2f}" if round(s * 100) % 10 else f"{v:.1f}"))
        ax.tick_params(axis="x", length=2.2, width=0.6, pad=2.2, labelsize=FS, color="#333333")
        ax.patch.set_alpha(0)
        fig.text(fx(x0 + pw / 2), fy(322), name, fontsize=6.5, fontweight="bold", ha="center", va="center_baseline")
    handles = [Line2D([], [], marker="o", ms=4.7, mfc=BLUE, mec="white", mew=0.55, ls="none", label="Full tensor"),
               Line2D([], [], marker="D", ms=3.35, mfc=ORANGE, mec="white", mew=0.55, ls="none", label="Own node only")]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(fx((PANEL_L + PANEL_R) / 2), fy(177)), ncol=2, frameon=False,
               fontsize=FS, handletextpad=0.5, columnspacing=1.6, borderpad=0)

    # panel b
    deltas = [r["delta"][k] for r in rows for k in ("P", "R", "F1")]
    vmax = max(0.05, np.ceil(max(abs(d) for d in deltas) / 0.05 - 1e-9) * 0.05)
    print(f"panel b: differences {min(deltas):+.4f} to {max(deltas):+.4f}; colour scale -{vmax:.2f} to +{vmax:.2f}")
    fig.text(fx(CELL_X[0] + 12), fy(96), "b", fontsize=8.5, fontweight="bold", ha="left", va="center_baseline")
    fig.text(fx(CELL_X[0] + 83), fy(96), "Δ (full tensor − own node)", fontsize=FS, ha="left", va="center_baseline")
    bar_l, bar_r = CELL_X[0] + 83, CELL_X[0] + 69 + 566
    cax = fig.add_axes([fx(bar_l), fy(177), (bar_r - bar_l) / W, 37 / H])
    cax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=PRGN)
    cax.set_axis_off()
    for px, lab in ((bar_l, signed(-vmax, 2)), ((bar_l + bar_r) / 2, "0"), (bar_r, signed(vmax, 2))):
        fig.text(fx(px), fy(221), lab, fontsize=FS, ha="center", va="center_baseline")
    for cx, head in zip(CELL_X, ("P", "R", "F1")):
        fig.text(fx(cx + CELL_W / 2), fy(322), head, fontsize=6.5, fontweight="bold", ha="center", va="center_baseline")
    bax = fig.add_axes([0, 0, 1, 1])
    bax.set_xlim(0, W); bax.set_ylim(H, 0); bax.set_axis_off(); bax.patch.set_alpha(0)
    for r, yy in zip(rows, ys):
        for cx, key in zip(CELL_X, ("P", "R", "F1")):
            d = r["delta"][key]
            rgba = PRGN(0.5 + 0.5 * max(-1.0, min(1.0, d / vmax)))
            bax.add_patch(plt.Rectangle((cx, yy - CELL_H / 2), CELL_W, CELL_H, fc=rgba, ec="none"))
            lum = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
            bax.text(cx + CELL_W / 2, yy, signed(d), ha="center", va="center_baseline", fontsize=5.3,
                     color="white" if lum < 0.62 else INK)
    out = Path(args.out)
    fig.savefig(out, dpi=DPI, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    from PIL import Image
    print(f"wrote {out} {Image.open(out).size} px at {DPI} dpi; {len(rows)} rows; canvas {W} x {H}")


if __name__ == "__main__":
    main()

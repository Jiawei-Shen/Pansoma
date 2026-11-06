"""Figures and table for R3.9 point 4 from $DATA/curves.json (collect.py).

fig_pr_{SNV,INDEL}_{pon,raw}.png/.pdf: precision-recall curves on chr1 (benchmark BED), tumor x platform grid; filled
  circle = PASS (the model's operating threshold; comparators: the tool's PASS), open diamond = p_somatic 0.5 (Pansoma);
  thick line = seed 0, thin = seeds 1-2. 'pon' = after our population filter (main protocol), 'raw' = before it.
fig_threshold_{SNV,INDEL}_{pon,raw}.png/.pdf: F1 against the p_somatic threshold (Pansoma only; the comparators' QUAL is
  on another scale).
F1_threshold_curve_v2/fig_threshold_{SNV,INDEL}_{pon,raw}.png/.pdf: the same for the four-line and six-tumor models only
  (seed 0), with vertical lines at p_somatic 0.5 and at each model's PASS threshold.
summary.md: per panel and model, three-seed mean (or the single run) of F1 at PASS, F1 at p_somatic 0.5, best F1 on
  the curve (test-tuned upper bound), and the PASS threshold range.
Run: python plot.py   (login node, a minute)
"""
import json
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/E27_pr_curves")
TUMORS = ("HG008T", "COLO829T", "HCC1395T")
PLATFORMS = ("Illumina", "PacBio", "ONT")
# series identity: (family, ref) -> label, colour (validated categorical order, light surface), dashed
SERIES = [
    (("single", "linear"), "GRCh38-only, single-tumor", "#2a78d6", False),
    (("single", "d9"), "HPRC v1.1, single-tumor", "#eb6834", False),
    (("single", "d46"), "HPRC v2.1, single-tumor", "#1baf7a", False),
    (("four-line", "d46"), "HPRC v2.1, four-line", "#eda100", False),
    (("six-tumor", "d46"), "HPRC v2.1, six-tumor", "#e87ba4", False),
    (("deepsomatic", "GRCh38"), "DeepSomatic-TO", "#008300", True),
    (("clairs_to", "GRCh38"), "ClairS-TO", "#6250d6", True),
]
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6, "figure.facecolor": SURFACE,
                     "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE})


def plat(tumor, p):
    return {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if tumor == "COLO829T" else "PacBio HiFi"}[p]


def at_p(curve, t):
    """The curve point at threshold t: the last point with p_somatic >= t (None when no record reaches t)."""
    sel = [c for c in curve if c[0] is not None and c[0] >= t]
    return sel[-1] if sel else None


def f1_at(curve, t):
    """F1 when calling at threshold t: 0 when no record reaches t (no calls)."""
    h = at_p(curve, t)
    return 0.0 if h is None else h[3]


def no_curve(panels, kind):
    """Pansoma series of `kind` without a curve, as 'family ref' labels."""
    lab = {k: l for (k, l, _, _) in SERIES}
    miss = {}
    for k, v in panels.items():
        if k.endswith(kind):
            for s in v:
                if s["pon"]["curve"] is None:
                    miss.setdefault(lab[(s["family"], s["ref"])], set()).add(s["seed"])
    return [l + ("" if 0 in seeds else " (seeds " + ", ".join(map(str, sorted(seeds))) + ")") for l, seeds in sorted(miss.items())]


def grid(kind, stage, xmax):
    fig, axes = plt.subplots(3, 3, figsize=(7.4, 7.9), sharex=True, sharey=True)
    return fig, axes


def iso_f1(ax, xmax):
    """Iso-F1 contours (dotted), labelled at the right edge."""
    r = np.linspace(0.005, xmax, 300)
    for f in (0.2, 0.4, 0.6, 0.8):
        with np.errstate(divide="ignore", invalid="ignore"):
            p = f * r / (2 * r - f)
        ok = (r > f / 2) & (p > 0) & (p <= 1.0)
        ax.plot(r[ok], p[ok], color=GRID, lw=0.6, ls=(0, (2, 2)), zorder=0)
        pe = f * xmax / (2 * xmax - f) if 2 * xmax != f else 2.0
        if 0 < pe <= 1:
            ax.text(xmax, pe + 0.012, f"F1 {f}", color="#9a9994", fontsize=5.5, ha="right", va="bottom")


def legend(fig, kind, extra):
    h = [Line2D([], [], color=c, lw=1.6, ls="--" if d else "-", label=l) for (_, l, c, d) in SERIES]
    h += [Line2D([], [], color=MUTED, lw=1.5, label="seed 0"), Line2D([], [], color=MUTED, lw=0.7, alpha=0.6, label="seeds 1–2")]
    h += [Line2D([], [], color=MUTED, marker="o", ls="", ms=5.5, mec=SURFACE, mew=0.8, label="PASS (operating threshold)")]
    if extra:
        h += [Line2D([], [], color=MUTED, marker="D", ls="", ms=4.5, mfc=SURFACE, mew=1.1, label="p_somatic = 0.5")]
    fig.legend(handles=h, loc="upper center", ncol=4, frameon=False, fontsize=6.8, bbox_to_anchor=(0.5, 0.995),
               handlelength=2.4, columnspacing=1.2)


def pr_figure(panels, kind, stage, out):
    def row_max(tumor):
        if kind == "SNV":
            return 1.0
        m = max(max([c[2] for c in (s[stage]["curve"] or [])] + [s[stage]["op"]["R"] if s[stage]["op"] else 0])
                for p in PLATFORMS for s in panels[f"{tumor}|{p}|{kind}"])
        return min(1.0, np.ceil((m + 0.02) * 10) / 10)
    fig, axes = plt.subplots(3, 3, figsize=(7.4, 8.0), sharex="row" if kind == "INDEL" else True, sharey=True)
    for i, tumor in enumerate(TUMORS):
        xmax = row_max(tumor)
        for j, p in enumerate(PLATFORMS):
            ax = axes[i, j]
            iso_f1(ax, xmax)
            series = panels[f"{tumor}|{p}|{kind}"]
            for (key, label, color, dashed) in SERIES:
                ss = [s for s in series if (s["family"], s["ref"]) == key]
                for s in sorted(ss, key=lambda s: -s["seed"]):
                    first = s["seed"] == 0
                    st = s[stage]
                    if st["curve"]:
                        c = np.array([[x[2], x[1]] for x in st["curve"]])
                        ax.plot(c[:, 0], c[:, 1], color=color, lw=1.5 if first else 0.7, alpha=1 if first else 0.5,
                                ls="--" if dashed else "-", zorder=3 if first else 2, solid_capstyle="round")
                        if not dashed:
                            h = at_p(st["curve"], 0.5)
                            if h:
                                ax.plot(h[2], h[1], marker="D", ms=4.2 if first else 3.2, mfc=SURFACE, mec=color,
                                        mew=1.1, ls="", zorder=5, alpha=1 if first else 0.7)
                    if st["op"]:
                        ax.plot(st["op"]["R"], st["op"]["P"], marker="o", ms=5.5 if first else 4.0, color=color,
                                mec=SURFACE, mew=0.8, ls="", zorder=6, alpha=1 if first else 0.75)
            ax.set_xlim(0, xmax)
            ax.set_ylim(0, 1.02)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
            ax.set_title(f"{tumor} · {plat(tumor, p)}", fontsize=8, color=INK, loc="left", pad=3)
            ax.grid(False)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            if i == 2 or kind == "INDEL":
                ax.set_xlabel("Recall (chr1)")
            if j == 0:
                ax.set_ylabel("Precision (chr1)")
    legend(fig, kind, True)
    bottom = 0.0
    miss = no_curve(panels, kind)
    if miss:
        fig.text(0.01, 0.005, "PASS points only (no chr1 curve): " + "; ".join(miss) + ".", fontsize=6, color=MUTED,
                 ha="left", va="bottom")
        bottom = 0.02
    fig.tight_layout(rect=(0, bottom, 1, 0.905), h_pad=1.2, w_pad=0.8)
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=220)
    plt.close(fig)


def threshold_figure(panels, kind, stage, out):
    fig, axes = plt.subplots(3, 3, figsize=(7.4, 8.0), sharex=True, sharey=True)
    ymax = 0
    for i, tumor in enumerate(TUMORS):
        for j, p in enumerate(PLATFORMS):
            ax = axes[i, j]
            ax.axvline(0.5, color=GRID, lw=0.8, zorder=0)
            series = panels[f"{tumor}|{p}|{kind}"]
            for (key, label, color, dashed) in SERIES:
                if dashed:
                    continue
                for s in sorted([s for s in series if (s["family"], s["ref"]) == key], key=lambda s: -s["seed"]):
                    first = s["seed"] == 0
                    st = s[stage]
                    if not st["curve"]:
                        continue
                    c = np.array([[x[0], x[3]] for x in st["curve"] if x[0] is not None])
                    ymax = max(ymax, c[:, 1].max())
                    ax.plot(c[:, 0], c[:, 1], color=color, lw=1.5 if first else 0.7, alpha=1 if first else 0.5,
                            zorder=3 if first else 2)
                    t = s.get("threshold")
                    if st["op"] and t is not None and "op_source" not in s:
                        ax.plot(t, st["op"]["F1"], marker="o", ms=5.5 if first else 4.0, color=color, mec=SURFACE,
                                mew=0.8, ls="", zorder=6, alpha=1 if first else 0.75)
            ax.set_xlim(0, 1)
            ax.set_title(f"{tumor} · {plat(tumor, p)}", fontsize=8, color=INK, loc="left", pad=3)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            if i == 2:
                ax.set_xlabel("p_somatic threshold")
            if j == 0:
                ax.set_ylabel("F1 (chr1)")
    for ax in axes.flat:
        ax.set_ylim(0, min(1.0, ymax + 0.05))
    h = [Line2D([], [], color=c, lw=1.6, label=l) for (_, l, c, d) in SERIES if not d]
    h += [Line2D([], [], color=MUTED, marker="o", ls="", ms=5.5, mec=SURFACE, mew=0.8, label="PASS threshold"),
          Line2D([], [], color=GRID, lw=0.8, label="p_somatic = 0.5")]
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 0.995),
               handlelength=2.4, columnspacing=1.4)
    fig.tight_layout(rect=(0, 0, 1, 0.925), h_pad=1.2, w_pad=0.8)
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=220)
    plt.close(fig)


def threshold_figure_v2(panels, kind, stage, out):
    """F1 against the p_somatic threshold, four-line and six-tumor models only (seed 0), with vertical lines at
    p_somatic 0.5 and at each model's PASS threshold (one model per platform, so the same line in every tumor)."""
    keep = [x for x in SERIES if x[0][0] in ("four-line", "six-tumor")]
    fig, axes = plt.subplots(3, 3, figsize=(7.4, 8.0), sharex=True, sharey=True)
    ymax = 0
    for i, tumor in enumerate(TUMORS):
        for j, p in enumerate(PLATFORMS):
            ax = axes[i, j]
            ax.axvline(0.5, color=GRID, lw=1.0, zorder=0)
            for (key, label, color, _) in keep:
                s = next(s for s in panels[f"{tumor}|{p}|{kind}"] if (s["family"], s["ref"]) == key and s["seed"] == 0)
                st, t = s[stage], s["threshold"]
                c = np.array([[x[0], x[3]] for x in st["curve"] if x[0] is not None])
                ymax = max(ymax, c[:, 1].max())
                ax.axvline(t, color=color, lw=1.0, ls=(0, (4, 2)), zorder=1)
                ax.plot(c[:, 0], c[:, 1], color=color, lw=1.5, zorder=3, solid_capstyle="round")
                ax.plot(t, st["op"]["F1"], marker="o", ms=5.5, color=color, mec=SURFACE, mew=0.8, ls="", zorder=6)
            ax.set_xlim(0, 1)
            ax.set_title(f"{tumor} · {plat(tumor, p)}", fontsize=8, color=INK, loc="left", pad=3)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            if i == 2:
                ax.set_xlabel("p_somatic threshold")
            if j == 0:
                ax.set_ylabel("F1 (chr1)")
    for ax in axes.flat:
        ax.set_ylim(0, min(1.0, ymax + 0.05))
    h = [Line2D([], [], color=c, lw=1.6, label=f"{l} (seed 0)") for (_, l, c, _) in keep]
    h += [Line2D([], [], color=c, lw=1.0, ls=(0, (4, 2)), label=f"{l.split(', ')[1]} PASS threshold") for (_, l, c, _) in keep]
    h += [Line2D([], [], color=GRID, lw=1.0, label="p_somatic = 0.5"),
          Line2D([], [], color=MUTED, marker="o", ls="", ms=5.5, mec=SURFACE, mew=0.8, label="F1 at PASS")]
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 0.995),
               handlelength=2.4, columnspacing=1.4)
    fig.tight_layout(rect=(0, 0, 1, 0.925), h_pad=1.2, w_pad=0.8)
    out.parent.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=220)
    plt.close(fig)


def fmt(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return "–"
    return f"{xs[0]:.3f}" if len(xs) == 1 else f"{statistics.mean(xs):.3f} ({min(xs):.3f}–{max(xs):.3f})"


def summary(panels):
    L = ["# R3.9 operating point: chr1 F1 at PASS, at p_somatic 0.5 and at the curve's best point", "",
         "Benchmark regions, chr1. Mean over seeds (range); single-run families (SNV four-line / six-tumor, comparators) one value.",
         "Best = highest F1 on the chr1 curve (threshold chosen on the test data: an upper bound). Comparators: tool PASS and",
         "the best point of their QUAL-ranked curve over all records except RefCall / NoCall. F1 at p 0.5 is 0 when no record",
         "reaches 0.5 (no calls). Curves start at p_somatic 0.01 (graph_vcf writes no record below it).",
         "Before our filter the comparators' curve also holds the records they filter out themselves (non-PASS), most of which",
         "our population filter removes: its best point lies below their PASS, so before the filter compare their PASS only.", ""]
    for kind in ("SNV", "INDEL"):
        miss = no_curve(panels, kind)
        if miss:
            L += [f"{kind}: PASS and best only (no chr1 curve): " + "; ".join(miss) + ".", ""]
    for kind in ("SNV", "INDEL"):
        for stage, sname in (("pon", "after our population filter"), ("raw", "before it")):
            L += [f"## {kind}, {sname}", "",
                  "| tumor | platform | model | seeds | PASS threshold | F1 at PASS | F1 at p 0.5 | best F1 | best − PASS |",
                  "|---|---|---|---:|---|---|---|---|---|"]
            for tumor in TUMORS:
                for p in PLATFORMS:
                    for (key, label, color, dashed) in SERIES:
                        ss = [s for s in panels[f"{tumor}|{p}|{kind}"] if (s["family"], s["ref"]) == key]
                        if not ss:
                            continue
                        pas = [s[stage]["op"]["F1"] if s[stage]["op"] else None for s in ss]
                        half = [f1_at(s[stage]["curve"], 0.5) if s[stage]["curve"] and not dashed else None for s in ss]
                        best = [s[stage]["best"]["F1"] if s[stage]["best"] else None for s in ss]
                        thr = [s.get("threshold") for s in ss if s.get("threshold") is not None and "op_source" not in s]
                        gap = [b - a for a, b in zip(pas, best) if a is not None and b is not None]
                        L.append(f"| {tumor} | {plat(tumor, p)} | {label} | {len(ss)} | "
                                 f"{'–' if not thr else (f'{thr[0]:.3f}' if len(thr) == 1 else f'{min(thr):.3f}–{max(thr):.3f}')} | "
                                 f"{fmt(pas)} | {fmt(half)} | {fmt(best)} | {fmt(gap)} |")
            L.append("")
    (HERE / "summary.md").write_text("\n".join(L))


def main():
    d = json.loads((DATA / "curves.json").read_text())
    panels = d["panels"]
    for kind in ("SNV", "INDEL"):
        for stage in ("pon", "raw"):
            pr_figure(panels, kind, stage, HERE / f"fig_pr_{kind}_{stage}")
            threshold_figure(panels, kind, stage, HERE / f"fig_threshold_{kind}_{stage}")
            threshold_figure_v2(panels, kind, stage, HERE / "F1_threshold_curve_v2" / f"fig_threshold_{kind}_{stage}")
    summary(panels)
    print("figures + summary.md written")


if __name__ == "__main__":
    main()

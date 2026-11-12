"""Supplementary Fig. R1.7 (R1.7: "provide evidence that the data volume supports the chosen capacity").

Training curves of the six-tumor models on the HPRC v2.1 graph (d46 filter), one run (seed 0) per platform, trained on
chr2-22 of HCC1395T, HCC1954T, H1437T, H2009T, HG008T and COLO829T (HCC1395T without chr20; chr1 of every tumor held
out), the models of the letter's six-tumor rows:
  (a) SNVs, 50.4 M-parameter network, 20 epochs:
      Illumina / PacBio HiFi  snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe (current SNV pipeline: site features,
                              AF 0.05 sets, non with AF < 0.1 sampled 25% per epoch at weight 1, selected by post-filter
                              validation F1)
      ONT                     multicancer6_ONT_SNV (the letter's six-tumor ONT SNV model, "normfix" predictions;
                              selected by validation somatic F1 before the filter). The ONT model on the current SNV
                              pipeline (snv_newpipe_20261109/multicancer6_ONT_SNV_newpipe, job 477091) replaces it via
                              --ont-snv newpipe once trained.
  (b) INDELs, 8.7 M-parameter network, 12 epochs: release_s0_<P>_INDEL (the final INDEL model, selected by post-filter
      validation F1).
Each panel shows epochs 1 to the selected epoch + 2 (user 2026-11-12); data.md keeps every epoch.
Per epoch (metrics.jsonl): training loss = mean class-weighted cross-entropy over that epoch's batches (stochastic depth
on, weights changing during the epoch); validation loss = the same class-weighted cross-entropy on the validation set
(5% of the 20,000-node blocks of the training chromosomes, disjoint from training) after the epoch; validation somatic
F1 = the metric that selects the checkpoint (args.json --best-by: postpon = after our population filter at the
validation best-F1 threshold; f1 = before the filter). Selected epoch = train.log "done; best epoch" (asserted equal to
the maximum of the plotted metric).
Writes fig_R1.7.png / .pdf and data.md (per model: parameters, training / validation examples per class, selected
epoch, losses; per-epoch values).
Size: 6.5 in wide (Word text width, insert at 100%), 300 dpi.
Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_figure.py [--ont-snv newpipe]   (login node, seconds)
"""
import argparse
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
PLATFORMS = ("Illumina", "PacBio", "ONT")
HEADS = ("Illumina", "PacBio HiFi (COLO829T: Fiber-seq)", "ONT")
KINDS = ("SNV", "INDEL")
AFTER = 2  # epochs drawn after the selected one (user 2026-11-12: "就画选中的epoch加上后续的两个epoch")
# same validated categorical slots as Supplementary Fig. R1.6 (light surface); validation also carries markers
TRAIN_C, VAL_C = "#2a78d6", "#eb6834"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
                     "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6, "figure.facecolor": SURFACE,
                     "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE})


def run_dir(kind, p, ont_snv):
    if kind == "INDEL":
        return RUNS / f"release_s0_{p}_INDEL"
    if p == "ONT" and ont_snv == "letter":
        return RUNS / "multicancer6_ONT_SNV"
    return RUNS / "snv_newpipe_20261109" / f"multicancer6_{p}_SNV_newpipe"


def counts(text, label):
    m = re.search(rf"^{label} ([\d,]+) \(\{{'non': (\d+), 'somatic': (\d+), 'germline': (\d+)\}}\)", text, re.M)
    if not m:
        m = re.search(rf"{label} ([\d,]+) \(\{{'non': (\d+), 'somatic': (\d+), 'germline': (\d+)\}}\)", text)
    total, non, som, germ = int(m.group(1).replace(",", "")), *map(int, m.groups()[1:])
    assert total == non + som + germ, (label, total, non, som, germ)
    return dict(total=total, somatic=som, germline=germ, non=non)


def load(d):
    args = json.loads((d / "args.json").read_text())
    select = args.get("best_by") or "f1"
    log = (d / "train.log").read_text()
    params = float(re.search(r"\}, ([\d.]+) M parameters", log).group(1))
    tr = counts(log, "training tensors")
    va = counts(log, r"\), validation")
    m = re.search(r"each epoch ([\d,]+) tensors", log)
    per_epoch = int(m.group(1).replace(",", ""))
    best = int(re.search(r"^done; best epoch (\d+)", log, re.M).group(1))
    rows = [json.loads(x) for x in (d / "metrics.jsonl").read_text().splitlines() if x.strip()]
    ep = np.array([r["epoch"] for r in rows])
    assert list(ep) == list(range(1, len(ep) + 1)), d
    trl = np.array([r["train_loss"] for r in rows])
    vl = np.array([r["val"]["loss"] for r in rows])
    f1 = np.array([r["val"]["postpon"]["f1"] if select == "postpon" else r["val"]["thresholded"]["somatic"]["f1"]
                   for r in rows])
    # train.py keeps the first epoch with the highest score (ties: lower validation loss)
    order = sorted(range(len(ep)), key=lambda i: (-f1[i], vl[i], i))
    assert ep[order[0]] == best, (d, ep[order[0]], best)
    return dict(dir=d, select=select, params=params, train=tr, val=va, per_epoch=per_epoch, best=best,
                epochs=args["epochs"], ep=ep, train_loss=trl, val_loss=vl, f1=f1)


def shown(d):
    """Plotted epochs: 1 .. selected epoch + AFTER (capped at the last epoch)."""
    return min(d["best"] + AFTER, len(d["ep"]))


def fmt_k(n):
    return f"{n / 1e6:.2f} M" if n >= 1e6 else f"{n / 1e3:.1f} k"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ont-snv", choices=("letter", "newpipe"), default="letter",
                    help="ONT SNV six-tumor run: the letter's multicancer6_ONT_SNV (default) or the current-pipeline retrain")
    ap.add_argument("--outdir", type=Path, default=HERE)
    ap.add_argument("--height", type=float, default=6.0)
    args = ap.parse_args()
    data = {(k, p): load(run_dir(k, p, args.ont_snv)) for k in KINDS for p in PLATFORMS}
    for k in KINDS:
        assert len({data[(k, p)]["params"] for p in PLATFORMS}) == 1, k
        assert len({data[(k, p)]["epochs"] for p in PLATFORMS}) == 1, k

    W, H = 6.5, args.height
    LEG, HEAD, RGAP, FOOT, BLOCK = 0.30, 0.36, 0.10, 0.30, 0.16
    LEFT, RIGHT, CGAP = 0.55, 0.08, 0.30
    frac = 0.62  # loss row share of a block's plot height
    ph = (H - LEG - BLOCK - 2 * (HEAD + RGAP + FOOT)) / 2
    ph_loss, ph_f1 = ph * frac, ph * (1 - frac)
    pw = (W - LEFT - RIGHT - 2 * CGAP) / 3
    fig = plt.figure(figsize=(W, H))
    y = H - LEG
    for bi, kind in enumerate(KINDS):
        y -= HEAD
        d0 = data[(kind, "PacBio")]
        fig.text(0.06 / W, (y + 0.19) / H, "ab"[bi], fontsize=10, fontweight="bold", color=INK, ha="left", va="baseline")
        fig.text(0.23 / W, (y + 0.19) / H,
                 f"{'SNVs' if kind == 'SNV' else 'INDELs'}: {d0['params']:.1f} M-parameter network",
                 fontsize=8.5, color=INK, ha="left", va="baseline")
        for j, head in enumerate(HEADS):
            fig.text((LEFT + j * (pw + CGAP) + pw / 2) / W, (y + 0.045) / H, head, fontsize=7.5, color=INK,
                     ha="center", va="bottom")
        top_loss = y - ph_loss
        top_f1 = top_loss - RGAP - ph_f1
        for j, p in enumerate(PLATFORMS):
            d = data[(kind, p)]
            x0 = (LEFT + j * (pw + CGAP)) / W
            axl = fig.add_axes([x0, top_loss / H, pw / W, ph_loss / H])
            axf = fig.add_axes([x0, top_f1 / H, pw / W, ph_f1 / H], sharex=axl)
            b = d["best"] - 1
            for ax in (axl, axf):
                ax.axvline(d["best"], color=MUTED, lw=0.7, ls=(0, (1, 1.6)), zorder=1)
                for side in ("top", "right"):
                    ax.spines[side].set_visible(False)
            n = shown(d)  # epochs 1 .. selected + AFTER (user 2026-11-12: later epochs not drawn)
            ep, trl, vl, f1 = d["ep"][:n], d["train_loss"][:n], d["val_loss"][:n], d["f1"][:n]
            axl.plot(ep, trl, color=TRAIN_C, lw=1.3, zorder=3, solid_capstyle="round")
            axl.plot(ep, vl, color=VAL_C, lw=1.3, zorder=4, marker="o", ms=2.6, mec=SURFACE, mew=0.4)
            axl.plot(d["best"], d["val_loss"][b], marker="o", ms=4.6, color=VAL_C, mec=INK, mew=0.7, ls="", zorder=6)
            top = max(trl.max(), vl.max())
            axl.set_ylim(0, top * 1.12)
            axl.yaxis.set_major_locator(MaxNLocator(nbins=3, steps=[1, 2, 2.5, 5, 10]))
            axl.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            axl.tick_params(labelbottom=False)
            axl.text(d["best"] - 0.2, top * 1.10, f"epoch {d['best']}", fontsize=6, color=MUTED, ha="right", va="top")
            axf.plot(ep, f1, color=INK, lw=1.1, zorder=3, marker="o", ms=2.2, mec=SURFACE, mew=0.4)
            axf.plot(d["best"], d["f1"][b], marker="o", ms=4.6, color=INK, mec=SURFACE, mew=0.7, ls="", zorder=6)
            lo, hi = f1.min(), f1.max()
            axf.set_ylim(np.floor((lo - 0.01) * 20) / 20, np.ceil((hi + 0.01) * 20) / 20)
            axf.yaxis.set_major_locator(MaxNLocator(nbins=2, steps=[1, 2, 5, 10]))
            axf.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.2f}"))
            axf.set_xlim(0.5, n + 0.5)
            axf.set_xticks(list(range(1, n + 1)) if n <= 9 else [1] + list(range(2, n + 1, 2)))
            if d["select"] == "f1":
                axf.text(0.98, 0.06, "before the filter", transform=axf.transAxes, fontsize=5.8, color=MUTED,
                         ha="right", va="bottom")
            if j == 0:
                axl.set_ylabel("Loss", fontsize=7.5, labelpad=4)
                axf.set_ylabel("Validation\nsomatic F1", fontsize=7, labelpad=4)
            if j == 1:
                axf.set_xlabel("Epoch", fontsize=7.5, labelpad=2)
        y = top_f1 - FOOT - BLOCK
    h = [Line2D([], [], color=TRAIN_C, lw=1.4, label="Training loss"),
         Line2D([], [], color=VAL_C, lw=1.4, marker="o", ms=2.8, mec=SURFACE, mew=0.4, label="Validation loss"),
         Line2D([], [], color=INK, lw=1.1, marker="o", ms=2.4, mec=SURFACE, mew=0.4,
                label="Validation somatic F1 (selects the checkpoint)"),
         Line2D([], [], color=MUTED, lw=0.8, ls=(0, (1, 1.6)), label="Selected epoch")]
    fig.legend(handles=h, loc="upper center", ncol=4, frameon=False, fontsize=6.8, bbox_to_anchor=(0.5, 1 - 0.03 / H),
               handlelength=2.4, columnspacing=1.4, handletextpad=0.5)
    tag = "" if args.ont_snv == "letter" else "_ont_newpipe"
    for ext in ("png", "pdf"):
        fig.savefig(args.outdir / f"fig_R1.7{tag}.{ext}", dpi=300)
    plt.close(fig)

    plat = {"Illumina": "Illumina", "PacBio": "PacBio HiFi + Fiber-seq", "ONT": "ONT"}

    def rng(k, f, scale=None):
        xs = [f(data[(k, p)]) for p in PLATFORMS]
        if scale == "M":
            return f"{min(xs) / 1e6:.2f}–{max(xs) / 1e6:.2f} M"
        if scale == "M1":
            return f"{min(xs) / 1e6:.2f}–{max(xs) / 1e6:.1f} M"
        return f"{min(xs):,}–{max(xs):,}"

    gap = [abs(data[(k, p)]["val_loss"][data[(k, p)]["best"] - 1] - data[(k, p)]["train_loss"][data[(k, p)]["best"] - 1])
           for k in KINDS for p in PLATFORMS]
    ratio_best = [data[(k, p)]["val_loss"][data[(k, p)]["best"] - 1] / data[(k, p)]["val_loss"].min() for k in KINDS for p in PLATFORMS]
    s, d_ = "SNV", "INDEL"
    caption = (
        "**Supplementary Fig. R1.7 | Training curves of the six-tumor models on the HPRC v2.1 graph.** One run per platform, "
        "trained on chr2–22 of HCC1395T, HCC1954T, H1437T, H2009T, HG008T and COLO829T (chr1 of every tumor and HCC1395T chr20 "
        "held out). **a**, SNV models: "
        f"{data[(s, 'PacBio')]['params']:.1f} M parameters, trained for {data[(s, 'PacBio')]['epochs']} epochs; "
        f"{rng(s, lambda d: d['train']['somatic'])} somatic, {rng(s, lambda d: d['train']['germline'], 'M')} germline and "
        f"{rng(s, lambda d: d['train']['non'], 'M1')} other training examples. **b**, INDEL models: "
        f"{data[(d_, 'PacBio')]['params']:.1f} M parameters, trained for {data[(d_, 'PacBio')]['epochs']} epochs; "
        f"{rng(d_, lambda d: d['train']['somatic'])} somatic, {rng(d_, lambda d: d['train']['germline'], 'M')} germline and "
        f"{rng(d_, lambda d: d['train']['non'], 'M1')} other. Top rows: class-weighted cross-entropy on the training set "
        "(mean over each epoch) and on the validation set (5% of the 20,000-node blocks of the training chromosomes, "
        "disjoint from training). Bottom rows: validation somatic F1 after our population filter, the metric that selects "
        f"the checkpoint{' (ONT SNVs: before the filter)' if data[('SNV', 'ONT')]['select'] == 'f1' else ''}. "
        f"Dotted line: selected epoch; curves are shown up to {AFTER} epochs after it. At the selected epoch the validation loss is "
        f"{min(ratio_best):.2f}–{max(ratio_best):.2f} times the lowest of the run. Later epochs may overfit the training "
        "data, so their checkpoints are not used.")
    L = [f"# Supplementary Fig. R1.7{' (ONT SNV on the current pipeline)' if tag else ''}: six-tumor models on HPRC v2.1, training curves",
         "", "## Caption (draft)", "", caption, "", f"Generated by `make_figure.py{' --ont-snv newpipe' if tag else ''}` from each run's `train.log`, `args.json`, "
         "`metrics.jsonl` (paths below). Losses are class-weighted cross-entropy (training: mean over the epoch's "
         "batches; validation: after the epoch). Selection metric: `postpon` = validation somatic F1 after our population "
         "filter; `f1` = before it. Examples = candidate sites (tensors) on chr2-22 (HCC1395T without chr20); validation "
         "= 5% of the 20,000-node blocks, disjoint from training; chr1 is in neither.", "",
         "## Per model", "",
         "| type | platform | run | parameters | training somatic | training germline | training other (non) | training total | "
         "per epoch | validation somatic | validation germline | validation other | selection | selected epoch | "
         "train loss | val loss | val loss min (epoch) | val / min | plotted epochs | val loss last epoch | val F1 at selected | val F1 last |",
         "|" + "---|" * 22]
    for k in KINDS:
        for p in PLATFORMS:
            d = data[(k, p)]
            b = d["best"] - 1
            i = int(np.argmin(d["val_loss"]))
            t, v = d["train"], d["val"]
            L.append(f"| {k} | {plat[p]} | `{d['dir'].relative_to(RUNS)}` | {d['params']:.1f} M | {t['somatic']:,} | "
                     f"{t['germline']:,} | {t['non']:,} | {t['total']:,} | {d['per_epoch']:,} | {v['somatic']:,} | "
                     f"{v['germline']:,} | {v['non']:,} | {d['select']} | {d['best']} / {d['epochs']} | "
                     f"{d['train_loss'][b]:.3f} | {d['val_loss'][b]:.3f} | {d['val_loss'][i]:.3f} ({i + 1}) | "
                     f"{d['val_loss'][b] / d['val_loss'][i]:.2f} | 1–{shown(d)} | {d['val_loss'][-1]:.3f} | {d['f1'][b]:.3f} | {d['f1'][-1]:.3f} |")
    L += ["", "`per epoch`: SNV runs with `--non-fraction 0.25` draw the non examples with AF < 0.1 at 25% each epoch "
          "(the others every epoch); all other runs use every example every epoch.", "", "## Per epoch", ""]
    for k in KINDS:
        for p in PLATFORMS:
            d = data[(k, p)]
            L += [f"### {k} {plat[p]} (`{d['dir'].relative_to(RUNS)}`, selection `{d['select']}`)", "",
                  "| epoch | train loss | val loss | val somatic F1 |", "|---|---|---|---|"]
            for e, a, c, f in zip(d["ep"], d["train_loss"], d["val_loss"], d["f1"]):
                L.append(f"| {e}{' (selected)' if e == d['best'] else ' (not plotted)' if e > shown(d) else ''} | {a:.4f} | {c:.4f} | {f:.4f} |")
            L.append("")
    (args.outdir / f"data{tag}.md").write_text("\n".join(L))
    print(f"fig_R1.7{tag}.png/.pdf and data{tag}.md written to {args.outdir}")


if __name__ == "__main__":
    main()

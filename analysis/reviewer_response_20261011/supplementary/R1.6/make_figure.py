"""Supplementary Fig. R1.6 (R1.6: "full precision-recall curves for every caller").

chr1 precision-recall curves after our population filter, benchmark regions, HG008T, COLO829T and HCC1395T x three
platforms, (a) SNVs and (b) INDELs, for the models whose recall R1.6 / Supplementary Table R1.6 report and for both comparators:
- Pansoma SNV: HPRC v2.1 (d46) `--norm dataset` model, seed 0 only (<S>_<P>_SNV_d46_normbq, RIS copies in
  ris_7fc7396_20261012; PASS = validation recall 0.9 before the filter). Same runs as E8 latest_models.py.
  HCC1395T has no single-tumor SNV model: its row uses the latest six-tumor SNV models (see --snv six-tumor), trained
  on HCC1395T's chr2-22 among others.
- Pansoma INDEL: the final model (six tumors, HPRC v2.1), seed 0 only (release_s0_<P>_INDEL,
  vcf_ppr8_chr1_<set>_d46; PASS = ppr8).
  Curve = report.json pon.curve_roc (rtg vcfeval --all-records ranked by p_somatic; records down to p_somatic 0.01).
- DeepSomatic-TO, ClairS-TO (E8 score report.json, eval_bed): curve = all_pon (every record except RefCall / NoCall,
  whatever its FILTER, after our filter, ranked by QUAL); point = pass_pon (the tool's PASS after the same filter).
The figure shows only the curves and the PASS points. The two recall bounds (chr1 truth in the benchmark regions;
bounds.py -> bounds.json, which reproduces E27 denominators/label_status.json for HG008T / COLO829T, i.e. Supplementary
Table R1.6, and adds HCC1395T) are not drawn (user 2026-11-12); they stay as the last two columns of points.md:
  truth that is not a population allele (the bound for any caller after this filter);
  truth with an exact Pansoma candidate that is not a population allele (Pansoma's bound on this platform).
Writes fig_R1.6.png / .pdf and points.md (the plotted operating points and curve ends, for checking against the table).
Size: 6.5 x 7.0 in (the text width of a Letter / A4 Word page with 1-in margins, caption below on the same page; insert
at 100%), 300 dpi PNG; platform = column header, tumor = row label at the right (--height changes the height).
--snv six-tumor: the same figure with the six-tumor SNV model in (a): the six-tumor SNV models of the current SNV
  pipeline (user 2026-11-13: R1.6 follows the updated results), one run per platform, trained on chr2-22 of HCC1395T,
  HCC1954T, H1437T, H2009T, HG008T and COLO829T: snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, vcf_chr1_<set>_d46
  (site features, AF <= 0.1 non sampled 25% at weight 1, best epoch by the post-filter validation F1, PASS = rule 2:
  validation callable recall 0.95 after the filter; AF floor 0.05 on the CASTLE lines and COLO829T for PacBio HiFi /
  Fiber-seq / Illumina, 0.08 for ONT and HG008T). Test sets: HG008T and every ONT set at AF 0.08 (d46); COLO829T and
  HCC1395T PacBio HiFi / Fiber-seq / Illumina on their AF 0.05 sets (d46_snvaf05; bounds.json has them).
  One threshold per platform, the same in every tumor.
  (b) unchanged. Writes fig_R1.6_six_tumor.png / .pdf and points_six_tumor.md.
--snv six-tumor --four-tumor (user 2026-11-13: "R1.6 的图里能加上 4-tumor 的 results么"): the six-tumor figure plus one more
  Pansoma curve and PASS point per panel, the four-cell-line model (HCC1395T, HCC1954T, H1437T, H2009T; HG008T and
  COLO829T are tumors it never saw, HCC1395T chr1 is held out):
  SNV PacBio HiFi / Illumina multicancer4_<P>_SNV_sitef_af05_lnf25w1 with the rule-2 VCFs (vcf_chr1_<set>_d46_r095,
  COLO829T on its AF 0.05 set: _d46_snvaf05_r095) = the current SNV pipeline (same training arguments as the six-tumor
  runs; the checkpoint's call recall 0.9 is replaced by 0.95 in the VCF job). SNV ONT:
  snv_newpipe_20261109/multicancer4_ONT_SNV_newpipe (job 477096) once its VCFs exist, else multicancer4_ONT_SNV normfix
  (the OLD recipe; points file says which). INDEL: unseen_s0_<P>_INDEL, the four-cell-line model of the final recipe.
  Writes fig_R1.6_four_six_tumor.png / .pdf and points_four_six_tumor.md.
INDEL PASS point: ppr8 for Illumina / PacBio HiFi / Fiber-seq, pp (post-filter validation best F1) for ONT (user
  2026-11-13, six-tumor and four-cell-line models; vcf_pp_chr1_<set>_d46). The curve is the same for both rules.
Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_figure.py [--snv six-tumor [--four-tumor]]   (login node, < 1 min)
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
RIS = RUNS / "ris_7fc7396_20261012"
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
BOUNDS = HERE / "bounds.json"

TUMORS = ("HCC1395T", "COLO829T", "HG008T")  # row order top to bottom (user 2026-11-13)
PLATFORMS = ("Illumina", "PacBio", "ONT")
KINDS = ("SNV", "INDEL")
# validated categorical slots 1-3 (dataviz palette, light surface); comparators dashed as a second encoding
CALLERS = (("pansoma", "Pansoma", "#2a78d6", "-"), ("deepsomatic", "DeepSomatic-TO", "#eb6834", (0, (5, 2))),
           ("clairs_to", "ClairS-TO", "#1baf7a", (0, (5, 2))))
FOUR = "#e87ba4"  # the four-cell-line Pansoma model (--four-tumor), solid like the six-tumor one; pink, far from the six-tumor blue (user 2026-11-13)
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
                     "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
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


def pansoma_reports(tumor, p, kind, snv="single"):
    if kind == "SNV" and (snv == "six-tumor" or tumor == "HCC1395T"):
        fs = [RUNS / "snv_newpipe_20261109" / f"multicancer6_{p}_SNV_newpipe" / f"vcf_chr1_{ds(tumor, p)}_d46" / "eval_bed" / "report.json"]
        assert (p == "ONT" or tumor == "HG008T") == ("snvaf05" not in json.loads(fs[0].read_text())["calls"]), fs[0]
    elif kind == "SNV":
        fs = [RIS / f"{ds(tumor, p, True)}_SNV_d46_normbq{x}" / "vcf_chr1" / "eval_bed" / "report.json" for x in ("",)]  # seed 0 only (user 2026-11-10)
    else:
        fs = [RUNS / f"release_s{s}_{p}_INDEL" / f"vcf_{indel_rule(p)}_chr1_{ds(tumor, p)}_d46" / "eval_bed" / "report.json" for s in (0,)]  # seed 0 only
    return [(f, json.loads(f.read_text())["pon"]) for f in fs]


def indel_rule(p):
    return "pp" if p == "ONT" else "ppr8"  # user 2026-11-13: ONT INDEL PASS = post-filter validation best F1


def four_report(tumor, p, kind):
    """The four-cell-line model's report (file, pon block, on the current pipeline?)."""
    cur = True
    if kind == "INDEL":
        f = RUNS / f"unseen_s0_{p}_INDEL" / f"vcf_{indel_rule(p)}_chr1_{ds(tumor, p)}_d46"
    elif p == "ONT":
        f = RUNS / "snv_newpipe_20261109" / "multicancer4_ONT_SNV_newpipe" / f"vcf_chr1_{ds(tumor, p)}_d46"
        if not (f / "eval_bed" / "report.json").exists():
            f, cur = RUNS / "multicancer4_ONT_SNV" / f"vcf_chr1_{ds(tumor, p)}_d46_normfix", False
    else:
        f = RUNS / f"multicancer4_{p}_SNV_sitef_af05_lnf25w1" / f"vcf_chr1_{ds(tumor, p)}_d46{'_snvaf05' if tumor == 'COLO829T' else ''}_r095"
    f = f / "eval_bed" / "report.json"
    r = json.loads(f.read_text())
    if kind == "SNV" and p != "ONT":  # same test sets as the six-tumor model
        assert (tumor == "HG008T") == ("snvaf05" not in r["calls"]), f
    return f, r["pon"], cur


def comparator_report(caller, tumor, p, kind):
    f = E8 / caller / ds(tumor, p) / kind / "report.json"
    e = json.loads(f.read_text())["eval_bed"]
    return f, e["all_pon"], e["pass_pon"]


def bounds(status, tumor, p, kind, snv="single"):
    """(not a population allele, exact candidate and not a population allele) as shares of chr1 truth in the BED."""
    # the six-tumor PacBio HiFi / Illumina SNV models are tested on the AF 0.05 sets of COLO829T and HCC1395T (d46_snvaf05)
    six = kind == "SNV" and (snv == "six-tumor" or tumor == "HCC1395T")
    sfx = "_d46_snvaf05" if six and p != "ONT" and tumor != "HG008T" else "_d46"
    c = status[f"{ds(tumor, p)}{sfx}|{kind}"]
    return c["bound_any"], c["bound_exact"]


def span(xs):
    return f"{xs[0]:.3f}" if len(xs) == 1 else f"{min(xs):.3f}–{max(xs):.3f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snv", choices=("single", "six-tumor"), default="single")
    ap.add_argument("--four-tumor", action="store_true", help="with --snv six-tumor: add the four-cell-line model's curve and PASS point")
    ap.add_argument("--outdir", type=Path, default=HERE)
    ap.add_argument("--height", type=float, default=7.0, help="figure height in inches (width is 6.5, the Word text width)")
    args = ap.parse_args()
    assert not args.four_tumor or args.snv == "six-tumor", "--four-tumor needs --snv six-tumor"
    tag = "" if args.snv == "single" else "_four_six_tumor" if args.four_tumor else "_six_tumor"
    status = json.loads(BOUNDS.read_text())
    data, rows = {}, []
    for kind in KINDS:
        for tumor in TUMORS:
            for p in PLATFORMS:
                d = dict(pansoma=pansoma_reports(tumor, p, kind, args.snv),
                         **{c: comparator_report(c, tumor, p, kind) for c, *_ in CALLERS[1:]},
                         bounds=bounds(status, tumor, p, kind, args.snv))
                if args.four_tumor:
                    d["four"] = four_report(tumor, p, kind)
                    n = lambda b: b["at_threshold"]["tp_baseline"] + b["at_threshold"]["fn"]
                    assert n(d["four"][1]) == n(d["pansoma"][0][1]), (kind, tumor, p)  # same truth
                data[(kind, tumor, p)] = d
    xmax = {}  # per row: SNV 0-1; INDEL to just past the furthest curve of that tumor
    for kind in KINDS:
        for tumor in TUMORS:
            m = 0
            for (k, t, _), d in data.items():
                if (k, t) != (kind, tumor):
                    continue
                m = max([m] + [curve(b)[:, 0].max() for _, b in d["pansoma"]] + ([curve(d["four"][1])[:, 0].max()] if "four" in d else []) +
                        [curve(d[c][1])[:, 0].max() for c, *_ in CALLERS[1:]])
            xmax[(kind, tumor)] = 1.0 if kind == "SNV" else min(1.0, np.ceil((m + 0.03) * 10) / 10)

    # Word page layout (Letter / A4 with 1-in margins: text block 6.5 x 9 in). 6.5 in wide at true size; 7.0 in tall
    # leaves room for the caption on the same page. Facet grid: platform = column header, tumor = row label at the
    # right, axis labels once per block, so no per-panel titles; positions are set in inches (no tight_layout).
    nt = len(TUMORS)
    W, H = 6.5, args.height
    LEG, HEAD, GAP, FOOT, BLOCK = 0.27, 0.36, 0.27, 0.34, 0.12  # legend, block header, gap between rows, x label, gap between blocks
    LEFT, RIGHT, CGAP = 0.50, 0.24, 0.25                        # y label + ticks, row labels, gap between columns
    ph = (H - LEG - BLOCK - 2 * (HEAD + (nt - 1) * GAP + FOOT)) / (2 * nt)
    pw = (W - LEFT - RIGHT - 2 * CGAP) / 3
    heads = ("Illumina", "PacBio HiFi (COLO829T: Fiber-seq)", "ONT")
    fig = plt.figure(figsize=(W, H))
    y = H - LEG
    for bi, kind in enumerate(KINDS):
        y -= HEAD
        fig.text(0.06 / W, (y + 0.19) / H, "ab"[bi], fontsize=10, fontweight="bold", color=INK, ha="left", va="baseline")
        fig.text(0.23 / W, (y + 0.19) / H, "SNVs" if kind == "SNV" else "INDELs", fontsize=8.5, color=INK, ha="left", va="baseline")
        for j, head in enumerate(heads):
            fig.text((LEFT + j * (pw + CGAP) + pw / 2) / W, (y + 0.045) / H, head, fontsize=7.5, color=INK, ha="center", va="bottom")
        for ti, tumor in enumerate(TUMORS):
            bottom = y - ph
            for j, p in enumerate(PLATFORMS):
                ax = fig.add_axes([(LEFT + j * (pw + CGAP)) / W, bottom / H, pw / W, ph / H])
                d = data[(kind, tumor, p)]
                X = xmax[(kind, tumor)]
                b_any, b_exact = d["bounds"]  # points.md only, not drawn
                key, label, color, ls = CALLERS[0]
                prs = []
                for f, blk in d["pansoma"]:
                    c = curve(blk)
                    ax.plot(c[:, 0], c[:, 1], color=color, lw=1.2, ls=ls, zorder=3, solid_capstyle="round")
                    r, pr = point(blk)
                    prs.append((r, pr, c[:, 0].max()))
                for r, pr, _ in prs:
                    ax.plot(r, pr, marker="o", ms=4.0, color=color, mec=SURFACE, mew=0.6, ls="", zorder=6)
                cells = [span([x[k] for x in prs]) for k in range(3)]
                if "four" in d:
                    _, blk, cur = d["four"]
                    c = curve(blk)
                    ax.plot(c[:, 0], c[:, 1], color=FOUR, lw=1.1, ls="-", zorder=2.5, solid_capstyle="round")
                    r, pr = point(blk)
                    ax.plot(r, pr, marker="o", ms=4.0, color=FOUR, mec=SURFACE, mew=0.6, ls="", zorder=5.5)
                    cells += [f"{r:.3f}", f"{pr:.3f}", f"{c[:, 0].max():.3f}" + ("" if cur else " †")]
                for key, label, color, ls in CALLERS[1:]:
                    _, allp, passp = d[key]
                    c = curve(allp)
                    ax.plot(c[:, 0], c[:, 1], color=color, lw=1.05, ls=(0, (4, 1.6)), zorder=4)
                    r, pr = point(passp)
                    ax.plot(r, pr, marker="o", ms=4.0, color=color, mec=SURFACE, mew=0.6, ls="", zorder=7)
                    cells += [f"{r:.3f}", f"{pr:.3f}", f"{c[:, 0].max():.3f}"]
                rows.append(f"| {kind} | {tumor} | {plat(tumor, p)} | " + " | ".join(cells) +
                            f" | {b_exact:.3f} | {b_any:.3f} |")
                ax.set_xlim(0, X)
                ax.set_ylim(0, 1.03)
                ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
                ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
                ax.set_yticks([0, 0.5, 1])
                ax.set_yticks([0.25, 0.75], minor=True)
                ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
                for side in ("top", "right"):
                    ax.spines[side].set_visible(False)
                if ti == nt - 1 and j == 1:
                    ax.set_xlabel("Recall (chr1)", fontsize=7.5, labelpad=3)
                if j == 0:
                    if ti == nt // 2:
                        ax.set_ylabel("Precision (chr1)", fontsize=7.5, labelpad=4)
                else:
                    ax.tick_params(labelleft=False)
            fig.text((W - RIGHT + 0.07) / W, (bottom + ph / 2) / H, tumor, fontsize=7.5, color=INK, rotation=270,
                     ha="left", va="center")
            y = bottom - (GAP if ti < nt - 1 else FOOT)
        y -= BLOCK
    h = [Line2D([], [], color=c, lw=1.3, ls="-" if k == "pansoma" else (0, (4, 1.6)), label=l) for (k, l, c, _) in CALLERS]
    if args.four_tumor:
        h[0].set_label("Pansoma, six tumors")
        h.insert(1, Line2D([], [], color=FOUR, lw=1.3, ls="-", label="Pansoma, four cell lines"))
    h += [Line2D([], [], color=MUTED, marker="o", ls="", ms=4, mec=SURFACE, mew=0.6, label="PASS (after the filter)")]
    fig.legend(handles=h, loc="upper center", ncol=len(h), frameon=False, fontsize=6.8, bbox_to_anchor=(0.5, 1 - 0.03 / H),
               handlelength=2.6 if len(h) == 4 else 2.0, columnspacing=1.6 if len(h) == 4 else 1.1, handletextpad=0.6, labelspacing=0.45)
    for ext in ("png", "pdf"):
        fig.savefig(args.outdir / f"fig_R1.6{tag}.{ext}", dpi=300)
    plt.close(fig)

    m6 = ("the six-tumor model of the current SNV pipeline, one run per platform (`snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe`, "
          "rule 2; COLO829T / HCC1395T PacBio HiFi / Fiber-seq / Illumina on their AF 0.05 sets)")
    snv = (f"SNV: {m6}" if args.snv == "six-tumor"
           else f"SNV: single-tumor `--norm dataset` HPRC v2.1 model, seed 0 (HCC1395T: {m6}; no single-tumor HCC1395T SNV model)")
    title = {"": "", "_six_tumor": " (six-tumor SNV model)", "_four_six_tumor": " (six-tumor and four-cell-line models)"}[tag]
    cmd = {"": "", "_six_tumor": " --snv six-tumor", "_four_six_tumor": " --snv six-tumor --four-tumor"}[tag]
    four = ("" if not args.four_tumor else " Four-cell-line model (HCC1395T, HCC1954T, H1437T, H2009T; HG008T and COLO829T unseen): SNV "
            "`multicancer4_<P>_SNV_sitef_af05_lnf25w1` at rule 2 (`_r095`; current SNV pipeline), ONT `multicancer4_ONT_SNV_newpipe` when "
            "its VCFs exist, else the old-recipe `multicancer4_ONT_SNV` normfix (marked †); INDEL `unseen_s0_<P>_INDEL`.")
    L = [f"# Supplementary Fig. R1.6{title}: plotted points (chr1, benchmark regions, after our population filter)", "",
         f"Generated by `make_figure.py{cmd}`. Pansoma: {snv}; INDEL: the final six-tumor model, seed 0, PASS at "
         f"ppr8 (ONT: pp, the post-filter validation best-F1 threshold).{four} Comparators: PASS after the same filter; curve end = highest recall of the "
         "QUAL-ranked curve over every record except RefCall / NoCall after the filter.", "",
         "| kind | tumor | platform | Pansoma PASS R | Pansoma PASS P | Pansoma curve end R | " +
         ("four-cell-line PASS R | P | curve end R | " if args.four_tumor else "") + "DeepSomatic-TO PASS R | P | "
         "curve end R | ClairS-TO PASS R | P | curve end R | exact candidate, not population | not population (any caller) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|" + ("---|---|---|" if args.four_tumor else "")] + rows
    (args.outdir / f"points{tag}.md").write_text("\n".join(L) + "\n")
    print(f"fig_R1.6{tag}.png/.pdf and points{tag}.md written to {args.outdir}")


if __name__ == "__main__":
    main()

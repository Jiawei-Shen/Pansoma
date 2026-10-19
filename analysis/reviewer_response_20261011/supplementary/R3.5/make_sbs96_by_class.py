"""Supplementary material of R3.5: SBS96 spectra of the graph-coordinate SNV calls, by validation class.

Graph-coordinate calls are the PASS SNV calls on nodes without a GRCh38 position (`*.linear.unplaced.vcf.gz`); their
trinucleotide context is taken from the graph: the REF graph-path context that E9c spelled from each call's tensor
(`ref_graph`, 50 bp each side; the ALT graph-path context `alt_graph` with the REF base at the site when the REF
context is missing). Validation classes are E9c / E9d's exclusive categories (`per_call.tsv.gz`, R3.3): population
(graph-native PoN: the ALT context is spelled by GRCh38 or an HPRC haplotype or walk) > germline (in the matched
normal's assembly) > tumor-only (HG008T: in HG008-T, not in HG008-N) / not in normal (COLO829T: no tumor assembly) >
unverifiable > unsupported.
Runs: the nine E9c in-sample SNV runs on a pangenome graph (HG008T on HPRC v1.1 d9 and v2.1 d46, COLO829T on d46,
three platforms each; COLO829T ONT is E9c's random 3,000 of 4,569). References: the chr1 truth and matched-normal
germline spectra of E18 (`experiments/E18_sbs96/sbs96.tsv`) and COSMIC v3.4 (E18's loader).
Writes sbs96_by_class.tsv (96 counts per run x class and per sample x class), cosine_by_class.tsv and
fig_R3.5_sbs96_by_class.png next to this script.
usage: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_sbs96_by_class.py
"""
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1] / "experiments"
sys.path.insert(0, str(EXP / "E18_sbs96"))
import e18_sbs96 as E18  # noqa: E402

PER_CALL = EXP / "E9c_offref_filter_validation" / "per_call.tsv.gz"
RUNS = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation/runs")
E18_TSV = EXP / "E18_sbs96" / "sbs96.tsv"
CLASS_ORDER = ["population", "germline", "tumor_only", "not_in_normal", "other"]
CLASS_LABEL = {"population": "population allele", "germline": "germline (normal assembly)",
               "tumor_only": "tumor-only (HG008-T, not HG008-N)", "not_in_normal": "not in the normal (COLO829BL)",
               "other": "unverifiable / unsupported"}
CHANNELS = E18.default_channels()


def sample_of(run):
    return "HG008T" if run.startswith("rcpu__HG008") else "COLO829T"


def tri_key(row):
    """SBS96 channel from the graph-path context of a call, or None."""
    ref, alt = row["evref"], row["evalt"]
    if len(ref) != 1 or len(alt) != 1 or ref == alt or alt not in "ACGT":
        return None
    for src in ("ref_graph", "alt_graph"):
        ctx, site = row.get(src) or "", row.get(f"{src}_site") or ""
        if not ctx or "-" not in site:
            continue
        s = int(site.split("-")[0])
        if s < 1 or s + 1 >= len(ctx):
            continue
        tri = (ctx[s - 1] + ref + ctx[s + 1]).upper()
        if src == "ref_graph" and ctx[s].upper() != ref:
            continue
        if any(b not in "ACGT" for b in tri):
            continue
        if ref in "GA":
            tri, ref2, alt2 = tri.translate(E18.COMP)[::-1], ref.translate(E18.COMP), alt.translate(E18.COMP)
        else:
            ref2, alt2 = ref, alt
        return f"{tri[0]}[{ref2}>{alt2}]{tri[2]}"
    return None


def main():
    cls_of = {}
    for r in csv.DictReader(gzip.open(PER_CALL, "rt"), delimiter="\t"):
        if r["run"].startswith("rcpu__") and r["group"] == "unplaced_pass" and r["kind"] == "SNP":
            c = r["category"] if r["category"] in CLASS_ORDER else "other"
            cls_of[(r["run"], r["candidate_id"])] = c
    runs = sorted({k[0] for k in cls_of})
    per_run, per_sample, skipped = defaultdict(Counter), defaultdict(Counter), Counter()
    for run in runs:
        with open(RUNS / run / "calls.tsv") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row["group"] != "unplaced_pass" or (run, row["candidate_id"]) not in cls_of:
                    continue
                c = cls_of[(run, row["candidate_id"])]
                k = tri_key(row)
                if k is None:
                    skipped[(run, c)] += 1
                    continue
                per_run[(run, c)][k] += 1
                per_sample[(sample_of(run), c)][k] += 1

    ref = {}
    with open(E18_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["spectrum"] in ("truth COLO829T", "truth HG008T", "germline COLO829T", "germline HG008T"):
                ref[row["spectrum"]] = Counter({k: int(row[k]) for k in CHANNELS if int(row[k])})
    _, _, sigs = E18.load_cosmic()

    with open(HERE / "sbs96_by_class.tsv", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["level", "sample", "run", "class", "n", "skipped_no_context"] + CHANNELS)
        for (run, c), cnt in sorted(per_run.items()):
            w.writerow(["run", sample_of(run), run.replace("rcpu__", "").replace("__vcf_chr1", ""), c,
                        sum(cnt.values()), skipped[(run, c)]] + [cnt.get(k, 0) for k in CHANNELS])
        for (smp, c), cnt in sorted(per_sample.items()):
            sk = sum(v for (run, cc), v in skipped.items() if cc == c and sample_of(run) == smp)
            w.writerow(["sample", smp, "all runs", c, sum(cnt.values()), sk] + [cnt.get(k, 0) for k in CHANNELS])

    rows = []
    for smp in ("HG008T", "COLO829T"):
        for c in CLASS_ORDER:
            cnt = per_sample.get((smp, c))
            if not cnt:
                continue
            h = E18.hallmarks(cnt)
            best = sorted(((E18.cosine(cnt, {k: v for k, v in s.items()}), name) for name, s in sigs.items()),
                          reverse=True)[:3] if sigs else []
            mix = E18.mix_truth_germline(cnt, ref[f"truth {smp}"], ref[f"germline {smp}"], CHANNELS)
            rows.append(dict(sample=smp, cls=c, n=h["n"], cos_truth=E18.cosine(cnt, ref[f"truth {smp}"]),
                             cos_germline=E18.cosine(cnt, ref[f"germline {smp}"]),
                             truth_like=mix[0] if mix else None, uv=h["uv"], cpg=h["cpg"],
                             best=", ".join(f"{n} {v:.3f}" for v, n in best)))
    for smp in ("HG008T", "COLO829T"):
        t, g = ref[f"truth {smp}"], ref[f"germline {smp}"]
        for lab, cnt in ((f"truth {smp} (reference)", t), (f"germline {smp} (reference)", g)):
            h = E18.hallmarks(cnt)
            rows.append(dict(sample=smp, cls=lab, n=h["n"], cos_truth=E18.cosine(cnt, t), cos_germline=E18.cosine(cnt, g),
                             truth_like=None, uv=h["uv"], cpg=h["cpg"], best=""))
    with open(HERE / "cosine_by_class.tsv", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "class", "n", "cosine_to_truth", "cosine_to_normal_germline", "truth_like_share_nnls",
                    "uv_share", "cpg_share", "best_cosmic"])
        for r in rows:
            f3 = lambda x: "" if x is None else f"{x:.3f}"  # noqa: E731
            w.writerow([r["sample"], r["cls"], r["n"], f3(r["cos_truth"]), f3(r["cos_germline"]), f3(r["truth_like"]),
                        f3(r["uv"]), f3(r["cpg"]), r["best"]])
    figure(per_sample, ref, rows)
    for r in rows:
        print(r["sample"], r["cls"], r["n"], "cos truth %s germline %s" % (r["cos_truth"], r["cos_germline"]))


def figure(per_sample, ref, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # COSMIC's conventional SBS96 class colours.
    col = {"C>A": "#1ebff0", "C>G": "#050708", "C>T": "#e62725", "T>A": "#cbcacb", "T>C": "#a1cf64", "T>G": "#edc8c5"}
    cos = {(r["sample"], r["cls"]): r for r in rows}
    panels = []
    for smp in ("COLO829T", "HG008T"):
        panels.append((smp, f"truth {smp}", f"{smp}: chr1 somatic truth", ref[f"truth {smp}"], None))
        panels.append((smp, f"germline {smp}", f"{smp}: matched-normal germline", ref[f"germline {smp}"], None))
        for c in CLASS_ORDER:
            cnt = per_sample.get((smp, c))
            if cnt and sum(cnt.values()) >= 10:
                panels.append((smp, c, f"{smp} graph-coordinate calls: {CLASS_LABEL[c]}", cnt, cos.get((smp, c))))
    fig, axs = plt.subplots(len(panels) + 1, 1, figsize=(11, 1.55 * len(panels) + 0.3), sharex=True,
                            gridspec_kw=dict(height_ratios=[0.18] + [1] * len(panels)))
    head, axes = axs[0], axs[1:]
    head.axis("off")
    for i, c in enumerate(E18.CLASSES):
        head.axvspan(16 * i, 16 * i + 15, ymin=0.0, ymax=0.35, color=col[c])
        head.text(16 * i + 7.5, 0.75, c, ha="center", va="center", fontsize=9, fontweight="bold")
    for ax, (smp, key, title, cnt, r) in zip(axes, panels):
        n = sum(cnt.values())
        y = [cnt.get(k, 0) / n * 100 for k in CHANNELS]
        ax.bar(range(96), y, color=[col[k[2:5]] for k in CHANNELS], width=0.8)
        sub = f"n = {n:,}"
        if r:
            sub += f"; cosine to truth {r['cos_truth']:.2f}, to normal germline {r['cos_germline']:.2f}"
        ax.set_title(f"{title} ({sub})", fontsize=8.5, loc="left")
        ax.set_ylabel("% of SNVs", fontsize=7.5)
        ax.tick_params(axis="y", labelsize=7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for i in range(1, 6):
            ax.axvline(16 * i - 0.5, color="#bbbbbb", lw=0.5)
    axes[-1].set_xticks(range(96))
    axes[-1].set_xticklabels([k[0] + k[2] + k[6] for k in CHANNELS], rotation=90, fontsize=5.5, family="monospace")
    head.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(HERE / "fig_R3.5_sbs96_by_class.png", dpi=200)


if __name__ == "__main__":
    main()

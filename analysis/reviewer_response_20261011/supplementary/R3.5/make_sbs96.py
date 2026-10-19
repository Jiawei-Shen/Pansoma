"""R3.5: SBS96 spectra for Figure 4a and the Supplement (chr1, one call set per caller, three platforms separately).

a) Figure 4a: the chr1 truth spectrum and the matched normal's germline spectrum (its germline truth on chr1 in the
   same benchmark regions, E18) beside the PASS SNV calls after population filtering of Pansoma (HPRC v2.1, d46
   filter), DeepSomatic-TO and ClairS-TO, in the benchmark regions. Counts from E8 (`experiments/E8_comparators/
   sbs96.tsv`: truth, and TP + FP of each call set after our PoN; the comparators are the user's DeepSomatic 1.9.0 and
   ClairS-TO 0.4.1 / 0.4.2 tumor-only call sets).
b) Supplement: Pansoma's graph-coordinate SNV PASS calls (no GRCh38 position, so no GRCh38 PoN and no benchmark region
   applies; whole chr1), same runs. Trinucleotide context from the graph: the bases flanking the site in the REF
   graph-path context that E9c spelled from each call's reads (`ref_graph`; `alt_graph` with the REF base when missing).
   COLO829T ONT: E9c's random 3,000 of its 4,569 calls.
Cosine similarities to the chr1 truth and, for b, to the matched normal's germline SNVs (E18 `sbs96.tsv`).
Writes sbs96.tsv, cosine.tsv, fig_R3.5a_<tumor>.png and fig_R3.5b_graph_coordinate.png next to this script.
usage: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_sbs96.py
"""
import csv
import math
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1] / "experiments"
E8 = EXP / "E8_comparators" / "sbs96.tsv"
E18 = EXP / "E18_sbs96" / "sbs96.tsv"
E9C_RUNS = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation/runs")
COMP = str.maketrans("ACGT", "TGCA")
CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
CHANNELS = [f"{a}[{s}]{b}" for s in CLASSES for a in "ACGT" for b in "ACGT"]
COLOR = {"C>A": "#1ebff0", "C>G": "#050708", "C>T": "#e62725", "T>A": "#cbcacb", "T>C": "#a1cf64", "T>G": "#edc8c5"}
TUMORS = {"HG008T": [("Illumina", "HG008T_Illumina", "HG008_Illumina_SNV_d46_7fc7396"),
                     ("PacBio HiFi", "HG008T_PacBio", "HG008_PacBio_SNV_d46_7fc7396"),
                     ("ONT", "HG008T_ONT", "HG008_ONT_SNV_d46_7fc7396")],
          "COLO829T": [("Illumina", "COLO829T_Illumina", "COLO829_Illumina_SNV_d46_7fc7396"),
                       ("Fiber-seq", "COLO829T_fiberseq", "COLO829_fiberseq_SNV_d46_7fc7396"),
                       ("ONT", "COLO829T_ONT", "COLO829_ONT_SNV_d46_7fc7396")]}
CALLERS = [("Pansoma", "Pansoma HPRC v2.1 d46"), ("DeepSomatic-TO", "DeepSomatic 1.9.0"), ("ClairS-TO", "ClairS-TO 0.4.1 / 0.4.2")]


def cosine(a, b):
    dot = sum(a.get(k, 0) * b.get(k, 0) for k in CHANNELS)
    na, nb = math.sqrt(sum(a.get(k, 0) ** 2 for k in CHANNELS)), math.sqrt(sum(b.get(k, 0) ** 2 for k in CHANNELS))
    return dot / (na * nb) if na and nb else float("nan")


def uv(c):
    n = sum(c.values())
    return sum(v for k, v in c.items() if k[2:5] == "C>T" and k[0] in "CT") / n if n else float("nan")


def e8_spectra():
    """dataset -> {'truth': Counter, caller: Counter (TP + FP after the PoN)}."""
    out = {}
    with open(E8) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            c = Counter({k: int(row[k]) for k in CHANNELS if int(row[k])})
            d = out.setdefault(row["dataset"], {})
            if row["call_set"] == "truth":
                d["truth"] = c
                continue
            for label, prefix in CALLERS:
                if row["call_set"].startswith(prefix) and row["call_set"].endswith((" TP", " FP")):
                    d.setdefault(label, Counter()).update(c)
    return out


def germline_spectra():
    out = {}
    with open(E18) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["spectrum"] in ("germline HG008T", "germline COLO829T"):
                out[row["spectrum"].split()[1]] = Counter({k: int(row[k]) for k in CHANNELS if int(row[k])})
    return out


def graph_coordinate(run):
    """SBS96 of the unplaced PASS SNV calls of one run, context from E9c's graph-path REF context."""
    c, skipped = Counter(), 0
    with open(E9C_RUNS / f"rcpu__{run}__vcf_chr1" / "calls.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["group"] != "unplaced_pass" or row["kind_event"] != "SNP":
                continue
            ref, alt, key = row["evref"], row["evalt"], None
            for src in ("ref_graph", "alt_graph"):
                ctx, site = row.get(src) or "", row.get(f"{src}_site") or ""
                if not ctx or "-" not in site:
                    continue
                s = int(site.split("-")[0])
                if s < 1 or s + 1 >= len(ctx) or (src == "ref_graph" and ctx[s].upper() != ref):
                    continue
                tri = (ctx[s - 1] + ref + ctx[s + 1]).upper()
                if any(b not in "ACGT" for b in tri) or alt not in "ACGT" or alt == ref:
                    continue
                if ref in "GA":
                    tri, ref2, alt2 = tri.translate(COMP)[::-1], ref.translate(COMP), alt.translate(COMP)
                else:
                    ref2, alt2 = ref, alt
                key = f"{tri[0]}[{ref2}>{alt2}]{tri[2]}"
                break
            if key:
                c[key] += 1
            else:
                skipped += 1
    return c, skipped


def panel(ax, c, title, pad=4):
    n = sum(c.values())
    ax.bar(range(96), [c.get(k, 0) / n * 100 if n else 0 for k in CHANNELS], color=[COLOR[k[2:5]] for k in CHANNELS],
           width=0.8)
    ax.set_title(title, fontsize=8, loc="left", pad=pad)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for i in range(1, 6):
        ax.axvline(16 * i - 0.5, color="#bbbbbb", lw=0.5)
    ax.tick_params(axis="y", labelsize=6.5)
    ax.set_xticks([])


def class_header(ax):
    for i, c in enumerate(CLASSES):
        ax.text(16 * i + 7.5, 1.0, c, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7,
                fontweight="bold", color=COLOR[c] if c not in ("T>A", "T>G") else "#555555")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    e8, germ = e8_spectra(), germline_spectra()
    tsv, cos = [], []
    for tumor, plats in TUMORS.items():
        fig, axes = plt.subplots(5, 3, figsize=(15, 9.0), sharey="row")
        for j, (plat, ds, run) in enumerate(plats):
            truth = e8[ds]["truth"]
            rows = [("chr1 somatic truth", truth), ("matched-normal germline", germ[tumor])] + \
                   [(lab, e8[ds][lab]) for lab, _ in CALLERS]
            for i, (lab, c) in enumerate(rows):
                n = sum(c.values())
                if i == 0:
                    t = f"{tumor} {plat}: {lab} (n = {n:,})"
                elif i == 1:
                    t = f"{lab} (n = {n:,}; cosine to truth {cosine(c, truth):.3f})"
                else:
                    t = f"{lab} (n = {n:,}; cosine to truth {cosine(c, truth):.3f}, to germline {cosine(c, germ[tumor]):.3f})"
                panel(axes[i, j], c, t, pad=14 if i == 0 else 4)
                tsv.append(["Figure 4a", tumor, plat, lab, n] + [c.get(k, 0) for k in CHANNELS])
                if i > 1:
                    cos.append(["Figure 4a", tumor, plat, lab, n, f"{cosine(c, truth):.3f}",
                                f"{cosine(c, germ[tumor]):.3f}", f"{uv(c):.3f}"])
            class_header(axes[0, j])
        for i in range(5):
            axes[i, 0].set_ylabel("% of SNVs", fontsize=7)
        for j in range(3):
            axes[4, j].set_xticks(range(96))
            axes[4, j].set_xticklabels([k[0] + k[2] + k[6] for k in CHANNELS], rotation=90, fontsize=3.8, family="monospace")
        fig.suptitle(f"{tumor}, chr1, PASS SNVs after population filtering, benchmark regions", fontsize=9, x=0.01, ha="left")
        fig.tight_layout()
        fig.savefig(HERE / f"fig_R3.5a_{tumor}.png", dpi=200)
        plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15, 4.2))
    for i, (tumor, plats) in enumerate(TUMORS.items()):
        truth = e8[plats[0][1]]["truth"]
        for j, (plat, ds, run) in enumerate(plats):
            c, skipped = graph_coordinate(run)
            n = sum(c.values())
            panel(axes[i, j], c, f"{tumor} {plat} (n = {n:,}; cosine to truth {cosine(c, truth):.2f}, "
                                 f"to normal germline {cosine(c, germ[tumor]):.2f})", pad=14 if i == 0 else 4)
            tsv.append(["Supplement graph-coordinate", tumor, plat, "Pansoma graph-coordinate", n] + [c.get(k, 0) for k in CHANNELS])
            cos.append(["Supplement graph-coordinate", tumor, plat, "Pansoma graph-coordinate", n,
                        f"{cosine(c, truth):.3f}", f"{cosine(c, germ[tumor]):.3f}", f"{uv(c):.3f}"])
            print(tumor, plat, "graph-coordinate", n, "skipped", skipped)
        axes[i, 0].set_ylabel("% of SNVs", fontsize=7)
    for j in range(3):
        class_header(axes[0, j])
        axes[1, j].set_xticks(range(96))
        axes[1, j].set_xticklabels([k[0] + k[2] + k[6] for k in CHANNELS], rotation=90, fontsize=3.8, family="monospace")
    fig.suptitle("Pansoma graph-coordinate SNV PASS calls, chr1 (HPRC v2.1, d46 filter); context from the graph path",
                 fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(HERE / "fig_R3.5b_graph_coordinate.png", dpi=200)
    plt.close(fig)
    for tumor in TUMORS:
        t = e8[TUMORS[tumor][0][1]]["truth"]
        cos.append(["reference", tumor, "-", "matched-normal germline", sum(germ[tumor].values()),
                    f"{cosine(germ[tumor], t):.3f}", "1.000", f"{uv(germ[tumor]):.3f}"])
        cos.append(["reference", tumor, "-", "chr1 somatic truth", sum(t.values()), "1.000",
                    f"{cosine(t, germ[tumor]):.3f}", f"{uv(t):.3f}"])

    with open(HERE / "sbs96.tsv", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["use", "tumor", "platform", "spectrum", "n"] + CHANNELS)
        w.writerows(tsv)
    with open(HERE / "cosine.tsv", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["use", "tumor", "platform", "spectrum", "n", "cosine_to_truth", "cosine_to_normal_germline", "uv_type_ct_share"])
        w.writerows(cos)
    for r in cos:
        print("\t".join(map(str, r)))


if __name__ == "__main__":
    sys.exit(main())

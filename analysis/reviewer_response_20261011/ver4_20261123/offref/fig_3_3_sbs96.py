"""ver4 offref (4): Response Figure 3.3, SBS96 spectra of Pansoma's graph-coordinate SNV PASS calls of the nine
six-tumor call sets (snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, HPRC v2.1, chr1, default output), one panel
per tumor and platform (3 x 3).

New script; the spectrum code (graph_coordinate, cosine, uv) is copied unchanged from supplementary/R3.5/make_sbs96.py
(which is NOT run: its main() would rewrite the Figure 4a files from the October models). Differences: nine panels
(HCC1395T added), the calls of the ver4 E9c runs (<X4>/DATA/runs/m6np__<SET>__vcf_chr1/calls.tsv, every call processed,
no sampling), truth and matched-normal germline spectra of all three tumors from experiments/E18_sbs96/sbs96.tsv
(`truth <T>`: chr1 somatic truth SNVs in the benchmark regions; `germline <T>`), and the panel title gives n and the
cosine similarity to the tumor's chr1 truth spectrum (the cosine to the matched normal's germline spectrum goes to
cosine.tsv). Trinucleotide context from the graph: the bases flanking the site in the REF graph-path context that E9c
spelled from each call's tensor (`ref_graph`; `alt_graph` with the REF base when missing). A call without a usable
context is left out and counted.

Writes fig_3.3_sbs96_graph_coordinate.png (15 x 9 in at 200 dpi = 3000 x 1800 px), fig_3.3_sbs96.tsv (counts per
channel), fig_3.3_cosine.tsv.
usage: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python fig_3_3_sbs96.py [--repro]
       --repro: the six ver3 panels' n and cosines recomputed with this script from the October E9c runs and compared
                with supplementary/R3.5/cosine.tsv (writes nothing)
"""
import csv
import math
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
A = HERE.parents[1]
E18 = A / "experiments" / "E18_sbs96" / "sbs96.tsv"
E9C_RUNS = Path("/scratch/jshen/data/reviewer_experiments_20261013/ver4_20261123/offref/DATA/runs")
OLD_RUNS = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation/runs")
OLD_COS = A / "supplementary" / "R3.5" / "cosine.tsv"
COMP = str.maketrans("ACGT", "TGCA")
CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
CHANNELS = [f"{a}[{s}]{b}" for s in CLASSES for a in "ACGT" for b in "ACGT"]
COLOR = {"C>A": "#1ebff0", "C>G": "#050708", "C>T": "#e62725", "T>A": "#cbcacb", "T>C": "#a1cf64", "T>G": "#edc8c5"}
INK, MUTED = "#222222", "#555555"
TUMORS = {"HG008T": [("Illumina", "m6np__HG008T_Illumina__vcf_chr1"), ("PacBio HiFi", "m6np__HG008T_PacBio__vcf_chr1"),
                     ("ONT", "m6np__HG008T_ONT__vcf_chr1")],
          "COLO829T": [("Illumina", "m6np__COLO829T_Illumina__vcf_chr1"), ("Fiber-seq", "m6np__COLO829T_fiberseq__vcf_chr1"),
                       ("ONT", "m6np__COLO829T_ONT__vcf_chr1")],
          "HCC1395T": [("Illumina", "m6np__HCC1395T_Illumina__vcf_chr1"), ("PacBio HiFi", "m6np__HCC1395T_PacBio__vcf_chr1"),
                       ("ONT", "m6np__HCC1395T_ONT__vcf_chr1")]}
OLD = {"HG008T": [("Illumina", "rcpu__HG008_Illumina_SNV_d46_7fc7396__vcf_chr1"), ("PacBio HiFi", "rcpu__HG008_PacBio_SNV_d46_7fc7396__vcf_chr1"),
                  ("ONT", "rcpu__HG008_ONT_SNV_d46_7fc7396__vcf_chr1")],
       "COLO829T": [("Illumina", "rcpu__COLO829_Illumina_SNV_d46_7fc7396__vcf_chr1"), ("Fiber-seq", "rcpu__COLO829_fiberseq_SNV_d46_7fc7396__vcf_chr1"),
                    ("ONT", "rcpu__COLO829_ONT_SNV_d46_7fc7396__vcf_chr1")]}
csv.field_size_limit(sys.maxsize)


def cosine(a, b):
    dot = sum(a.get(k, 0) * b.get(k, 0) for k in CHANNELS)
    na, nb = math.sqrt(sum(a.get(k, 0) ** 2 for k in CHANNELS)), math.sqrt(sum(b.get(k, 0) ** 2 for k in CHANNELS))
    return dot / (na * nb) if na and nb else float("nan")


def uv(c):
    n = sum(c.values())
    return sum(v for k, v in c.items() if k[2:5] == "C>T" and k[0] in "CT") / n if n else float("nan")


def e18_spectra():
    """tumor -> {'truth': Counter, 'germline': Counter} from E18 (chr1, benchmark regions)."""
    out = {}
    with open(E18) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            kind, _, tumor = row["spectrum"].partition(" ")
            if kind in ("truth", "germline") and tumor in TUMORS and row["kind"] == kind:
                out.setdefault(tumor, {})[kind] = Counter({k: int(row[k]) for k in CHANNELS if int(row[k])})
    return out


def graph_coordinate(calls_tsv):
    """SBS96 of the unplaced PASS SNV calls of one run, context from E9c's graph-path REF context.
    Returns (Counter, calls without a usable context, calls read)."""
    c, skipped, total = Counter(), 0, 0
    with open(calls_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["group"] != "unplaced_pass" or row["kind_event"] != "SNP":
                continue
            total += 1
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
    return c, skipped, total


def panel(ax, c, title, pad=4):
    n = sum(c.values())
    ax.bar(range(96), [c.get(k, 0) / n * 100 if n else 0 for k in CHANNELS], color=[COLOR[k[2:5]] for k in CHANNELS],
           width=0.8)
    ax.set_title(title, fontsize=8.5, loc="left", pad=pad, color=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for i in range(1, 6):
        ax.axvline(16 * i - 0.5, color="#bbbbbb", lw=0.5)
    ax.tick_params(axis="y", labelsize=6.5)
    ax.set_xlim(-1.5, 96.5)
    ax.set_xticks([])


def class_header(ax):
    """Substitution class above each block of 16 channels: a strip in the class colour, the label in ink."""
    tr = ax.get_xaxis_transform()
    for i, c in enumerate(CLASSES):
        ax.plot([16 * i - 0.1, 16 * i + 15.1], [1.012, 1.012], transform=tr, color=COLOR[c], lw=3, solid_capstyle="butt",
                clip_on=False)
        ax.text(16 * i + 7.5, 1.03, c, transform=tr, ha="center", va="bottom", fontsize=7, fontweight="bold", color=INK)


def repro():
    ref = {(r["tumor"], r["platform"]): r for r in csv.DictReader(open(OLD_COS), delimiter="\t")
           if r["use"] == "Supplement graph-coordinate"}
    e18 = e18_spectra()
    ok = True
    for tumor, plats in OLD.items():
        for plat, run in plats:
            c, skipped, total = graph_coordinate(OLD_RUNS / run / "calls.tsv")
            got = (str(sum(c.values())), f"{cosine(c, e18[tumor]['truth']):.3f}", f"{cosine(c, e18[tumor]['germline']):.3f}", f"{uv(c):.3f}")
            r = ref[(tumor, plat)]
            want = (r["n"], r["cosine_to_truth"], r["cosine_to_normal_germline"], r["uv_type_ct_share"])
            ok &= got == want
            print(tumor, plat, "this script", got, "skipped", skipped, "| R3.5/cosine.tsv", want, "OK" if got == want else "DIFFERENT")
    print("reproduction", "OK" if ok else "FAILED")
    return 0 if ok else 1


def main():
    if "--repro" in sys.argv:
        return repro()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    e18 = e18_spectra()
    assert set(e18) == set(TUMORS) and all(set(v) == {"truth", "germline"} for v in e18.values()), e18.keys()
    tsv, cos, n_all, skipped_all = [], [], 0, 0
    fig, axes = plt.subplots(3, 3, figsize=(15, 9))
    for i, (tumor, plats) in enumerate(TUMORS.items()):
        truth, germ = e18[tumor]["truth"], e18[tumor]["germline"]
        for j, (plat, run) in enumerate(plats):
            c, skipped, total = graph_coordinate(E9C_RUNS / run / "calls.tsv")
            n = sum(c.values())
            assert n + skipped == total
            n_all += total
            skipped_all += skipped
            panel(axes[i, j], c, f"{tumor} {plat} (n = {n:,}; cosine similarity to chr1 truth {cosine(c, truth):.2f})",
                  pad=16 if i == 0 else 5)
            tsv.append([tumor, plat, run, n, skipped] + [c.get(k, 0) for k in CHANNELS])
            cos.append([tumor, plat, run, total, n, skipped, f"{cosine(c, truth):.3f}", f"{cosine(c, germ):.3f}", f"{uv(c):.3f}"])
            print(tumor, plat, "graph-coordinate PASS", total, "in the spectrum", n, "no usable context", skipped)
        axes[i, 0].set_ylabel("% of SNVs", fontsize=7.5, color=MUTED)
    for j in range(3):
        class_header(axes[0, j])
        axes[2, j].set_xticks(range(96))
        axes[2, j].set_xticklabels([k[0] + k[2] + k[6] for k in CHANNELS], rotation=90, fontsize=3.8, family="monospace")
    fig.suptitle("Pansoma graph-coordinate SNV PASS calls, chr1 (six-tumor models, HPRC v2.1, default output); "
                 "trinucleotide context from the graph path", fontsize=10, x=0.01, ha="left", color=INK)
    fig.tight_layout()
    out = HERE / "fig_3.3_sbs96_graph_coordinate.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    for tumor in TUMORS:
        t, g = e18[tumor]["truth"], e18[tumor]["germline"]
        cos.append([tumor, "-", "reference: chr1 somatic truth (benchmark regions)", sum(t.values()), sum(t.values()), 0, "1.000", f"{cosine(t, g):.3f}", f"{uv(t):.3f}"])
        cos.append([tumor, "-", "reference: matched-normal germline", sum(g.values()), sum(g.values()), 0, f"{cosine(g, t):.3f}", "1.000", f"{uv(g):.3f}"])
    with open(HERE / "fig_3.3_sbs96.tsv", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["tumor", "platform", "run", "n", "no_usable_context"] + CHANNELS)
        w.writerows(tsv)
    with open(HERE / "fig_3.3_cosine.tsv", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["tumor", "platform", "run", "graph_coordinate_pass", "n_in_spectrum", "no_usable_context", "cosine_to_truth",
                    "cosine_to_normal_germline", "uv_type_ct_share"])
        w.writerows(cos)
    print(f"total graph-coordinate PASS calls {n_all}, without a usable context {skipped_all} -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

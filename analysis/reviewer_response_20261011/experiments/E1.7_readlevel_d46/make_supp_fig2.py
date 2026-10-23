"""E1.7 (R2.3): the revised Supplementary Fig. 2, read-level classes of the truth alleles without a candidate on the
HPRC v2.1 (d46) tensor sets, and the table behind it.

Reads <SET>/summary.json (report.py) of the six sets; writes supp_fig2.png / .pdf and classes.tsv / classes.md here.
Classes of report.py, merged into five groups shared by SNVs and INDELs (colour = group in both panels):
  graph   ALT already in the graph            S3; I1 (fully) + I2 (partly, residual edit with another spelling)
  repeat  repeat: absorbed or written apart   S4a / S4b / S4c; I3 (reads on GRCh38, another spelling)
  alt     fewer than 3 ALT reads / low AF     S2; I4
  mapq    no read with MAPQ > 10              S1; I5
  other   other (indel > 50 bp, unclear)      I6, I7, I9 and anything else
usage: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_supp_fig2.py   (base anaconda matplotlib is built against NumPy 1.x)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
SETS = [("HG008T_Illumina_d46", "HG008T Illumina"), ("HG008T_PacBio_d46", "HG008T PacBio HiFi"),
        ("HG008T_ONT_d46", "HG008T ONT"), ("COLO829T_Illumina_d46", "COLO829T Illumina"),
        ("COLO829T_fiberseq_d46", "COLO829T Fiber-seq"), ("COLO829T_ONT_d46", "COLO829T ONT")]
GROUPS = [("graph", "ALT already in the graph", "#2a78d6"),
          ("repeat", "Repeat: absorbed or written differently", "#eb6834"),
          ("alt", "Fewer than 3 ALT reads or low allele fraction", "#1baf7a"),
          ("mapq", "No read with MAPQ > 10", "#eda100"),
          ("other", "Other (indel > 50 bp, unclear)", "#e87ba4")]
INK, MUTED, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"


def group(cls):
    code = cls.split("=", 1)[-1].split("_", 1)[0]
    if code in ("S3", "I1", "I2"):
        return "graph"
    if code.startswith("S4") or code == "I3":
        return "repeat"
    if code in ("S2", "I4"):
        return "alt"
    if code in ("S1", "I5"):
        return "mapq"
    return "other"


def counts(summary, kinds):
    out = {g: 0 for g, _, _ in GROUPS}
    for kind in kinds:
        for cls, v in summary["by_kind"].get(f"miss:{kind}", {}).get("groups", {}).items():
            if cls.startswith("final_class="):
                out[group(cls)] += v["n"]
    return out


def main():
    rows = []
    for key, label in SETS:
        s = json.loads((HERE / key / "summary.json").read_text())
        for kind, kinds in (("SNV", ("SNP",)), ("INDEL", ("INS", "DEL"))):
            c = counts(s, kinds)
            rows.append(dict(set=key, label=label, kind=kind, total=sum(c.values()), **c))
    cols = ["set", "kind", "total"] + [g for g, _, _ in GROUPS]
    (HERE / "classes.tsv").write_text("\t".join(cols) + "\n" + "".join(
        "\t".join(str(r[c]) for c in cols) + "\n" for r in rows))
    md = ["# Read-level classes of truth alleles without a candidate (HPRC v2.1, chr1-22; make_supp_fig2.py)", "",
          "| set | kind | no candidate | " + " | ".join(n for _, n, _ in GROUPS) + " |", "|---|---|---:|" + "---:|" * len(GROUPS)]
    for r in rows:
        md.append(f"| {r['label']} | {r['kind']} | {r['total']:,} | " + " | ".join(
            f"{r[g]:,} ({100 * r[g] / r['total']:.1f}%)" for g, _, _ in GROUPS) + " |")
    (HERE / "classes.md").write_text("\n".join(md) + "\n")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "text.color": INK,
                         "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK})
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True, facecolor=SURFACE)
    for ax, kind in zip(axes, ("SNV", "INDEL")):
        ax.set_facecolor(SURFACE)
        sub = [r for r in rows if r["kind"] == kind]
        for y, r in enumerate(sub):
            left = 0.0
            for g, _, colour in GROUPS:
                w = 100 * r[g] / r["total"] if r["total"] else 0
                if w > 0:
                    ax.barh(y, w, left=left, height=0.62, color=colour, edgecolor=SURFACE, linewidth=1.5)
                    if w >= 7:
                        ax.text(left + w / 2, y, f"{w:.0f}", ha="center", va="center", fontsize=8, color=INK)
                left += w
            ax.text(101, y, f"n = {r['total']:,}", va="center", fontsize=8, color=MUTED)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels([r["label"] for r in sub])
        ax.set_xlim(0, 100)
        ax.set_xlabel("% of truth alleles without a candidate")
        ax.set_title(f"{kind}s", loc="left", fontsize=10, color=INK)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].invert_yaxis()  # shared y axis: invert once (twice would undo it)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in GROUPS]
    fig.legend(handles, [n for _, n, _ in GROUPS], loc="lower center", ncol=3, frameon=False, fontsize=8.5)
    fig.tight_layout(rect=(0, 0.13, 0.97, 1), w_pad=4)
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"supp_fig2.{ext}", dpi=200, facecolor=SURFACE)
    print("\n".join(md))


if __name__ == "__main__":
    main()

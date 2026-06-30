#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def chr_sort_key(chrom):
    if chrom == "chrX":
        return 23
    if chrom == "chrY":
        return 24
    return int(chrom.replace("chr", ""))


def build_all_row(df_autosomes):
    rows = []

    for variant_type in ["SNV", "INDEL"]:
        sub = df_autosomes[df_autosomes["variant_type"] == variant_type].copy()

        rows.append({
            "chr": "all",
            "variant_type": variant_type,
            "COLO829_only_count": sub["COLO829_only_count"].sum(),
            "HPRC_only_count": sub["HPRC_only_count"].sum(),
            "overlap_count_COLO829rep": sub["overlap_count_COLO829rep"].sum(),
            "overlap_count_HPRCrep": sub["overlap_count_HPRCrep"].sum(),
            "COLO829_total_count": sub["COLO829_total_count"].sum(),
            "HPRC_total_count": sub["HPRC_total_count"].sum(),
        })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", required=True, help="Input TSV file")
    parser.add_argument("--out", default=None, help="Output figure path, e.g. plot.png or plot.pdf")
    parser.add_argument(
        "--ca",
        action="store_true",
        help="Only plot chr1 and all (sum of chr1-chr22); exclude chrX and chrY"
    )
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df["variant_type"] = df["variant_type"].str.upper()

    required_cols = [
        "chr",
        "variant_type",
        "COLO829_only_count",
        "HPRC_only_count",
        "overlap_count_COLO829rep",
        "overlap_count_HPRCrep",
        "COLO829_total_count",
        "HPRC_total_count",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if args.ca:
        autosomes = [f"chr{i}" for i in range(1, 23)]
        df_auto = df[df["chr"].isin(autosomes)].copy()

        df_all = build_all_row(df_auto)

        df_chr1 = df_auto[df_auto["chr"] == "chr1"].copy()

        df = pd.concat([df_chr1, df_all], ignore_index=True)

        chrom_order = ["chr1", "all"]
    else:
        chrom_order = sorted(df["chr"].unique(), key=chr_sort_key)

    df_snv = df[df["variant_type"] == "SNV"].copy()
    df_indel = df[df["variant_type"] == "INDEL"].copy()

    df_snv = df_snv.set_index("chr").loc[chrom_order].reset_index()
    df_indel = df_indel.set_index("chr").loc[chrom_order].reset_index()

    df_full = pd.DataFrame({
        "Chromosome": chrom_order,
        "SNVs_Overlap": df_snv["overlap_count_COLO829rep"].values,
        "SNVs_A_unique": df_snv["HPRC_only_count"].values,
        "SNVs_B_unique": df_snv["COLO829_only_count"].values,
        "INDELs_Overlap": df_indel["overlap_count_COLO829rep"].values,
        "INDELs_A_unique": df_indel["HPRC_only_count"].values,
        "INDELs_B_unique": df_indel["COLO829_only_count"].values,
    })

    x = np.arange(len(df_full))
    bar_width = 0.38

    color_overlap = "#4daf4a"
    color_A = "#e6ab02"
    color_B = "#377eb8"

    plt.rcParams.update({
        "font.size": 18,
        "axes.titlesize": 22,
        "axes.labelsize": 20,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 18
    })

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(24, 14), sharex=True)

    # SNVs
    ax1.bar(x - bar_width/2, df_full["SNVs_Overlap"], width=bar_width, color=color_overlap, label="Overlap")
    ax1.bar(
        x - bar_width/2,
        df_full["SNVs_A_unique"],
        width=bar_width,
        bottom=df_full["SNVs_Overlap"],
        color=color_A,
        label="HPRC Human Pangenome Graph AF-Filtered Full Graph"
    )
    ax1.bar(x + bar_width/2, df_full["SNVs_Overlap"], width=bar_width, color=color_overlap)
    ax1.bar(
        x + bar_width/2,
        df_full["SNVs_B_unique"],
        width=bar_width,
        bottom=df_full["SNVs_Overlap"],
        color=color_B,
        label="COLO829BL Germline Variants"
    )
    ax1.set_title("SNVs per Chromosome: COLO829BL Germline variants overlapping the HPRC Human Pangenome Graph")
    ax1.set_ylabel("Variant Count")
    ax1.legend(ncol=3, loc="upper right")

    # INDELs
    ax2.bar(x - bar_width/2, df_full["INDELs_Overlap"], width=bar_width, color=color_overlap, label="Overlap")
    ax2.bar(
        x - bar_width/2,
        df_full["INDELs_A_unique"],
        width=bar_width,
        bottom=df_full["INDELs_Overlap"],
        color=color_A,
        label="HPRC Human Pangenome Graph AF-Filtered Full Graph"
    )
    ax2.bar(x + bar_width/2, df_full["INDELs_Overlap"], width=bar_width, color=color_overlap)
    ax2.bar(
        x + bar_width/2,
        df_full["INDELs_B_unique"],
        width=bar_width,
        bottom=df_full["INDELs_Overlap"],
        color=color_B,
        label="COLO829BL Germline Variants"
    )
    ax2.set_title("INDELs per Chromosome: COLO829BL Germline variants overlapping the HPRC Human Pangenome Graph")
    ax2.set_ylabel("Variant Count")
    ax2.set_xticks(x)
    ax2.set_xticklabels(df_full["Chromosome"], rotation=45)
    ax2.legend(ncol=3, loc="upper right")

    fig.text(0.5, 0.04, "Chromosome", ha="center", fontsize=20)

    plt.tight_layout(rect=[0, 0.05, 1, 1])

    if args.out:
        plt.savefig(args.out, dpi=300, bbox_inches="tight")
        print(f"[INFO] Saved figure to: {args.out}")

    plt.show()


if __name__ == "__main__":
    main()
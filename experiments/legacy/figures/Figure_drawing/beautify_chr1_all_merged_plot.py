
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch

hg008_path = "reads_on_hg38_other_paths_summary_HG008-T_p23_BCM_IlluminaWGS_20240313.tsv"
colo_path = "reads_on_hg38_other_paths_summary_WGS_COLO829T.tsv"

def summarize_chr1_and_all(df, sample_name):
    chr1 = df[df["chr"].astype(str) == "chr1"]
    if chr1.empty:
        raise ValueError(f"{sample_name}: chr1 not found in TSV")
    chr1 = chr1.iloc[0]
    all_sum = df.select_dtypes(include=["number"]).sum()

    rows = []
    for label, src in [("chr1", chr1), ("all", all_sum)]:
        rows.append({
            "group": f"{sample_name}\n{label}",
            "GRCh38_total": src["GRCh38_path_total_nodes"],
            "Other_total": src["Other_paths_total_nodes"],
            "GRCh38_perfect": src["GRCh38_path_nodes_perfect_gt0"],
            "Other_perfect": src["Other_paths_nodes_perfect_gt0"],
            "GRCh38_notperfect": src["GRCh38_path_nodes_not_perfect_gt0"],
            "Other_notperfect": src["Other_paths_nodes_not_perfect_gt0"],
        })
    return pd.DataFrame(rows)

def compact_num(x, pos=None):
    if x >= 1e6:
        return f"{x/1e6:.0f}M"
    if x >= 1e3:
        return f"{x/1e3:.0f}K"
    return f"{int(x)}"

def compact_label(x):
    if x >= 1e6:
        return f"{x/1e6:.2f}M"
    if x >= 1e3:
        return f"{x/1e3:.1f}K"
    return str(int(x))

hg = pd.read_csv(hg008_path, sep="\t")
co = pd.read_csv(colo_path, sep="\t")

plot_df = pd.concat([
    summarize_chr1_and_all(hg, "HG008T"),
    summarize_chr1_and_all(co, "COLO829T"),
], ignore_index=True)

x = list(range(len(plot_df)))
bar_width = 0.22
pos_total = [p - bar_width for p in x]
pos_perfect = x
pos_notperfect = [p + bar_width for p in x]

g_total = plot_df["GRCh38_total"].tolist()
o_total = plot_df["Other_total"].tolist()
g_perfect = plot_df["GRCh38_perfect"].tolist()
o_perfect = plot_df["Other_perfect"].tolist()
g_not = plot_df["GRCh38_notperfect"].tolist()
o_not = plot_df["Other_notperfect"].tolist()

fig, ax = plt.subplots(figsize=(15, 8), dpi=220)

default_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
c_total = default_colors[0]
c_perfect = default_colors[2]
c_not = default_colors[1]

ax.bar(pos_total, g_total, width=bar_width, color=c_total)
ax.bar(pos_total, o_total, width=bar_width, bottom=g_total, color=c_total, alpha=0.42, hatch="///", linewidth=0)

ax.bar(pos_perfect, g_perfect, width=bar_width, color=c_perfect)
ax.bar(pos_perfect, o_perfect, width=bar_width, bottom=g_perfect, color=c_perfect, alpha=0.42, hatch="///", linewidth=0)

ax.bar(pos_notperfect, g_not, width=bar_width, color=c_not)
ax.bar(pos_notperfect, o_not, width=bar_width, bottom=g_not, color=c_not, alpha=0.42, hatch="///", linewidth=0)

ax.set_yscale("log")
ax.set_ylabel("Node count (log scale)", fontsize=12)
ax.set_xlabel("Sample and scope", fontsize=12)
ax.set_title("GRCh38 vs Other paths node counts", fontsize=17, pad=22)

ax.set_xticks(x)
ax.set_xticklabels(plot_df["group"], fontsize=11)
ax.yaxis.set_major_formatter(FuncFormatter(compact_num))

ax.grid(axis="y", which="major", alpha=0.18, linewidth=0.8)
ax.grid(axis="y", which="minor", alpha=0.06, linewidth=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.axvline(1.5, color="0.82", linewidth=1.0, linestyle="--", zorder=0)

ymax = max((plot_df["GRCh38_total"] + plot_df["Other_total"]).max(),
           (plot_df["GRCh38_perfect"] + plot_df["Other_perfect"]).max(),
           (plot_df["GRCh38_notperfect"] + plot_df["Other_notperfect"]).max())
ax.text(0.5, ymax * 1.12, "HG008T", ha="center", va="bottom", fontsize=12)
ax.text(2.5, ymax * 1.12, "COLO829T", ha="center", va="bottom", fontsize=12)

stack_tops = [
    [gt + ot for gt, ot in zip(g_total, o_total)],
    [gp + op for gp, op in zip(g_perfect, o_perfect)],
    [gn + on for gn, on in zip(g_not, o_not)],
]
positions = [pos_total, pos_perfect, pos_notperfect]

for pos_list, top_vals in zip(positions, stack_tops):
    for xpos, y in zip(pos_list, top_vals):
        ax.text(xpos, y * 1.05, compact_label(y), ha="center", va="bottom", fontsize=8)

legend1 = ax.legend(
    handles=[
        Patch(facecolor=c_total, label="Total"),
        Patch(facecolor=c_perfect, label="Perfect > 0"),
        Patch(facecolor=c_not, label="Not perfect > 0"),
    ],
    loc="upper left",
    frameon=False,
    ncol=3,
    bbox_to_anchor=(0.0, 1.03),
    borderaxespad=0.0
)
ax.add_artist(legend1)

ax.legend(
    handles=[
        Patch(facecolor="0.35", label="GRCh38"),
        Patch(facecolor="0.35", alpha=0.42, hatch="///", label="Other paths"),
    ],
    loc="upper right",
    frameon=False,
    ncol=2,
    bbox_to_anchor=(1.0, 1.03),
    borderaxespad=0.0
)

plt.subplots_adjust(top=0.80, bottom=0.16, left=0.08, right=0.98)
fig.savefig("HG008T_COLO829T_chr1_all_merged_beautified_v2.png", bbox_inches="tight")
fig.savefig("HG008T_COLO829T_chr1_all_merged_beautified_v2.pdf", bbox_inches="tight")
plt.close(fig)

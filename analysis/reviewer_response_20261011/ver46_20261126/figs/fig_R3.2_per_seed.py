#!/usr/bin/env python
"""Bar plot of Response Table 3.2, letter ver4.6: per-seed chr1 results after population filtering.

Rows F1 / Precision / Recall, columns = tumors, three platforms on x, three bars per platform: the final INDEL
models (six-tumor, four-tumor) and the six-tumor SNV models on HPRC v2.1. Bar = mean of the three training seeds
(the value printed in the table), points = seeds 0 / 1 / 2. Layout of fig_R3.2_per_seed_ver4.5.png without its
fourth column (single-tumor COLO829T SNV models, no longer in the table); the SNV bars are now in every tile.

Data: letter_tables.json ("Response Table 3.2" of ../build/Reviewer_Letter_Response_ver4.6.docx, written by
      extract_letter_tables.py): cells "s0 / s1 / s2 (mean ± sd)"; the nine SNV rows are asserted equal to
      ../tables/tables46.json t32_snv.rows.
Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python fig_R3.2_per_seed.py
"""
import json
import re

import figstyle as S
from figstyle import plt

SERIES = [("6 tumors, INDEL", "INDEL, six-tumor", S.C_V21), ("4 tumors, INDEL", "INDEL, four-tumor", S.C_FOUR),
          ("6 tumors, SNV", "SNV, six-tumor (HPRC v2.1)", S.C_SNV6)]
METRICS = [(5, "F1"), (3, "Precision"), (4, "Recall")]          # column of the table
CELL = re.compile(r"^(\d\.\d{3}) / (\d\.\d{3}) / (\d\.\d{3}) \((\d\.\d{3}) ± (\d\.\d{3})\)$")

rows = json.loads((S.HERE / "letter_tables.json").read_text())["Response Table 3.2"]
assert rows[0] == ["tumor", "platform", "training set", "precision", "recall", "F1"] and len(rows) == 28
snv = json.loads((S.V46 / "tables/tables46.json").read_text())["t32_snv"]["rows"]
assert [r for r in rows[1:] if r[2] == "6 tumors, SNV"] == snv
val = {}
for r in rows[1:]:
    for col, _ in METRICS:
        m = CELL.match(r[col])
        assert m, r
        s = [float(x) for x in m.groups()]
        assert abs(sum(s[:3]) / 3 - s[3]) < 0.0011, r          # printed mean = mean of the printed seeds (rounding)
        val[r[0], r[1], r[2], col] = (s[:3], s[3])
assert len(val) == 27 * 3

fig = plt.figure(figsize=(6.5, 5.2))
left, right, bottom, top, gx, gy = 0.058, 0.994, 0.062, 0.925, 0.017, 0.018
wd = (right - left - 2 * gx) / 3
ht = (top - bottom - 2 * gy) / 3
for r, (col, ylab) in enumerate(METRICS):
    for i, tumor in enumerate(S.TUMORS):
        ax = fig.add_axes([left + i * (wd + gx), top - (r + 1) * ht - r * gy, wd, ht])
        S.classic(ax, (0, 1.14), [0, 0.2, 0.4, 0.6, 0.8, 1.0], left=(i == 0), bottom=(r == 2))
        if r < 2:
            ax.axhline(0, color=S.GRID, lw=S.GRID_LW, zorder=0)
        plats = S.platforms(tumor)
        bw = 0.84 / 3
        for j, p in enumerate(plats):
            for k, (key, _, colr) in enumerate(SERIES):
                seeds, mean = val[tumor, p, key, col]
                x = j - 0.42 + bw * (k + 0.5)
                ax.bar(x, mean, width=bw, color=colr, edgecolor="white", linewidth=0.5, zorder=2)
                ax.scatter([x] * 3, seeds, s=2.6, color=S.INK, linewidths=0, zorder=4)
                ax.text(x, max(seeds + [mean]) + 0.022, f"{mean:.3f}", rotation=90, ha="center", va="bottom",
                        fontsize=S.PT - 1.6, color=S.INK, zorder=3)
        ax.set_xlim(-0.55, 2.55)
        ax.set_xticks(range(3))
        ax.set_xticklabels(plats if r == 2 else [], fontsize=S.PT - 0.5, color=S.INK)
        if r == 0:
            ax.set_title(tumor, fontsize=S.PT + 0.5, fontweight="bold", color=S.INK, pad=3)
        if i == 0:
            ax.set_ylabel(ylab, fontsize=S.PT, color=S.INK, labelpad=2)
fig.text(0.006, 0.012, "Bars, mean of three training seeds; points, the individual seeds (chr1, after population filtering).",
         ha="left", va="center", fontsize=S.PT - 0.5, color=S.MUTED)
S.top_legend(fig, [(c, n) for _, n, c in SERIES], x=left, y=0.99)
S.save(fig, "fig_R3.2_per_seed", "fig_R3.2_per_seed", "fig_R3.2_per_seed")

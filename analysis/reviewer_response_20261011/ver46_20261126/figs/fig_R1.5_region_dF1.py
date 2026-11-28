#!/usr/bin/env python
"""Bar plot of the SNV columns of Response Table 1.5 (= 2.5 = 3.15), letter ver4.6: SNV F1 gain of the HPRC graphs
over the GRCh38-only graph by region, chr1, before population filtering; six-tumor models without site features.

Rows = strata (SMaHT extreme; GIAB all difficult, low mappability or segdup, tandem repeats and homopolymers, outside
difficult), columns = tumors, three platforms on x, two bars per platform (HPRC v1.1, HPRC v2.1) with the paired 95 %
block-bootstrap interval. Layout of the ver4.5 figure fig_R1.5_region_dF1_ver4.5.png.

Data: ../regions/cells_six_tumor.json (region_nositef.py), cells["before|<dataset>|<d9|d46>"].values[<stratum>]
      = {d, lo, hi}; the printed label is asserted equal to the table cell text of the same file.
Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python fig_R1.5_region_dF1.py
"""
import json

import figstyle as S
from figstyle import plt

STRATA = [("extreme", "SMaHT\nextreme"), ("alldifficultregions", "GIAB\nall difficult"),
          ("alllowmapandsegdupregions", "GIAB low map\nor segdup"),
          ("AllTandemRepeatsandHomopolymers_slop5", "GIAB tandem\nrepeats and\nhomopolymers"),
          ("notinalldifficultregions", "GIAB outside\ndifficult")]
GRAPHS = [("d9", "Pansoma, HPRC v1.1", S.C_V11), ("d46", "Pansoma, HPRC v2.1", S.C_V21)]
KEY = {"Illumina": "Illumina", "PacBio HiFi": "PacBio", "Fiber-seq": "fiberseq", "ONT": "ONT"}

cells = json.loads((S.V46 / "regions/cells_six_tumor.json").read_text())["cells"]
assert sum(k.startswith("before|") for k in cells) == 18

fig = plt.figure(figsize=(6.51, 6.2))
left, right, bottom, top, gx, gy = 0.050, 0.885, 0.030, 0.940, 0.017, 0.012
wd = (right - left - 2 * gx) / 3
ht = (top - bottom - 4 * gy) / 5
allhi = []
for r, (st, stname) in enumerate(STRATA):
    for i, tumor in enumerate(S.TUMORS):
        ax = fig.add_axes([left + i * (wd + gx), top - (r + 1) * ht - r * gy, wd, ht])
        S.classic(ax, (0, 1.0), [0, 0.2, 0.4, 0.6, 0.8], left=(i == 0), bottom=(r == 4))
        if r < 4:                                   # the old figure has a faint baseline in the upper rows
            ax.axhline(0, color=S.GRID, lw=S.GRID_LW, zorder=0)
        plats = S.platforms(tumor)
        bw = 0.8 / 2
        for j, p in enumerate(plats):
            for k, (g, _, col) in enumerate(GRAPHS):
                c = cells[f"before|{tumor}_{KEY[p]}|{g}"]
                v = c["values"][st]
                d, lo, hi = v["d"], v["lo"], v["hi"]
                assert c["text"][st] == f"{d:+.3f} [{lo:.3f}, {hi:.3f}]", (c["text"][st], d, lo, hi)
                assert d > 0 and lo > 0
                x = j - 0.4 + bw * (k + 0.5)
                ax.bar(x, d, width=bw, color=col, edgecolor="white", linewidth=0.5, zorder=2)
                ax.errorbar(x, d, yerr=[[d - lo], [hi - d]], fmt="none", ecolor=S.INK, elinewidth=0.6, capsize=2.2, capthick=0.6, zorder=3)
                ax.text(x, hi + 0.018, f"{d:.3f}", ha="center", va="bottom", fontsize=S.PT - 1.6, color=S.INK, zorder=4)
                allhi.append(hi)
        ax.set_xlim(-0.55, 2.55)
        ax.set_xticks(range(3))
        ax.set_xticklabels(plats if r == 4 else [], fontsize=S.PT - 0.5, color=S.INK)
        if r == 0:
            ax.set_title(tumor, fontsize=S.PT + 0.5, fontweight="bold", color=S.INK, pad=3)
        if i == 2:
            ax.text(1.025, 0.5, stname, transform=ax.transAxes, ha="left", va="center", fontsize=S.PT, color=S.INK, linespacing=1.05)
assert max(allhi) < 0.93, max(allhi)                # labels stay inside the panel
fig.text(0.012, (top + bottom) / 2, "SNV F1 gain over GRCh38-only (chr1, before filtering; bars, 95% bootstrap interval)",
         rotation=90, ha="center", va="center", fontsize=S.PT, color=S.INK)
S.top_legend(fig, [(c, n) for _, n, c in GRAPHS], x=left, y=0.992)
S.save(fig, "fig_R1.5_region_dF1", "fig_R1.5_region_dF1", "fig_R1.5_region_dF1")

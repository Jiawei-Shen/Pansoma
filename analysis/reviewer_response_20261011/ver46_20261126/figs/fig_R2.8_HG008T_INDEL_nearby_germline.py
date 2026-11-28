#!/usr/bin/env python
"""Bar plot of Response Table 2.8, letter ver4.6: INDEL calling on HG008T chr1 by the number of matched-normal
germline variants within 50 bp (0 / 1 / >=2), after population filtering.

Rows F1 / Precision / Recall / False positives, columns = platforms, x = the three groups, THREE bars per group:
Pansoma on GRCh38-only, HPRC v1.1 and HPRC v2.1 (six-tumor INDEL models, seed 0, none with site features). The
ver4.5 figure had two HPRC v2.1 bars (best seed, mean of three seeds); the table and the figure now have one.

Data: ../strata/strata_indel_nositef.json (strata_indel_nositef.py), rows[*].groups[<0|1|2+>] = {p, r, f, fp}.
Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python fig_R2.8_HG008T_INDEL_nearby_germline.py
"""
import json

import figstyle as S
from figstyle import plt

REFS = [("GRCh38-only", "Pansoma, GRCh38-only", S.C_LINEAR), ("HPRC v1.1", "Pansoma, HPRC v1.1", S.C_V11),
        ("HPRC v2.1", "Pansoma, HPRC v2.1", S.C_V21)]
PLATS = ("Illumina", "PacBio HiFi", "ONT")
GROUPS = (("0", "0"), ("1", "1"), ("2+", "≥2"))
# metric key, y label, y limit, ticks, tick format, value format
ROWS = [("f", "F1", (0, 0.42), [0, 0.1, 0.2, 0.3, 0.4], "{:.1f}", "{:.3f}"),
        ("p", "Precision", (0, 1.22), [0, 0.2, 0.4, 0.6, 0.8, 1.0], "{:.1f}", "{:.3f}"),
        ("r", "Recall", (0, 0.335), [0, 0.1, 0.2, 0.3], "{:.1f}", "{:.3f}"),
        ("fp", "False positives", (0, 430), [0, 100, 200, 300, 400], "{:.0f}", "{:.0f}")]

J = json.loads((S.V46 / "strata/strata_indel_nositef.json").read_text())
assert J["truth"] == {"0": 452, "1": 84, "2+": 51}
data = {(r["platform"], r["reference"]): r["groups"] for r in J["rows"]}
assert len(data) == 9 and set(data) == {(p, ref) for p in PLATS for ref, _, _ in REFS}

fig = plt.figure(figsize=(6.5, 6.4))
left, right, bottom, top, gx, gy = 0.062, 0.994, 0.052, 0.945, 0.017, 0.014
wd = (right - left - 2 * gx) / 3
ht = (top - bottom - 3 * gy) / 4
for r, (key, ylab, ylim, ticks, tf, vf) in enumerate(ROWS):
    for i, p in enumerate(PLATS):
        ax = fig.add_axes([left + i * (wd + gx), top - (r + 1) * ht - r * gy, wd, ht])
        S.classic(ax, ylim, ticks, left=(i == 0), fmt=tf, bottom=(r == 3))
        if r < 3:
            ax.axhline(0, color=S.GRID, lw=S.GRID_LW, zorder=0)
        bw = 0.63 / 3
        for j, (g, _) in enumerate(GROUPS):
            for k, (ref, _, col) in enumerate(REFS):
                v = data[p, ref][g][key]
                assert v <= ticks[-1] * 1.0001 or key == "fp", (p, ref, g, key, v)
                x = j - 0.315 + bw * (k + 0.5)
                ax.bar(x, v, width=bw, color=col, edgecolor="white", linewidth=0.5, zorder=2)
                ax.text(x, v + (ylim[1] - ylim[0]) * 0.014, vf.format(v), rotation=90, ha="center", va="bottom",
                        fontsize=S.PT - 1.6, color=S.INK, zorder=3)
        ax.set_xlim(-0.6, 2.6)
        ax.set_xticks(range(3))
        ax.set_xticklabels([n for _, n in GROUPS] if r == 3 else [], fontsize=S.PT - 0.5, color=S.INK)
        if r == 0:
            ax.set_title(p, fontsize=S.PT + 0.5, fontweight="bold", color=S.INK, pad=3)
        if i == 0:
            ax.set_ylabel(ylab, fontsize=S.PT, color=S.INK, labelpad=2)
fig.align_ylabels()
fig.text((left + right) / 2, 0.012, "Germline variants near the site (HG008T chr1 INDELs)", ha="center", va="center",
         fontsize=S.PT, color=S.INK)
S.top_legend(fig, [(c, n) for _, n, c in REFS], x=left, y=0.992)
S.save(fig, "fig_R2.8_HG008T_INDEL_nearby_germline", "fig_R2.8_HG008T_INDEL_nearby_germline", "fig_R2.8_HG008T_INDEL_nearby_germline")

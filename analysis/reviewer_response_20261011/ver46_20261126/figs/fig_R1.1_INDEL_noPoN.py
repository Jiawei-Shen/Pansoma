#!/usr/bin/env python
"""Bar plot of Response Table 1.2 (= 2.10 = 3.4), letter ver4.6: chr1 INDELs with no population filter.

One tile per test tumor; three platforms on x; five bars per platform: Pansoma on GRCh38-only / HPRC v1.1 / HPRC v2.1
(six-tumor INDEL models, seed 0, none with site features), DeepSomatic-TO and ClairS-TO. Layout of the ver4.5 figure
fig_R1.1_INDEL_noPoN_<metric>.png (ggplot script fig_R1.1_SNV_noPoN.R).

Data: Pansoma  ../tables/tables46.json, t12.cells[<dataset>].before.{linear,d9,d46}
      comparators (unchanged)  Response Table 1.2 of ../build/Reviewer_Letter_Response_ver4.6.docx (ClairS-TO and
      DeepSomatic-TO columns; letter_tables.json, written by extract_letter_tables.py); the Pansoma cells of that
      table are asserted equal to the JSON at three decimals.
Output: fig_R1.1_INDEL_noPoN_{f1,precision,recall}.png / .pdf here and ..._ver4.6.* in Pansoma_Resubmit/fig_R1.1_SNV_noPoN.

Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python fig_R1.1_INDEL_noPoN.py
"""
import json

import figstyle as S
from figstyle import plt

MODELS = [("linear", "Pansoma, GRCh38-only", S.C_LINEAR), ("d9", "Pansoma, HPRC v1.1", S.C_V11), ("d46", "Pansoma, HPRC v2.1", S.C_V21),
          ("ds", "DeepSomatic-TO", S.C_DS), ("clairs", "ClairS-TO", S.C_CLAIRS)]
METRICS = {"f1": ("f1", "F1 score"), "precision": ("precision", "Precision"), "recall": ("recall", "Recall")}

cells = json.loads((S.V46 / "tables/tables46.json").read_text())["t12"]["cells"]
assert len(cells) == 9
rows = json.loads((S.HERE / "letter_tables.json").read_text())["Response Table 1.2"]     # extract_letter_tables.py
assert rows[0][1:] == ["Pansoma GRCh38-only", "Pansoma HPRC v1.1", "Pansoma HPRC v2.1", "ClairS-TO", "DeepSomatic-TO"], rows[0]
assert len(rows) == 10
val = {}
for r in rows[1:]:
    ds = r[0]
    b = cells[ds]["before"]
    trip = lambda s: dict(zip(("precision", "recall", "f1"), (float(x) for x in s.split(" / "))))  # noqa: E731
    for arm, txt in zip(("linear", "d9", "d46"), r[1:4]):      # the letter's Pansoma cells = the JSON
        assert txt == " / ".join(f"{b[arm][k]:.3f}" for k in ("precision", "recall", "f1")), (ds, arm, txt)
    val[ds] = {"linear": b["linear"], "d9": b["d9"], "d46": b["d46"], "clairs": trip(r[4]), "ds": trip(r[5])}

for metric, (key, ylab) in METRICS.items():
    fig = plt.figure(figsize=(6.5, 2.8))
    left, right, bottom, top, gap = 0.052, 0.994, 0.072, 0.865, 0.017
    wd = (right - left - 2 * gap) / 3
    for i, tumor in enumerate(S.TUMORS):
        ax = fig.add_axes([left + i * (wd + gap), bottom, wd, top - bottom])
        S.classic(ax, (0, 1.14), [0, 0.2, 0.4, 0.6, 0.8, 1.0], left=(i == 0))
        plats = S.platforms(tumor)
        bw = 0.86 / len(MODELS)
        for j, p in enumerate(plats):
            for k, (arm, _, col) in enumerate(MODELS):
                v = val[f"{tumor} {p}"][arm][key]
                x = j - 0.43 + bw * (k + 0.5)
                ax.bar(x, v, width=bw, color=col, edgecolor="white", linewidth=0.5, zorder=2)
                ax.text(x, v + 0.012, f"{v:.3f}", rotation=90, ha="center", va="bottom", fontsize=S.PT - 1.6, color=S.INK, zorder=3)
        ax.set_xlim(-0.55, 2.55)
        ax.set_xticks(range(3))
        ax.set_xticklabels(plats, fontsize=S.PT - 0.5, color=S.INK)
        ax.set_title(tumor, fontsize=S.PT + 0.5, fontweight="bold", color=S.INK, pad=3)
        if i == 0:
            ax.set_ylabel(f"{ylab} (chr1 INDELs, no population filter)", fontsize=S.PT, color=S.INK, labelpad=2)
    S.top_legend(fig, [(c, n) for _, n, c in MODELS], x=left, y=0.985)
    S.save(fig, f"fig_R1.1_INDEL_noPoN_{metric}", "fig_R1.1_SNV_noPoN", f"fig_R1.1_INDEL_noPoN_{metric}")

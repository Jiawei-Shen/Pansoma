#!/usr/bin/env python
"""Bar plot of Response Table 3.9, letter ver4.6: INDEL false positives before population filtering on the
GRCh38-only graph by germline / graph class, and the total on HPRC v2.1 (chr1, benchmark regions).

Layout of fig_R3.9_INDEL_FP_classes_ver4.5.png (= fig_R2.6_SNV_FP_classes_ver4.5.png); its geometry was measured
on that image (3900 x 2148 px at 600 dpi) and is kept in pixels below.
  a  class of the GRCh38-graph false positives (matched-normal germline in / not in the HPRC v2.1 graph, non-germline)
  b  share of the germline alleles of each class that still form a candidate (tensor) on HPRC v2.1; a class with
     fewer than 20 alleles is labelled with its count ("8 of 12")
  c  total false positives before filtering, GRCh38-only and HPRC v2.1 (six-tumor INDEL models without site
     features on both references, seed 0), log scale

Data: ../fpdecomp/fp_decomp.json (fp_decomp_nositef.py; kind INDEL, stage raw) for the bar lengths and
      ../fpdecomp/table_INDEL_before_filter.md (the letter's table) for the printed labels; the two are asserted
      consistent.
Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python fig_R3.9_INDEL_FP_classes.py
"""
import json
import re

import figstyle as S
from figstyle import plt

W, H = 3900, 2148                                # px of the ver4.5 image
ROWS_PX = [485, 635, 785, 1008, 1158, 1310, 1532, 1682, 1832]
TOP_PX, AXIS_PX = 395, 1922
BAND = 90                                        # grid lines run this far above / below the first / last row of a group
PANELS = {"a": (750, 1880), "b": (1980, 2760), "c": (2845, 3810)}

# ---- data
js = [x for x in json.loads((S.V46 / "fpdecomp/fp_decomp.json").read_text()) if x["kind"] == "INDEL"]
assert [x["dataset"] for x in js] == [f"{t} {p}" for t in S.TUMORS for p in S.platforms(t)]
md = [[c.strip() for c in ln.strip().strip("|").split("|")]
      for ln in (S.V46 / "fpdecomp/table_INDEL_before_filter.md").read_text().splitlines() if ln.startswith("|")]
md = [r for r in md if not set("".join(r)) <= set("-: ")]
assert len(md) == 10 and md[0][0] == "Dataset" and md[0][5].startswith("HPRC v2.1 FP before")


def pct(s):
    """'86%' -> 86.0, '<0.1%' -> None (only an upper bound is printed)."""
    return None if s.startswith("<") else float(s.rstrip("%"))


rows = []
for x, m in zip(js, md[1:]):
    assert m[0] == x["dataset"]
    r = x["stages"]["raw"]
    n = r["lin_fp_keys"]
    gi, go, ng = r["germline_in_graph"], r["germline_not_in_graph"], r["non_germline"]
    assert gi["n"] + go["n"] + ng["n"] == n == int(m[1].replace(",", "")) and r["d46_fp_records"] == int(m[5].replace(",", ""))
    s_in, k_in = m[2].split("; ")
    s_out, k_out = m[3].split("; ")
    d = {"dataset": x["dataset"], "tumor": x["tumor"], "n": n, "d46": r["d46_fp_records"],
         "share": [100 * gi["n"] / n, 100 * go["n"] / n, 100 * ng["n"] / n], "share_lab": [s_in, s_out, m[4]],
         "keep": [100 * gi["d46_tensor"] / gi["n"], 100 * go["d46_tensor"] / go["n"]], "keep_lab": [k_in, k_out],
         "n_lab": m[1], "d46_lab": m[5]}
    for v, lab in zip(d["share"], d["share_lab"]):                          # printed labels = the JSON numbers
        assert pct(lab) is not None and abs(v - pct(lab)) <= 0.5 + 1e-9, (x["dataset"], v, lab)
    for v, lab, c in zip(d["keep"], d["keep_lab"], (gi, go)):
        mm = re.match(r"^(\d+) of (\d+)$", lab)
        if mm:
            assert (c["d46_tensor"], c["n"]) == (int(mm.group(1)), int(mm.group(2))) and c["n"] < 20
        elif lab.startswith("<"):
            assert v < float(lab[1:].rstrip("%"))
        else:
            assert abs(v - pct(lab)) <= 0.5 + 1e-9, (x["dataset"], v, lab)
    rows.append(d)
assert len(rows) == 9

# ---- figure
fig = plt.figure(figsize=(W / S.DPI, H / S.DPI))
ov = fig.add_axes([0, 0, 1, 1])
ov.set_xlim(0, W)
ov.set_ylim(H, 0)
ov.axis("off")
axes = {}
for k, (x0, x1) in PANELS.items():
    ax = fig.add_axes([x0 / W, 1 - AXIS_PX / H, (x1 - x0) / W, (AXIS_PX - TOP_PX) / H])
    ax.set_ylim(AXIS_PX, TOP_PX)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(S.MUTED)
    ax.spines["bottom"].set_linewidth(S.LW)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=1.5, width=S.LW, color=S.MUTED, pad=1.5, labelsize=S.PT - 0.5, labelcolor=S.INK)
    ax.set_facecolor("none")
    axes[k] = ax
a, b, c = axes["a"], axes["b"], axes["c"]
a.set_xlim(0, 100 * (PANELS["a"][1] - PANELS["a"][0]) / 1125)
b.set_xlim(0, 100 * (PANELS["b"][1] - PANELS["b"][0]) / 660)
c.set_xscale("log")
c.set_xlim(10, 10 ** (1 + (PANELS["c"][1] - PANELS["c"][0]) / 192.5))
a.set_xticks([0, 25, 50, 75, 100])
b.set_xticks([0, 25, 50, 75, 100])
c.set_xticks([10, 100, 1000, 10000, 100000])
c.set_xticklabels(["10", "100", "1,000", "10,000", "100,000"])
c.minorticks_off()
for ax in (a, b, c):
    for g in range(3):
        y0, y1 = ROWS_PX[3 * g] - BAND, ROWS_PX[3 * g + 2] + BAND
        ax.vlines(ax.get_xticks(), y0, y1, color=S.GRID, lw=0.75, zorder=0)
a.set_xlabel("Share of GRCh38-graph FP (%)", fontsize=S.PT, color=S.INK, labelpad=2)
b.set_xlabel("Candidates retained (%)", fontsize=S.PT, color=S.INK, labelpad=2)
c.set_xlabel("FP before filtering (log scale)", fontsize=S.PT, color=S.INK, labelpad=2)

FS = S.PT - 1.0
COLS = (S.C_GERM_IN, S.C_GERM_OUT, S.C_NONGERM)
for d, y in zip(rows, ROWS_PX):
    plat = d["dataset"].split(" ", 1)[1]
    ov.text(735, y, plat, ha="right", va="center_baseline", fontsize=S.PT, color=S.INK)
    if plat == "Illumina":
        ov.text(95, y - 42, d["tumor"], ha="left", va="center_baseline", fontsize=S.PT + 0.5, fontweight="bold", color=S.INK)
    # a: stacked shares
    x = 0.0
    for v, lab, col, tc in zip(d["share"], d["share_lab"], COLS, ("white", S.INK, S.INK)):
        a.barh(y, v, left=x, height=104, color=col, edgecolor="white", linewidth=0.5, zorder=2)
        if v >= 6.5:
            a.text(x + v / 2, y, lab, ha="center", va="center_baseline", fontsize=FS, color=tc, zorder=3)
        x += v
    # b: retained on HPRC v2.1 (in graph above, not in graph below)
    for v, lab, col, dy in zip(d["keep"], d["keep_lab"], COLS[:2], (-29, 29)):
        b.barh(y + dy, v, height=54, color=col, edgecolor="none", zorder=2)
        b.text(v + 1.8, y + dy, lab, ha="left", va="center_baseline", fontsize=FS, color=S.INK, zorder=3)
    # c: totals
    for v, lab, col, dy in ((d["n"], d["n_lab"], S.C_LINEAR, -29), (d["d46"], d["d46_lab"], S.C_V21, 29)):
        c.barh(y + dy, v - 10, left=10, height=54, color=col, edgecolor="none", zorder=2)
        c.text(v * 1.13, y + dy, lab, ha="left", va="center_baseline", fontsize=FS, color=S.INK, zorder=3)

# titles and legends (positions of the ver4.5 image)
for k, title in (("a", "a  Class of GRCh38-graph FP"), ("b", "b  Retained on HPRC v2.1"), ("c", "c  Total FP by reference")):
    ov.text(PANELS[k][0], 95, title, ha="left", va="center_baseline", fontsize=S.PT + 0.5, fontweight="bold", color=S.INK)


def key(x, y, col, text):
    ov.add_patch(plt.Rectangle((x, y - 22), 44, 44, facecolor=col, edgecolor="none"))
    ov.text(x + 74, y, text, ha="left", va="center_baseline", fontsize=S.PT, color=S.INK)


for y, col, text in ((178, S.C_GERM_IN, "Germline, in HPRC v2.1 graph"), (240, S.C_GERM_OUT, "Germline, not in HPRC v2.1 graph"),
                     (302, S.C_NONGERM, "Non-germline")):
    key(755, y, col, text)
for y, col, text in ((208, S.C_LINEAR, "GRCh38-only"), (270, S.C_V21, "HPRC v2.1, no site features")):
    key(2850, y, col, text)
S.save(fig, "fig_R3.9_INDEL_FP_classes", "fig_R3.9_INDEL_FP_classes", "fig_R3.9_INDEL_FP_classes")

#!/usr/bin/env python
"""Check of the ver4 scripts on the OLD data: reproduce the letter's ver3 Response Figure 1.3.

No script of the ver3 figure exists. This check
  1. parses the HPRC v2.1 columns of Table R1.1 of the ver3 letter (the values the figure was drawn from), computes
     the 24 deltas of panel b and compares them with the numbers printed in the image (typed below from the image,
     the only place they exist);
  2. reads the same six October pairs (single-tumor SNV models on HPRC v2.1, recipe 7fc7396, seed 0) from their
     report.json with pathcount_data.load() and counts how many of the 48 values equal the letter table;
  3. redraws the figure from the letter-table values with pathcount_fig.draw() (same code path as the ver4 figure,
     without the seed rug) as check_october_redraw.png, and compares the colours of the heat cells and marks with
     the pixels of the original.

Run: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python check_october.py
"""
from pathlib import Path

import re

import numpy as np
from PIL import Image

import pathcount_data as D
import pathcount_fig as F

HERE = Path(__file__).resolve().parent
ORIG = HERE.parent / "ver3_inputs/image3.png"      # Response Figure 1.3 of Reviewer_Letter_Response_ver3.docx (word/media/image3.png)
R = D.R
ROWS = (("HG008", "Illumina", "HG008T Illumina"), ("HG008", "PacBio", "HG008T PacBio"), ("HG008", "ONT", "HG008T ONT"),
        ("COLO829", "Illumina", "COLO829T Illumina"), ("COLO829", "fiberseq", "COLO829T Fiber-seq"), ("COLO829", "ONT", "COLO829T ONT"))
# panel b of the ver3 image, read from the image: P, R, F1, F1 pre
IMAGE = ((0.031, 0.069, 0.042, 0.221), (0.007, 0.020, 0.011, 0.239), (-0.015, 0.018, -0.009, 0.220),
         (0.505, 0.037, 0.311, 0.420), (-0.017, 0.022, 0.008, 0.261), (0.085, 0.043, 0.062, 0.255))

# the figure was drawn from the v2.1 columns of Table R1.1 of the ver3 letter (3 decimals): parse that table
LETTER = Path("/wanglab/jshen/scratch_jshen/Pansoma_papers/ver1/Jiawei-Smaht-Paper/Pansoma_Resubmit/Reviewer_Letter_Response_ver3.md")
TABLE = {}
for line in LETTER.read_text().splitlines():
    c = [x.strip() for x in line.strip().strip("|").split("|")]
    key = c[0].replace("\\*", "").replace("*", "").strip() if c else ""
    if len(c) == 5 and key in {r[2] for r in ROWS} and key not in TABLE:
        cell = lambda x: [float(v) for v in re.findall(r"\d\.\d{3}", x)[:4]]
        TABLE[key] = (cell(c[3]), cell(c[4]))             # HPRC v2.1 with, without: P, R, F1, F1 before
assert len(TABLE) == 6, TABLE.keys()

rows, n_same, n_cells, differ = [], 0, 0, []
for (t, p, label), img in zip(ROWS, IMAGE):
    tw, two = TABLE[label]
    dt = [x - y for x, y in zip(tw, two)]
    assert all(abs(d - i) < 5e-4 for d, i in zip(dt, img)), (label, dt, img)   # letter table -> image, 24 of 24
    # the same pairs read today from the reports with pathcount_data.load()
    w = D.load(R / f"reviewer_cpu_20261013/{t}_{p}_SNV_d46_7fc7396/vcf_chr1/eval_bed/report.json")
    wo = D.load(R / f"ris_7fc7396_20261012/{t}_{p}_SNV_nopc_d46_7fc7396/vcf_chr1/eval_bed/report.json")
    assert w["pon"]["truth"] == wo["pon"]["truth"] == w["raw"]["truth"] == wo["raw"]["truth"]
    a = [w[st][k] for k, st, _ in D.METRICS]
    b = [wo[st][k] for k, st, _ in D.METRICS]
    same = [abs(round(x, 3) - y) < 5e-4 for x, y in zip(a + b, tw + two)]
    n_same += sum(same)
    n_cells += len(same)
    if not all(same):
        differ.append(label)
    print(f"{label:19s} truth in BED {w['pon']['truth']:5d}  reports today: with {' / '.join(f'{x:.3f}' for x in a)}  without "
          f"{' / '.join(f'{x:.3f}' for x in b)}  | letter table: with {' / '.join(f'{x:.3f}' for x in tw)}  without "
          f"{' / '.join(f'{x:.3f}' for x in two)}  {'same' if all(same) else 'DIFFERENT'}")
    rows.append({"label": label, "group": t, "with": tw, "without": two, "seeds": None, "delta": dt})
print("letter Table R1.1 (v2.1 columns) -> panel b of the image: 24 of 24 deltas reproduced")
print(f"pathcount_data.load() on today's October reports -> letter Table R1.1: {n_same} of {n_cells} values equal at 3 decimals; "
      f"rows with a different value: {differ or 'none'}")

outs, size, geo = F.draw(
    rows, HERE / "check_october_redraw", [(0.32, 0.925), (0.633, 0.875), (0.447, 0.808), (0.015, 0.735)],
    [[0.4, 0.6, 0.8], [0.7, 0.8], [0.5, 0.6, 0.7, 0.8], [0.2, 0.4, 0.6]], 0.5, "0.5",
    legend=[("circle", "With haplotype-count channel"), ("diamond", "Without haplotype-count channel")],
    delta_title=("Δ (with − without", "haplotype-count channel)"), formats=("png",))
print("redraw", outs[0], size, geo)

a = np.asarray(Image.open(ORIG).convert("RGB")).astype(int)
b = np.asarray(Image.open(outs[0]).convert("RGB")).astype(int)
assert a.shape == b.shape, (a.shape, b.shape)
print("same size as the original:", a.shape[1], "x", a.shape[0])
ys = geo["rows_y"]
worst = 0
for j, cx in enumerate(F.CELL_X):
    for y in ys:
        pa, pb = a[int(y) - 35, cx + 12], b[int(y) - 35, cx + 12]
        worst = max(worst, int(np.abs(pa - pb).max()))
print("heat cells: largest RGB channel difference to the original over 24 cells =", worst)
for name, hexc in (("blue", F.BLUE), ("orange", F.ORANGE), ("connector", F.CONNECT), ("guide", F.GUIDE)):
    c = np.array([int(hexc[i:i + 2], 16) for i in (1, 3, 5)])
    na, nb = int((np.abs(a - c).max(axis=2) == 0).sum()), int((np.abs(b - c).max(axis=2) == 0).sum())
    print(f"pixels of exact colour {name} {hexc}: original {na}, redraw {nb}")

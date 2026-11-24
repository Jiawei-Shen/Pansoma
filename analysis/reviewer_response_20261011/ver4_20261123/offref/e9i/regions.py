"""ver4 offref: merged GRCh38 regions (locus +- 600 bp) of every SNV call (unplaced PASS, ctrl_tp, ctrl_fp) of the three
HG008T six-tumor call sets (m6np__HG008T_*), for the remote HG008-N-D fetch. Copy of tmp/r33_hcc1395_20261104/e9i/regions.py."""
import csv
from pathlib import Path
D = Path("/scratch/jshen/data/reviewer_experiments_20261013/ver4_20261123/offref/DATA/runs")
iv = []
for run in sorted(D.glob("m6np__HG008T_*")):
    for r in csv.DictReader(open(run / "calls.tsv"), delimiter="\t"):
        if r["group"] in ("unplaced_pass", "ctrl_tp", "ctrl_fp") and r["kind_event"] == "SNP" and r.get("locus_chrom"):
            iv.append((r["locus_chrom"], max(0, int(r["locus_lo"]) - 600), int(r["locus_hi"]) + 600))
iv.sort(); m = []
for c, a, b in iv:
    if m and m[-1][0] == c and a <= m[-1][2]: m[-1][2] = max(m[-1][2], b)
    else: m.append([c, a, b])
Path("/scratch/jshen/data/reviewer_experiments_20261013/ver4_20261123/offref/e9i/HG008T/regions.bed").write_text("".join(f"{c}\t{a}\t{b}\n" for c, a, b in m))
print(len(iv), "loci ->", len(m), "regions")

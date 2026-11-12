#!/usr/bin/env python3
"""summary.md from overlap_long.tsv: per scope, normal, graph and type, chr1-22 and chr1 totals (normal, graph, shared,
normal-only, graph-only, shared as % of the normal's alleles)."""
import csv
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAMES = {"GRCh38": "GRCh38-only", "HPRC_v1.1_d9": "HPRC v1.1 (d9, >= 9 of 89)",
         "HPRC_v2.1_d46": "HPRC v2.1 (d46, >= 46 of 463)", "HPRC_v1.1_d9_VCF": "check: HPRC v1.1 d9 deconstruct VCF"}
t = defaultdict(lambda: [0, 0, 0])
for r in csv.DictReader(open(HERE / "overlap_long.tsv"), delimiter="\t"):
    for reg in ("chr1-22",) + (("chr1",) if r["chr"] == "chr1" else ()):
        k = (r["scope"], reg, r["normal"], r["graph"], r["variant_type"])
        for i, c in enumerate(("normal_total", "graph_total", "overlap")):
            t[k][i] += int(r[c])
md = ["# Germline alleles of the normals vs the graph alleles (distinct normalized keys; SNV, INDEL <= 50 bp)", ""]
for sc, what in (("all", "every chr1-22 allele of both sides"),
                 ("bed", "normal: FILTER PASS in its germline BED; graph: REF span in that BED")):
    for reg in ("chr1-22", "chr1"):
        md += [f"## scope {sc} ({what}), {reg}", "",
               "| normal | graph | type | normal alleles | graph alleles | shared | normal-only | graph-only | shared, % of normal |",
               "|---|---|---|---|---|---|---|---|---|"]
        for n in ("HG008N", "COLO829BL", "HCC1395BL"):
            for g in NAMES:
                for v in ("SNV", "INDEL"):
                    a, b, c = t[(sc, reg, n, g, v)]
                    md.append(f"| {n} | {NAMES[g]} | {v} | {a:,} | {b:,} | {c:,} | {a - c:,} | {b - c:,} | {100 * c / a:.1f}% |")
        md.append("")
(HERE / "summary.md").write_text("\n".join(md))
print("\n".join(md[:24]))

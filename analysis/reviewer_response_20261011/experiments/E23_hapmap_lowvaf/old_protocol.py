#!/usr/bin/env python3
"""E23 recall in the VAF bins of the submitted Supplementary Table 3 (1-2, 2-5, 5-10, 10-16.5 %), chr1, benchmark BED,
before population filtering (graph-carried truth included), over all truth SNVs and over those that reached a tensor.
The submitted table's recall denominator (96,393 SNVs) was the same for Pansoma and Mutect2, so it was a fixed truth
subset, not the truth with a Pansoma candidate (README, Results 6).
usage: python old_protocol.py RUN_DIR OUT.md
"""
import gzip
import sys
from collections import Counter
from pathlib import Path

from tensor_vaf import read_vcf_af

SET = Path("/scratch/jshen/data/pansoma_v2_tensors/HapMap_Illumina_nonHG005")
VCF = "/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/HapMap_somatic_snv_indel.vcf.gz"
BINS = [(0.01, 0.02, "1-2 %"), (0.02, 0.05, "2-5 %"), (0.05, 0.10, "5-10 %"), (0.10, 0.1651, "10-16.5 %")]
run, out = Path(sys.argv[1]), Path(sys.argv[2])
af = read_vcf_af(VCF)
table = {}
with open(SET / "truth" / "somatic.graph.tsv") as f:
    head = f.readline().rstrip("\n").split("\t"); col = {k: i for i, k in enumerate(head)}
    last = max(col[k] for k in ("truth_id", "chrom", "vcf_pos", "vcf_ref", "vcf_alt", "kind", "in_bed")) + 1
    for line in f:
        r = line.split("\t", last)
        if r[col["chrom"]] == "chr1" and r[col["kind"]] == "SNP" and r[col["in_bed"]] == "True":
            table[r[col["truth_id"]]] = ("chr1", int(r[col["vcf_pos"]]), r[col["vcf_ref"]], r[col["vcf_alt"]])
tp = set()
with gzip.open(run / "vcf_chr1" / "eval_bed" / "vcfeval" / "raw_calls" / "tp-baseline.vcf.gz", "rt") as f:
    for line in f:
        if not line.startswith("#"):
            c = line.split("\t", 5); tp.add((c[0], int(c[1]), c[3], c[4]))
n = Counter()
with open(SET / "tensors" / "somatic.recall.tsv") as f:
    head = f.readline().rstrip("\n").split("\t"); col = {k: i for i, k in enumerate(head)}
    for line in f:
        r = line.rstrip("\n").split("\t")
        key = table.get(r[col["truth_id"]])
        if key is None or r[col["passed"]] != "True":
            continue
        a = af.get(key)
        b = next((name for lo, hi, name in BINS if a is not None and lo <= a < hi), None)
        if b is None:
            continue
        has_tensor = r[col["status"]].startswith("tensor_")
        n[b, "all"] += 1; n[b, "all_tp"] += key in tp
        if has_tensor:
            n[b, "tensor"] += 1; n[b, "tensor_tp"] += key in tp
lines = [f"# E23 recall in the submitted Supplementary Table 3 bins: {run.name}, chr1 SNVs, BED, before population filtering", "",
         "| truth VAF | truth | recall | truth with a tensor | recall of those |", "|---|---|---|---|---|"]
for _, _, b in BINS:
    lines.append(f"| {b} | {n[b, 'all']:,} | {n[b, 'all_tp'] / max(n[b, 'all'], 1):.3f} | {n[b, 'tensor']:,} | "
                 f"{n[b, 'tensor_tp'] / max(n[b, 'tensor'], 1):.3f} |")
tot = {k: sum(n[b, k] for _, _, b in BINS) for k in ("all", "all_tp", "tensor", "tensor_tp")}
lines.append(f"| 1-16.5 % | {tot['all']:,} | {tot['all_tp'] / tot['all']:.3f} | {tot['tensor']:,} | {tot['tensor_tp'] / tot['tensor']:.3f} |")
lines += ["", "Submitted Supplementary Table 3 (Pansoma SNV, recall evaluation): 1-2 % 0.877 (17,147 truth), 2-5 % 0.966 "
          "(54,462), 5-10 % 0.935 (22,641), 10-16.5 % 0.738, 1-16.5 % 0.938 (96,393)."]
out.write_text("\n".join(lines) + "\n"); print("\n".join(lines))

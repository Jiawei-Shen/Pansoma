#!/bin/bash
# R1.m2 Supplementary Fig.: one chr1 SNV tensor whose rows follow their reads across node boundaries (HG008T PacBio
# HiFi, HPRC v1.1 graph with the d9 filter; exact somatic truth SNV chr1:242,529,453 G>A, candidate 3130240:147:SNP:G>A,
# chr1_labels.ndjson label 1, reason representative_allele_is_somatic_truth; 38 ALT, 40 REF, 0 other reads, AF 0.49).
# At column 71 (21 columns from the site, column 50) the reads cross a single-base bubble: all 38 ALT rows traverse the
# branch carried by 18 of the 90 haplotypes (channel 6 code 18), all 40 REF rows the GRCh38 branch carried by 72 (code 82).
# The earlier example (3133114:8:SNP:G>A) was only a residual partial match to the insertion chr1:242,680,268 G>GCTCGTTCCCA
# and was replaced on 2026-10-22. Renders with the repository's scripts/visualize_tensor.py (Python >= 3.10).
set -euo pipefail
D=/scratch/jshen/data/pansoma_v2_tensors/HG008T_PacBio/tensors/SNV
PY=${PY:-/wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python}
HERE=$(cd "$(dirname "$0")" && pwd)
$PY "$HERE/../../../../scripts/visualize_tensor.py" $D/chr1_shard_00000_data.npy --sample-index 32257 \
  --summary-path $D/chr1_variant_summary.ndjson --output "$HERE/fig_R1.m2_tensor_across_nodes.png" \
  --title "HG008T PacBio HiFi, HPRC v1.1 (d9): somatic SNV chr1:242,529,453 G>A (candidate 3130240:147)"

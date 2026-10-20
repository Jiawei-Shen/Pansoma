#!/bin/bash
# R1.m2 Supplementary Fig.: one chr1 SNV tensor whose rows follow their reads across node boundaries (HG008T PacBio
# HiFi, HPRC v1.1 graph with the d9 filter; somatic truth SNV chr1:242,680,268 G>A, candidate 3133114:8:SNP:G>A, AF 0.50).
# Chosen as the somatic tensor (chr1_labels.ndjson label 1) whose reads split most evenly between two nodes inside the
# window: in columns 30-37 the 49 ALT reads traverse a node carried by 10-12 of the 90 haplotypes (channel 6 codes 10 / 12),
# the 47 REF reads one carried by 89 (code 99). Renders with the repository's scripts/visualize_tensor.py (Python >= 3.10).
set -euo pipefail
D=/scratch/jshen/data/pansoma_v2_tensors/HG008T_PacBio/tensors/SNV
PY=${PY:-/wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python}
HERE=$(cd "$(dirname "$0")" && pwd)
$PY "$HERE/../../../../scripts/visualize_tensor.py" $D/chr1_shard_00000_data.npy --sample-index 32328 \
  --summary-path $D/chr1_variant_summary.ndjson --output "$HERE/fig_R1.m2_tensor_across_nodes.png" \
  --title "HG008T PacBio HiFi, HPRC v1.1 (d9): somatic SNV chr1:242,680,268 G>A (candidate 3133114:8)"

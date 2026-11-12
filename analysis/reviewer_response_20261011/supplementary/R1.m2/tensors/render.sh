#!/bin/bash
# Redraw fig_R1.m2_tensor_across_nodes.png and indel/fig_R1.m2_INDEL_span_01.png from the tensors in this directory
# (extract.py), with the same titles and panel order as ../make_figure.sh and ../indel/make_figures.sh.
# Usage: render.sh [output directory, default: this directory]. Python >= 3.10 with numpy + matplotlib.
set -euo pipefail
PY=${PY:-/wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python}
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-$HERE}
VIS="$HERE/../../../../../scripts/visualize_tensor.py"
draw() {  # name, output png, title
  $PY "$VIS" "$HERE/$1.npy" --sample-index 0 --shard-index 0 --summary-path "$HERE/$1.summary.ndjson" \
    --panel-order 0,5,6,1,2,3,4,7 --output "$OUT/$2" --title "$3" > /dev/null
}
draw R1.m2_tensor_across_nodes fig_R1.m2_tensor_across_nodes.png \
  "HG008T PacBio HiFi, HPRC v1.1 (d9): somatic SNV chr1:242,529,453 G>A (candidate 3130240:147)"
draw R1.m2_INDEL_span_01 fig_R1.m2_INDEL_span_01.png \
  "HG008T Illumina, HPRC v2.1 (d46): somatic 8-bp DEL chr1:23,017,177 TTATATATA>T (candidate 1811743:0)"

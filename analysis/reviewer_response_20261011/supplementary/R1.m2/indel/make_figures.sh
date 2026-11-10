#!/bin/bash
# R1.m2 Supplementary Figs. (INDEL): ten chr1 somatic truth INDEL tensors whose windows span several graph nodes.
# Steps (read-only on the tensor sets; outputs in this directory):
#   select_examples.py -> candidates.tsv  every exact, PASS, in-BED somatic truth INDEL tensor of the twelve HPRC sets
#   choose.py          -> picks.tsv       the ten examples (rules in choose.py)
#   annotate.py        -> annotations.md  caption facts (nodes under the window, decoded path counts per branch column)
#   visualize_tensor.py -> fig_R1.m2_INDEL_NN.png, panels in the order 0,5,6,1,2,3,4,7 (graph base second, path count third)
# Renders with the repository's scripts/visualize_tensor.py (Python >= 3.10). SKIP_SELECT=1 reuses candidates.tsv.
set -euo pipefail
PY=${PY:-/wanglab/jshen/anaconda3/envs/hunyuanvideo15/bin/python}
HERE=$(cd "$(dirname "$0")" && pwd)
D=/scratch/jshen/data/pansoma_v2_tensors
cd "$HERE"
[ "${SKIP_SELECT:-0}" = 1 ] || $PY select_examples.py candidates.tsv 2> select_examples.log
$PY choose.py candidates.tsv picks.tsv
$PY annotate.py picks.tsv > annotations.md
tail -n +2 picks.tsv | while IFS=$'\t' read -r fig set label cid shard index truth event len af rest; do
  $PY ../../../../../scripts/visualize_tensor.py "$D/$set/tensors/INDEL/$shard" --sample-index "$index" \
    --summary-path "$D/$set/tensors/INDEL/chr1_variant_summary.ndjson" --panel-order 0,5,6,1,2,3,4,7 \
    --output "$HERE/fig_$fig.png" --title "$label: somatic $event $truth (candidate $cid)" > /dev/null
  echo "$fig $set $cid"
done

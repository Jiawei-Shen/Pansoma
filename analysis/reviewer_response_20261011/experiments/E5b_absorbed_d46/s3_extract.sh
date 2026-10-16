#!/bin/bash
# E5b step 3: cut the subgraph around every truth window out of a graph's GBZ (nodes within 250 bp of a GRCh38 window
# node, and the edges among them), as GFA. Usage: s3_extract.sh d9|d46   (Slurm: tmp/reviewer_cpu_20261013/E5b_absorbed_d46)
set -euo pipefail
g=$1
D=/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46
case $g in
  d9)  VG=/scratch/wzhang/bin/vg_v1.65.0; GBZ=/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz ;;   # vg 1.77 find crashes on the d9 GBZ
  d46) VG=/scratch/wzhang/bin/vg_v1.65.0; GBZ=/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/hprc-v2.1-mc-grch38.d46.gbz ;;  # vg 1.77 find crashed on d46 too (job 387310, get_position_of_step)
  *) echo "graph d9|d46" >&2; exit 2 ;;
esac
$VG version | head -1
/usr/bin/time -v $VG find -x $GBZ -N $D/window_nodes_$g.txt -c 250 -L > $D/subgraph_$g.vg.tmp
mv $D/subgraph_$g.vg.tmp $D/subgraph_$g.vg
$VG convert -f $D/subgraph_$g.vg | awk -F'\t' '$1=="S"||$1=="L"' > $D/subgraph_$g.gfa.tmp
mv $D/subgraph_$g.gfa.tmp $D/subgraph_$g.gfa
echo "S $(grep -c '^S' $D/subgraph_$g.gfa) L $(grep -c '^L' $D/subgraph_$g.gfa)"

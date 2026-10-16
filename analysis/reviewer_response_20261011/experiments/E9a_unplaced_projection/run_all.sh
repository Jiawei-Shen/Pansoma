#!/bin/bash
# E9a/E9e over every finished CPU re-prediction: each RUN/vcf_chr1* under reviewer_cpu_20261013 whose eval_nobed/report.json
# exists, on the d9 or d46 graph, without <scratch>/runs/RUN/VCFDIR/metrics.json yet -> one Slurm array (at most 6 tasks
# at a time). Linear-graph runs are skipped (no off-reference nodes; collect.py reads their own reports).
# usage: bash run_all.sh [--dry-run]     (E9A_RUNS=<dir> / E9A_OUT=<dir>: another runs / output directory, e.g. for tests)
set -euo pipefail
R=${E9A_RUNS:-/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013}
S=/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E9a_unplaced_projection
J=/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E9a_unplaced_projection
for g in d9 d46; do
  # CHM13 (E9b) dropped by the user 2026-10-13: only the edges are required; projection writes chm13_class not_built
  [ -s $S/graph/$g/edges/meta.json ] || { echo "graph $g not prepared ($S/graph/$g): run prep_graph.sh first"; exit 1; }
done
LIST=$J/tasks.$(date +%Y%m%d_%H%M%S).txt
: > $LIST
for d in $R/*/vcf_chr1*; do
  [ -s $d/eval_nobed/report.json ] || continue
  RUN=$(basename $(dirname $d)); VDIR=$(basename $d)
  [ -s ${E9A_OUT:-$S/runs}/$RUN/$VDIR/metrics.json ] && continue
  G=$(zcat $d/*.graph.vcf.gz | head -200 | grep -m1 '^##pansoma_graph_index=' || true)
  case $G in *d9.graph_index*|*d46.graph_index*) echo "$RUN $VDIR" >> $LIST ;; esac
done
N=$(wc -l < $LIST)
echo "$N run/vcf dirs to process ($LIST)"
[ $N = 0 ] && exit 0
cat $LIST
[ "${1:-}" = "--dry-run" ] && exit 0
cd $J && sbatch --array=0-$((N - 1))%6 -J E9a_task -o $J/logs/task.%A_%a.log --export=ALL $E/task.sh $LIST

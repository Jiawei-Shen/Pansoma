#!/bin/bash
# E21: finish every table once the index-build array (build_index.sbatch, job 387137) is done. Login node, about 10 min.
#   bash /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E21_training_composition/collect.sh
set -euo pipefail
HERE=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E21_training_composition
PY=/wanglab/jshen/anaconda3/bin/python
E13=/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/e13
cd "$HERE"
grep -h "tensors,\|Error\|Traceback\|Killed" /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E21_training_composition/logs/index_*.log | sort | uniq -c | awk '$1>0' | wc -l
$PY collect_runs.py > /dev/null           # T1: runs_training.tsv / .md
$PY summarize.py                          # T2-T6 TSVs and runs_replay_check.tsv
for s in HG008T_PacBio_d46 HG008T_PacBio; do  # E13 learning-curve subsets (exclude files for train.py)
  $PY make_e13_subsets.py $s SNV "$E13" 0.10 0.25 0.50 0.75 --replicates 2
done
$PY make_tables.py                        # tables.md

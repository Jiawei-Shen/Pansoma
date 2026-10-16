#!/bin/bash
# E16_bootstrap collect: block-bootstrap CIs (P / R / F1 after the PoN) and paired F1 differences
# for every finished chr1 call set under /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013 (discovery:
# ../E16_bootstrap/callsets.py; comparators: ../E16_bootstrap/comparators.tsv). Rerun whenever more runs finish; it
# processes whatever exists and rewrites: tables.md ci.tsv paired.tsv summary.json.
# Light (1 core, < 4G, a few minutes): run on the login node. Read-only on the run directories.
#   bash collect.sh            here
#   bash collect.sh --slurm    as one small Slurm job (-p general -c 2 --mem 4G); log in tmp/reviewer_cpu_20261013/E16_bootstrap/
set -euo pipefail
HERE=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E16_bootstrap
T=/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E16_bootstrap
PY=/wanglab/jshen/anaconda3/bin/python
export PYTHONDONTWRITEBYTECODE=1
if [ "${1:-}" = --slurm ]; then
  mkdir -p $T/logs
  sbatch -p general -c 2 --mem=4G -t 1:00:00 -J E16_collect -o $T/logs/E16_collect-%j.out --wrap "bash $HERE/collect.sh"
  exit 0
fi
cd $HERE
/usr/bin/time -f "E16_bootstrap: %es, max RSS %MKB" $PY e16_bootstrap.py

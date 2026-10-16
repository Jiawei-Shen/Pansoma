#!/bin/bash
# E17a_indel_strata collect: INDEL TP / FP / FN, recall and precision by INS / DEL, length and repeat context
# for every finished chr1 call set under /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013 (discovery:
# ../E16_bootstrap/callsets.py; comparators: ../E16_bootstrap/comparators.tsv). Rerun whenever more runs finish; it
# processes whatever exists and rewrites: tables.md strata.tsv summary.json.
# Light (1 core, < 2G, a few minutes): run on the login node. Read-only on the run directories.
#   bash collect.sh            here
#   bash collect.sh --slurm    as one small Slurm job (-p general -c 2 --mem 2G); log in tmp/reviewer_cpu_20261013/E17a_indel_strata/
set -euo pipefail
HERE=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E17a_indel_strata
T=/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E17a_indel_strata
PY=/wanglab/jshen/anaconda3/bin/python
export PYTHONDONTWRITEBYTECODE=1
if [ "${1:-}" = --slurm ]; then
  mkdir -p $T/logs
  sbatch -p general -c 2 --mem=2G -t 1:00:00 -J E17a_collect -o $T/logs/E17a_collect-%j.out --wrap "bash $HERE/collect.sh"
  exit 0
fi
cd $HERE
/usr/bin/time -f "E17a_indel_strata: %es, max RSS %MKB" $PY e17a_indel_strata.py

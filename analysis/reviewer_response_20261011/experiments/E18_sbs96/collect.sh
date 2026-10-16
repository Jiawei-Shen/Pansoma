#!/bin/bash
# E18_sbs96 collect: SBS96 spectra of truth, TP, FP, FN and germline; cosine to truth / COSMIC
# for every finished chr1 call set under /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013 (discovery:
# ../E16_bootstrap/callsets.py; comparators: ../E16_bootstrap/comparators.tsv). Rerun whenever more runs finish; it
# processes whatever exists and rewrites: tables.md sbs96.tsv summary.json.
# Light (1 core, < 4G, a few minutes): run on the login node. Read-only on the run directories.
#   bash collect.sh            here
#   bash collect.sh --slurm    as one small Slurm job (-p general -c 2 --mem 4G); log in tmp/reviewer_cpu_20261013/E18_sbs96/
set -euo pipefail
HERE=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E18_sbs96
T=/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E18_sbs96
PY=/wanglab/jshen/anaconda3/bin/python
export PYTHONDONTWRITEBYTECODE=1
if [ "${1:-}" = --slurm ]; then
  mkdir -p $T/logs
  sbatch -p general -c 2 --mem=4G -t 1:00:00 -J E18_collect -o $T/logs/E18_collect-%j.out --wrap "bash $HERE/collect.sh"
  exit 0
fi
cd $HERE
/usr/bin/time -f "E18_sbs96: %es, max RSS %MKB" $PY e18_sbs96.py

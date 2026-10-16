#!/bin/bash
# E1.6 recall funnel: rewrites tables.md and summary.json in this folder and
# /scratch/jshen/data/reviewer_experiments_20261013/E1.6_recall_funnel/per_truth.tsv.gz, from the in-sample call sets under
# /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013 (plus the COLO829T fiberseq joint d9 model at P100).
# Light (about 30 s and under 1 GB once the E2.4 truth-AF caches exist): run on the login node (bash collect.sh) or as a
# Slurm job (-p general, 2 CPUs):
#   sbatch -J E1.6_recall_funnel -o /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E1.6/logs/%x-%j.out \
#     /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E1.6/job.sh
set -euo pipefail
P=/wanglab/jshen/anaconda3/bin/python
export PYTHONDONTWRITEBYTECODE=1
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments
/usr/bin/time -f "%es %MKB" $P $E/E1.6_recall_funnel/e16_recall_funnel.py

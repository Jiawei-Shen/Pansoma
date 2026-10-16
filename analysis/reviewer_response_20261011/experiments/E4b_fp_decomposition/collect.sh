#!/bin/bash
# E4b_fp_decomposition: classify the FPs that survive the PoN for every finished call set under
# /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013 and rewrite tables.md / summary.json here.
# Per-call outputs are cached in /scratch/jshen/data/reviewer_experiments_20261013/call_level/ and recomputed only when
# the call set's rtg files are newer, so a rerun processes the new call sets only.
# Run: on the login node when few call sets are new (bash collect.sh), else as a Slurm job:
#   sbatch -J E4b_fp_decomposition -o /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E4b_fp_decomposition/logs/%x-%j.out /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E4b_fp_decomposition/job.sh
# or all four analyses at once: bash ../call_level_common/collect_all.sh
set -euo pipefail
P=/wanglab/jshen/anaconda3/bin/python
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments
$P $E/call_level_common/prep_hprc.py   # HPRC v1.1 chr1 allele table (skipped when it exists)
/usr/bin/time -f "%es %MKB" $P $E/E4b_fp_decomposition/fp_decomposition.py

#!/bin/bash
# E19_germline_density: P / R / F1 by germline density for every finished call set under
# /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013 and rewrite tables.md / summary.json here.
# Per-call outputs are cached in /scratch/jshen/data/reviewer_experiments_20261013/call_level/ and recomputed only when
# the call set's rtg files are newer, so a rerun processes the new call sets only.
# Run: on the login node when few call sets are new (bash collect.sh), else as a Slurm job:
#   sbatch -J E19_germline_density -o /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E19_germline_density/logs/%x-%j.out /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E19_germline_density/job.sh
# or all four analyses at once: bash ../call_level_common/collect_all.sh
set -euo pipefail
P=/wanglab/jshen/anaconda3/bin/python
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments
/usr/bin/time -f "%es %MKB" $P $E/E19_germline_density/germline_density.py

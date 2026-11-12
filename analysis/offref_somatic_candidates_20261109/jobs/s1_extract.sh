#!/bin/bash
# s1: read-consensus contexts of the off-reference tensors (and on-reference controls) of one (set, kind, chromosome).
# $1 = task list (one "SET KIND CHROM [BASE_SET]" per line), task = SLURM_ARRAY_TASK_ID (0-based line).
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 12:00:00
set -euo pipefail
A=/scratch/jshen/Github/Pansoma/analysis/offref_somatic_candidates_20261109
D=/scratch/jshen/data/pansoma_net_v2_runs/offref_somatic_candidates_20261109/s1
read S K C B < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$1")
echo "task $SLURM_ARRAY_TASK_ID $S $K $C host $(hostname) job ${SLURM_JOB_ID:-}"
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python -I -B $A/s1_extract.py $S $K $C $D $B

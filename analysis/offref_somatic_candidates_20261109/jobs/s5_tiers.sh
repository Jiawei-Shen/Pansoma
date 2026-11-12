#!/bin/bash
# s5: tiers of one tumor ($1).
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH -t 4:00:00
set -euo pipefail
A=/scratch/jshen/Github/Pansoma/analysis/offref_somatic_candidates_20261109
echo "host $(hostname) job ${SLURM_JOB_ID:-}"
/usr/bin/time -f "s5_tiers %es %MKB" /wanglab/jshen/anaconda3/bin/python -I -B $A/s5_tiers.py $1

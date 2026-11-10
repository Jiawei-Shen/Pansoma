#!/bin/bash
# s2a: unique queries of one tumor ($1) from the 132 s1 tables.
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 4:00:00
set -euo pipefail
A=/scratch/jshen/Github/Pansoma/analysis/offref_somatic_candidates_20261109
echo "host $(hostname) job ${SLURM_JOB_ID:-}"
/usr/bin/time -f "s2_queries %es %MKB" /wanglab/jshen/anaconda3/bin/python -I -B $A/s2_queries.py $1

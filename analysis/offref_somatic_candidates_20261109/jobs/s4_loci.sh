#!/bin/bash
# s4a: loci and k-mer pairs of one tumor ($1).
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 12:00:00
set -euo pipefail
A=/scratch/jshen/Github/Pansoma/analysis/offref_somatic_candidates_20261109
echo "host $(hostname) job ${SLURM_JOB_ID:-}"
/usr/bin/time -f "s4_loci %es %MKB" /wanglab/jshen/anaconda3/bin/python -I -B $A/s4_loci.py $1

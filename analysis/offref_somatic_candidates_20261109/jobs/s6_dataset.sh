#!/bin/bash
# s6: off-reference training labels -> <SET>_offref tensor sets; $1 TUMOR [--dry-run] [--relabel-test].
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 8:00:00
set -euo pipefail
A=/scratch/jshen/Github/Pansoma/analysis/offref_somatic_candidates_20261109
echo "host $(hostname) job ${SLURM_JOB_ID:-} args $*"
/usr/bin/time -f "s6_dataset %es %MKB" /wanglab/jshen/anaconda3/bin/python -I -B $A/s6_dataset.py "$@"

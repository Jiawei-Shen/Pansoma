#!/bin/bash
# s4d: fetch regions of the loci a label has not counted yet; $1 TUMOR $2 LABEL $3 OUT_BED [$4 SELECT].
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 1:00:00
set -euo pipefail
A=/scratch/jshen/Github/Pansoma/analysis/offref_somatic_candidates_20261109
/wanglab/jshen/anaconda3/bin/python -I -B $A/s4_regions.py "$@"

#!/bin/bash
# s9_integrate.py: per-locus integration (16 workers, one sqlite connection each); outputs in $D
#SBATCH -p general,gpu
#SBATCH -c 16
#SBATCH --mem=12G
#SBATCH -t 2:00:00
#SBATCH -J s9_integrate
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s9_integrate.%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
/wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/s9_integrate.py 16

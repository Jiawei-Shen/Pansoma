#!/bin/bash
#SBATCH -J miss_ont
#SBATCH -p general
#SBATCH -c 40
#SBATCH --mem=60G
#SBATCH -t 12:00:00
#SBATCH -o reads_part%a_%A.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd /scratch/jshen/Github/Pansoma/tmp/somatic_miss_analysis_ont
PART="${SLURM_ARRAY_TASK_ID}/6" OUTPUT="reads_part${SLURM_ARRAY_TASK_ID}.jsonl" /wanglab/jshen/anaconda3/bin/python reads.py

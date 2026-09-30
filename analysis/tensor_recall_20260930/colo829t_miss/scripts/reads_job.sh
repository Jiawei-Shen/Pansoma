#!/bin/bash
# sbatch --array=0-(N-1) -c CPUS --mem=MEM -o ../<platform>/reads_part%a_%A.out --export=ALL,SAMPLE=COLO829T_<platform>,PARTS=N reads_job.sh
#SBATCH -J miss_colo
#SBATCH -p general
#SBATCH -t 12:00:00
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd /scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930/colo829t_miss/scripts
PART="${SLURM_ARRAY_TASK_ID}/${PARTS}" OUTPUT="reads_part${SLURM_ARRAY_TASK_ID}.jsonl" /wanglab/jshen/anaconda3/bin/python reads.py

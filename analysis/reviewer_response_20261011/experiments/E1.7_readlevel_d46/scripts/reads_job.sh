#!/bin/bash
# E1.7 read-level classes on the current tensor sets (port of tensor_recall_20260930/colo829t_miss/scripts/reads_job.sh).
# sbatch --array=0-(N-1) -c CPUS --mem=MEM -o <tmp>/logs/reads_<SAMPLE>_%a_%A.out --export=ALL,SAMPLE=<set>,PARTS=N reads_job.sh
#SBATCH -J e17_reads
#SBATCH -p general
#SBATCH -t 24:00:00
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E1.7_readlevel_d46/scripts
PART="${SLURM_ARRAY_TASK_ID}/${PARTS}" OUTPUT="reads_part${SLURM_ARRAY_TASK_ID}.jsonl" /wanglab/jshen/anaconda3/bin/python reads.py

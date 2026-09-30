#!/bin/bash
# sbatch -c CPUS --mem=MEM -o ../<platform>/snv_%j.out --export=ALL,SAMPLE=COLO829T_<platform> snv_job.sh
#SBATCH -J snv_colo
#SBATCH -p general
#SBATCH -t 12:00:00
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd /scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930/colo829t_miss/scripts
/wanglab/jshen/anaconda3/bin/python snv_nodes.py
/wanglab/jshen/anaconda3/bin/python snv_detail.py

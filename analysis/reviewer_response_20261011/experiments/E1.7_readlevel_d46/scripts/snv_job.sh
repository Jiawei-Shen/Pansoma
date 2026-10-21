#!/bin/bash
# sbatch -c CPUS --mem=MEM -o <tmp>/logs/snv_<SAMPLE>_%j.out --export=ALL,SAMPLE=<set> snv_job.sh   (after reads_job.sh)
#SBATCH -J e17_snv
#SBATCH -p general
#SBATCH -t 12:00:00
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E1.7_readlevel_d46/scripts
/wanglab/jshen/anaconda3/bin/python snv_nodes.py
/wanglab/jshen/anaconda3/bin/python snv_detail.py

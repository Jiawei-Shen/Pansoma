#!/bin/bash
#SBATCH -J ont_snv
#SBATCH -p general
#SBATCH -c 24
#SBATCH --mem=30G
#SBATCH -t 8:00:00
#SBATCH -o snv_%j.out
cd /scratch/jshen/Github/Pansoma/tmp/somatic_miss_analysis_ont
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
/wanglab/jshen/anaconda3/bin/python snv_nodes.py
/wanglab/jshen/anaconda3/bin/python snv_detail.py

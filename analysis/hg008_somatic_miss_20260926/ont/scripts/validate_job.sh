#!/bin/bash
#SBATCH -J ont_trimval
#SBATCH -p general
#SBATCH -c 24
#SBATCH --mem=32G
#SBATCH -t 6:00:00
#SBATCH -o validate_trim_%j.out
cd /scratch/jshen/Github/Pansoma/tmp/somatic_miss_analysis_ont
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
/wanglab/jshen/anaconda3/bin/python validate_trim.py 24

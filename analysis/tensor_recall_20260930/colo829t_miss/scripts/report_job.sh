#!/bin/bash
# sbatch -o ../<platform>/report_%j.out --export=ALL,SAMPLE=COLO829T_<platform> report_job.sh
#SBATCH -J report_colo
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 2:00:00
cd /scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930/colo829t_miss/scripts
/wanglab/jshen/anaconda3/bin/python report.py

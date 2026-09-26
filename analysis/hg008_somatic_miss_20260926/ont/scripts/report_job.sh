#!/bin/bash
#SBATCH -J ont_report
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH -o report_%j.out
cd /scratch/jshen/Github/Pansoma/tmp/somatic_miss_analysis_ont
/wanglab/jshen/anaconda3/bin/python report.py > report.log 2>&1 && python3 tables.py > tables.md 2>tables.err && echo REPORT_DONE

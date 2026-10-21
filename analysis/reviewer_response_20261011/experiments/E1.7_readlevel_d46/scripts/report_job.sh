#!/bin/bash
# sbatch -o <tmp>/logs/report_<SAMPLE>_%j.out --export=ALL,SAMPLE=<set> report_job.sh   (after snv_job.sh)
#SBATCH -J e17_report
#SBATCH -p general
#SBATCH -t 4:00:00
#SBATCH -c 2
#SBATCH --mem=8G
set -euo pipefail
cd /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E1.7_readlevel_d46/scripts
/wanglab/jshen/anaconda3/bin/python report.py

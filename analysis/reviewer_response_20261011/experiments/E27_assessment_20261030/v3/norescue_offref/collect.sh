#!/bin/bash
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 1:00:00
set -euo pipefail
/wanglab/jshen/anaconda3/bin/python -u /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_norescue/code/collect.py
/wanglab/jshen/anaconda3/bin/python -u /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_norescue/code_e9i/report.py

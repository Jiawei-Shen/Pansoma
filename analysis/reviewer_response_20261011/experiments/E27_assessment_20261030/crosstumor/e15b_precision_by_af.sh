#!/bin/bash
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -J e15b_prec_af
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/crosstumor/logs/e15b_prec_af-%j.out
cd /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/crosstumor
/wanglab/jshen/anaconda3/bin/python -I e15b_precision_by_af.py > e15b_precision_by_af.md

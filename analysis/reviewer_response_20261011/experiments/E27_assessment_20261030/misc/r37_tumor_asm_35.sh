#!/bin/bash
#SBATCH -J r37_asm35
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=28G
#SBATCH -t 2:00:00
#SBATCH -w meow
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/misc/r37_tumor_asm_35.%j.log
set -euo pipefail
/wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/misc/r37_tumor_asm_35.py /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/misc/r37_tumor_asm_35.tsv

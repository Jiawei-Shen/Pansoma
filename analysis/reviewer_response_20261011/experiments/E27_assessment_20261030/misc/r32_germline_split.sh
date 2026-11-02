#!/bin/bash
#SBATCH -J r32_split
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 3:00:00
#SBATCH -w miao
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/misc/r32_germline_split.%j.log
set -euo pipefail
/wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/misc/r32_germline_split.py /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/misc/r32_germline_split.tsv

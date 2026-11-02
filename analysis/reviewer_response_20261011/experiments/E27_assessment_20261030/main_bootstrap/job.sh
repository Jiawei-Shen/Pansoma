#!/bin/bash
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -J assess_bootstrap
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/main_bootstrap/logs/%x-%j.out
set -euo pipefail
cd /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/main_bootstrap
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python pairs.py

#!/bin/bash
#SBATCH -J assess_label_status
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 2:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators/logs/label_status-%j.out
set -euo pipefail
cd /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators
/usr/bin/time -v /wanglab/jshen/anaconda3/bin/python label_status.py

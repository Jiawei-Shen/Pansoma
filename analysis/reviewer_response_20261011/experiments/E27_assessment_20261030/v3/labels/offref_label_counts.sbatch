#!/bin/bash
#SBATCH -J assess_offref_counts
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/v3/labels/logs/offref_label_counts-%j.out
set -euo pipefail
cd /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/v3/labels
PYTHONDONTWRITEBYTECODE=1 /wanglab/jshen/anaconda3/bin/python offref_label_counts.py

#!/bin/bash
#SBATCH -J assess_label_arms
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators/logs/label_status_arms-%j.out
set -euo pipefail
cd /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators
ARM=linear /wanglab/jshen/anaconda3/bin/python label_status.py
ARM=d9 /wanglab/jshen/anaconda3/bin/python label_status.py

#!/bin/bash
#SBATCH -J s8_dipcall
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 1:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s8_dipcall-%j.out
set -euo pipefail
/usr/bin/time -v /wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/s8_dipcall.py

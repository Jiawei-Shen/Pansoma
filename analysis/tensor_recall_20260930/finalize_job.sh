#!/bin/bash
# After the COLO829T ONT class report and the filtered scan: regenerate tables.md, per_truth_*.tsv and pon/pon_tables.md.
#SBATCH -J recall_tables
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 1:00:00
#SBATCH -o finalize_%j.out
cd /scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930
/wanglab/jshen/anaconda3/bin/python recall_tables.py > /dev/null && (cd pon && /wanglab/jshen/anaconda3/bin/python pon_tables.py > /dev/null) && echo "tables regenerated $(date)"

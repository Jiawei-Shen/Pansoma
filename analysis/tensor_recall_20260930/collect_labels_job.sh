#!/bin/bash
#SBATCH -J recall_labels
#SBATCH -p general
#SBATCH -c 6
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -o labels/collect_labels_%j.out
cd /scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930
/wanglab/jshen/anaconda3/bin/python collect_labels.py

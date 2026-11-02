#!/bin/bash
#SBATCH -J q3_valsplit
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q3/q3_%j.out
cd /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q3
/wanglab/jshen/anaconda3/envs/sambcfenv/bin/python -I val_split.py /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q3/val_split.json

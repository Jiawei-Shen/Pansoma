#!/bin/bash
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J q4_e16_all
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q4/q4_%j.out
cd /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q4
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python e16_pairs_q4.py

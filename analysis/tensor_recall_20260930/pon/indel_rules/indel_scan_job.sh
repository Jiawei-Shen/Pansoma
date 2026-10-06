#!/bin/bash
# sbatch --array=0-21 indel_scan_job.sh   (chr1-22)
#SBATCH -J pon_indel
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 4:00:00
#SBATCH -o scan/indel_scan_%a_%A.out
cd /scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930/pon/indel_rules
/wanglab/jshen/anaconda3/bin/python indel_scan.py chr$((SLURM_ARRAY_TASK_ID + 1))

#!/bin/bash
# sbatch --array=0-23 pon_scan_job.sh   (chr1-22, chrX, chrY; chr1 also writes scan/chr1.keys.tsv)
#SBATCH -J pon_scan
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 6:00:00
#SBATCH -o scan/pon_scan_%a_%A.out
cd /scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930/pon
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
C=${CHROMS[$SLURM_ARRAY_TASK_ID]}
[ "$C" = chr1 ] && export DUMP_KEYS=1
/wanglab/jshen/anaconda3/bin/python pon_scan.py $C

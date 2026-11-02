#!/bin/bash
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 3:00:00
set -euo pipefail
echo "rcpu__COLO829_Illumina_SNV_7fc7396_ddp4_b1024__vcf_chr1 $SLURM_ARRAY_TASK_ID/6 host $(hostname) job ${SLURM_JOB_ID:-}"
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_ep8/code_e9i/normal_reads.py rcpu__COLO829_Illumina_SNV_7fc7396_ddp4_b1024__vcf_chr1 /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_ep8/E9i_out $SLURM_ARRAY_TASK_ID 6

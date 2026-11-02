#!/bin/bash
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 12:00:00
set -euo pipefail
echo "host $(hostname) job ${SLURM_JOB_ID:-}"
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python -u /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_ep8/code/e9c_run.py /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/COLO829_Illumina_SNV_7fc7396_ddp4_b1024/vcf_chr1 rcpu__COLO829_Illumina_SNV_7fc7396_ddp4_b1024__vcf_chr1

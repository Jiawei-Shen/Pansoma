#!/bin/bash
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 6:00:00
#SBATCH --array=0-11
set -euo pipefail
RUN=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_norescue/runs.txt)
echo "$RUN host $(hostname) job ${SLURM_JOB_ID:-}"
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python -u /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_norescue/code/e9c_run.py /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/${RUN}__norescue/vcf_chr1 rcpu__${RUN}__norescue__vcf_chr1

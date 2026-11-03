#!/bin/bash
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 3:00:00
#SBATCH --array=0-35
set -euo pipefail
O=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_norescue
R=$((6 + SLURM_ARRAY_TASK_ID / 6)); I=$((SLURM_ARRAY_TASK_ID % 6))
RUN=$(sed -n "$((R+1))p" $O/runs.txt)
echo "$RUN shard $I/6 host $(hostname) job ${SLURM_JOB_ID:-}"
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python $O/code_e9i/normal_reads.py rcpu__${RUN}__norescue__vcf_chr1 $O/E9i_out $I 6

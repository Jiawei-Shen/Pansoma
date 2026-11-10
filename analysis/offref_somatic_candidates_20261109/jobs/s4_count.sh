#!/bin/bash
# s4b: k-mer read counts; $1 TUMOR, $2 LABEL, $3 BAM[,BAM], $4 NSHARDS [$5 SELECT $6 TAG]; array task = shard.
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -t 6:00:00
set -euo pipefail
A=/scratch/jshen/Github/Pansoma/analysis/offref_somatic_candidates_20261109
echo "host $(hostname) job ${SLURM_JOB_ID:-} $1 $2 shard ${SLURM_ARRAY_TASK_ID}/$4"
/usr/bin/time -f "s4_count %es %MKB" /wanglab/jshen/anaconda3/bin/python -I -B $A/s4_count.py $1 $2 $3 $SLURM_ARRAY_TASK_ID $4 ${5:-} ${6:-}

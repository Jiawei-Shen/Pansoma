#!/bin/bash
# E9h step 1: A1-read consensus queries + event windows (build_queries.py), read-only on the tensor sets.
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 2:00:00
set -euo pipefail
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E9h_offref_key_reads
D=/scratch/jshen/data/reviewer_experiments_20261013/E9h_offref_key_reads
echo "host $(hostname) job ${SLURM_JOB_ID:-} start $(date)"
/usr/bin/time -f "build %es %MKB" /wanglab/jshen/anaconda3/bin/python $E/build_queries.py $D
echo "end $(date)"

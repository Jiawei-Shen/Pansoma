#!/bin/bash
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -J q6_side
#SBATCH --array=1-12%4
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/side/logs/%x_%A_%a.out
set -euo pipefail
Q=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/side
read S C < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $Q/sets.tsv)
echo "host $(hostname) $S $(date)"
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python -u $Q/build_side.py \
  /scratch/jshen/data/pansoma_v2_tensors/$S/tensors/SNV $Q/$S.SNV.side.npz /scratch/jshen/data/pansoma_net_v2_runs/index_cache/$C

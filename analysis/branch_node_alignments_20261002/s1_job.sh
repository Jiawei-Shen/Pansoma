#!/bin/bash
# s1_branch_node_stats.py for one set: sbatch --array=0-5 s1_job.sh
#SBATCH -p general,gpu
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 2:00:00
#SBATCH -J branch_stats
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/branch_node_alignments_20261002/s1_%a_%A.out
set -euo pipefail
SETS=(HG008T_PacBio HG008T_ONT HG008T_Illumina COLO829T_fiberseq COLO829T_ONT COLO829T_Illumina)
PYTHONDONTWRITEBYTECODE=1 /usr/bin/time -f "%e s %M KB" /wanglab/jshen/anaconda3/bin/python \
  /scratch/jshen/Github/Pansoma/analysis/branch_node_alignments_20261002/s1_branch_node_stats.py "${SETS[$SLURM_ARRAY_TASK_ID]}"

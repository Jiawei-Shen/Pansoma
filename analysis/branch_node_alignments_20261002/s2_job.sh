#!/bin/bash
# s2_read_sample.py for one set: sbatch --array=0-5 s2_job.sh
#SBATCH -p general,gpu
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -J branch_reads
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/branch_node_alignments_20261002/s2_%a_%A.out
set -euo pipefail
SETS=(HG008T_PacBio HG008T_ONT HG008T_Illumina COLO829T_fiberseq COLO829T_ONT COLO829T_Illumina)
PYTHONDONTWRITEBYTECODE=1 /usr/bin/time -f "%e s %M KB" /wanglab/jshen/anaconda3/bin/python \
  /scratch/jshen/Github/Pansoma/analysis/branch_node_alignments_20261002/s2_read_sample.py "${SETS[$SLURM_ARRAY_TASK_ID]}"

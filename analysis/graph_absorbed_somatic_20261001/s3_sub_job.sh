#!/bin/bash
# s3_read_paths.py run, sub-parts: sbatch --array=0-(M-1) s3_sub_job.sh PLATFORM K N M  (array task id = sub-part j of M of part K of N)
#SBATCH -p general
#SBATCH -t 3:00:00
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
SUB="${SLURM_ARRAY_TASK_ID}/$4" /wanglab/jshen/anaconda3/bin/python \
  /scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/s3_read_paths.py run "$1" "$2" "$3"

#!/bin/bash
# s3_read_paths.py run: sbatch [--array=0-(N-1)] s3_job.sh PLATFORM N [LIMIT]  (array task id = part K of N)
#SBATCH -p general
#SBATCH -t 6:00:00
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
/wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/s3_read_paths.py \
  run "$1" "${SLURM_ARRAY_TASK_ID:-0}" "$2" ${3:-}

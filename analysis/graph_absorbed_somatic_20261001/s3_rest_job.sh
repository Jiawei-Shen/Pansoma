#!/bin/bash
# s3_read_paths.py rest: sbatch --array=0-(N-1) s3_rest_job.sh PLATFORM TAG N   (array task id = K of N over rest{TAG}.ids)
# make the list first: python s3_read_paths.py rest-list PLATFORM TAG
#SBATCH -p general,gpu
#SBATCH -t 4:00:00
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
/wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/s3_read_paths.py \
  rest "$1" "$2" "${SLURM_ARRAY_TASK_ID}" "$3"

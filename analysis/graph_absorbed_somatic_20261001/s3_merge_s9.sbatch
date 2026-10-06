#!/bin/bash
# merge the PacBio / ONT read paths written so far (s3_read_paths.py merge; safe while rest jobs still run: a locus counts
# only once its loci line is complete), then re-run s9_integrate.py
#SBATCH -p general,gpu
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH -t 2:00:00
#SBATCH -J s3merge_s9
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s3merge_s9.%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
A=/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001
for p in PacBio ONT; do /usr/bin/time -v /wanglab/jshen/anaconda3/bin/python $A/s3_read_paths.py merge $p 2>&1 | grep -v "^\s\s*[A-Z].*: [0-9]" ; done
/wanglab/jshen/anaconda3/bin/python $A/s9_integrate.py 16

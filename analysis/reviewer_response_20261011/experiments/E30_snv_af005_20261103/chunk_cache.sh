#!/bin/bash
# The zstd chunk cache (and index cache) of one _snvaf05 set, built ahead of the GPU job (chunks.py main; train --cache-dir
# R/index_cache, chunk cache R/index_cache/chunks); the training then finds every shard cached.
#SBATCH -p general
#SBATCH -c 16
#SBATCH --mem=12G
#SBATCH -t 12:00:00
set -euo pipefail
R=/scratch/jshen/data/pansoma_net_v2_runs
cd /scratch/jshen/Github/Pansoma/tmp/snv_af05_20261103/code_ml/machine_learning
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python -m pansoma_net_v2.chunks --tensors /scratch/jshen/data/pansoma_v2_tensors/$1/tensors \
  --kinds SNV --index-cache $R/index_cache --cache $R/index_cache/chunks --threads 16

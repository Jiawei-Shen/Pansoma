#!/bin/bash
# RIS INDEL ppr8 single-set checkpoint (ris_ppr8_20261026, 9919d4a recipe) on one set / chromosome, with this repo's
# pansoma_net_v2 and the local site-feature files (sha256 compared with the checkpoint's by predict).
# usage: sbatch predict.sh RUN SET CHROM
#SBATCH -p gpu
#SBATCH -w tequila
#SBATCH --gres=gpu:12gb:1
#SBATCH -c 8
#SBATCH --mem=11G
#SBATCH -t 4:00:00
set -euo pipefail
RUN=$1; SET=$2; CHROM=$3
R=/scratch/jshen/data/pansoma_net_v2_runs; D=/scratch/jshen/data/pansoma_v2_tensors; P=/wanglab/jshen/anaconda3/bin/python
echo "host $(hostname) job ${SLURM_JOB_ID:-} gpus ${CUDA_VISIBLE_DEVICES:-} repo $(git -C /scratch/jshen/Github/Pansoma rev-parse --short HEAD); $RUN $SET $CHROM"
cd /scratch/jshen/Github/Pansoma/machine_learning
$P -m pansoma_net_v2.predict --checkpoint $R/ris_ppr8_20261026/$RUN/best.pth --tensors $D/$SET/tensors --kinds INDEL \
  --chroms $CHROM --output $R/ris_ppr8_20261026/$RUN/test_${CHROM}_$SET --cache-dir $R/index_cache --num-workers 8 \
  --site-features /scratch/jshen/data/indel_gap_20261020/site_features 2>&1 | grep --line-buffered -v 'UserWarning\|warnings.warn'

#!/bin/bash
# E28 (2): predict one test set / chromosome with a chosen checkpoint of a selcheck rerun (best.pth = the original
# selection, best_excl_offpath.pth = the selection without the off-path validation negatives), as
# tmp/indel_gap_20261020/jobs/predict_set.sh does for best.pth, with the frozen code_final (77cee0e).
# usage: sbatch predict_ckpt.sh RUN(under pansoma_net_v2_runs) CKPT SET CHROM TESTDIR
#SBATCH -p gpu
#SBATCH -w tequila
#SBATCH --gres=gpu:12gb:1
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 4:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/rerun_20261101/indel_select/logs/%x_%j.out
set -euo pipefail
RUN=$1; CKPT=$2; SET=$3; CHROM=$4; TEST=$5
R=/scratch/jshen/data/pansoma_net_v2_runs; D=/scratch/jshen/data/pansoma_v2_tensors; P=/wanglab/jshen/anaconda3/bin/python
CODE=/scratch/jshen/Github/Pansoma/tmp/indel_gap_20261020/code_final/machine_learning
echo "host $(hostname) job ${SLURM_JOB_ID:-} gpus ${CUDA_VISIBLE_DEVICES:-} code $CODE; $RUN/$CKPT $SET $CHROM -> $TEST"
cd $CODE
$P -m pansoma_net_v2.predict --checkpoint $R/$RUN/$CKPT --tensors $D/$SET/tensors --kinds INDEL --chroms $CHROM \
  --output $R/$RUN/$TEST --cache-dir $R/index_cache --num-workers 8 2>&1 | grep --line-buffered -v 'UserWarning\|warnings.warn'

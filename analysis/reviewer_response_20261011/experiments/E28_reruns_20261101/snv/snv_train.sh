#!/bin/bash
#SBATCH -p gpu
#SBATCH -w tequila
#SBATCH --gres=gpu:h100:1
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -t 3-00:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/rerun_20261101/snv/logs/%x_%j.out
# Reruns of 2026-11-01 (user: run all), one H100 each, the SNV recipe of Table R3.1 (50.4 M, 20 epochs, pow:0.75,
# lr 1e-4, --norm train, global batch 4096 = 1 GPU x 4096 in 4 micro-batches; the RIS controls: 4 GPUs x 1024), then
# chr1 predictions; the VCF + rtg step is a separate CPU job (submit.sh).
#   MODE mask    R1.1: --mask-off-node (cells outside the candidate's node visit zeroed, q6 offnode.py), node-block split
#   MODE nomask  R1.1: the same setup without the mask (control)
#   MODE shared  R3.1: validation = chr20-22 for every reference (a split shared by GRCh38-only, v1.1 and v2.1)
# usage: sbatch -J NAME snv_train.sh NAME SET MODE SEED
set -euo pipefail
NAME=$1 SET=$2 MODE=$3 SEED=$4
W=/scratch/jshen/Github/Pansoma/tmp/rerun_20261101/snv
SIDE=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/side
R=/scratch/jshen/data/pansoma_net_v2_runs; D=/scratch/jshen/data/pansoma_v2_tensors; P=/wanglab/jshen/anaconda3/bin/python
case $MODE in
  mask)   ROOT=offnode_retrain_20261101; EXTRA="--mask-off-node $SIDE"; SPLIT="" ;;
  nomask) ROOT=offnode_retrain_20261101; EXTRA=""; SPLIT="" ;;
  shared) ROOT=sharedsplit_20261101; EXTRA=""; SPLIT="--val-chroms chr20 chr21 chr22" ;;
  *) echo "MODE mask|nomask|shared"; exit 1 ;;
esac
O=$R/$ROOT/$NAME
cd $W/code
echo "host $(hostname) job ${SLURM_JOB_ID:-} $NAME $SET $MODE seed $SEED $(date)"; nvidia-smi --query-gpu=index,name,memory.total --format=csv || true
if [ ! -s $O/best.pth ] || ! grep -q "^done; best epoch" $O/train.log 2>/dev/null; then
  $P -m pansoma_net_v2.train --tensors $D/$SET/tensors --kinds SNV --output $O --cache-dir $R/index_cache \
    --test-chroms chr1 $SPLIT \
    --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 4096 --micro-batches 4 --lr 1e-4 \
    --class-weights pow:0.75 --non-fraction 1.0 --block v2 \
    --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth \
    --norm train --seed $SEED --read-threads 8 --val-workers 3 $EXTRA 2>&1 \
    | grep --line-buffered -v 'UserWarning\|warnings.warn\|x = torch.from_numpy\|isnan'
fi
grep -q "^done; best epoch" $O/train.log
[ -s $O/test_chr1/$SET.tensors.SNV.predictions.ndjson.gz ] || $P -m pansoma_net_v2.predict --checkpoint $O/best.pth \
  --tensors $D/$SET/tensors --kinds SNV --chroms chr1 --output $O/test_chr1 --cache-dir $R/index_cache \
  --num-workers 6 --batch-size 256 $EXTRA
echo "done $(date)"

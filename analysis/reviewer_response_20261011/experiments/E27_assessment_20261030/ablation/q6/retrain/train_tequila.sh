#!/bin/bash
#SBATCH -p gpu
#SBATCH -w tequila
#SBATCH --gres=gpu:h100:2
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/retrain/logs/%x_%j.out
# NOT SUBMITTED unless the user chooses tequila for the retraining (else RIS, 4 x H100 as the controls).
# R1.1 cross-node context, retraining: the v2.1 (d46) SNV recipe of <S>_<P>_SNV_d46_7fc7396 (50.4 M, 20 epochs,
# pow:0.75, lr 1e-4, --norm train, seed SEED, global batch 4 x 1024 = here 2 GPUs x 2048 with 2 micro-batches) with
# MODE=mask (--mask-off-node) or MODE=nomask (a control on the same 2-GPU setup), then chr1 predict + jobs/vcf.sh.
# usage: sbatch -J q6_<NAME> train_tequila.sh NAME SET MODE SEED
set -euo pipefail
NAME=$1 SET=$2 MODE=$3 SEED=$4
Q=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6
R=/scratch/jshen/data/pansoma_net_v2_runs; D=/scratch/jshen/data/pansoma_v2_tensors; P=/wanglab/jshen/anaconda3/bin/python
OUTRUN=offnode_mask_20261030/$NAME; O=$R/$OUTRUN
EXTRA=""; [ $MODE = mask ] && EXTRA="--mask-off-node $Q/side"
cd $Q/code
echo "host $(hostname) job ${SLURM_JOB_ID:-} $NAME $SET $MODE seed $SEED $(date)"; nvidia-smi --query-gpu=index,name,memory.total --format=csv || true
if [ ! -s $O/best.pth ]; then
  $P -m torch.distributed.run --standalone --nproc_per_node=2 -m pansoma_net_v2.train --ddp --tensors $D/$SET/tensors --kinds SNV \
    --output $O --cache-dir $R/index_cache --test-chroms chr1 \
    --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 2048 --micro-batches 2 --lr 1e-4 --class-weights pow:0.75 \
    --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth \
    --norm train --seed $SEED --read-threads 8 --val-workers 3 $EXTRA
fi
[ -s $O/test_chr1/$SET.tensors.SNV.predictions.ndjson.gz ] || $P -m pansoma_net_v2.predict --checkpoint $O/best.pth \
  --tensors $D/$SET/tensors --kinds SNV --chroms chr1 --output $O/test_chr1 --cache-dir $R/index_cache --num-workers 6 --batch-size 256 $EXTRA
[ -s $O/vcf_chr1/eval_nobed/report.json ] || SET_NAME=$SET bash /scratch/jshen/Github/Pansoma/machine_learning/pansoma_net_v2/jobs/vcf.sh $OUTRUN SNV $O/vcf_chr1
echo "done $(date)"

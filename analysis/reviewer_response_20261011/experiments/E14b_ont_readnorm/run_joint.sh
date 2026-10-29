#!/bin/bash
# COLO829T ONT purity joint SNV model (P100 + P75 + P50 + P25; chr2-22 train, the same 5 % node blocks validate in every
# set, chr1 of each set tested), with the recipe of the RIS run COLO829_ONT_SNV_joint_7fc7396_gpu4 (ris/recipes.sh SNV,
# --read-threads 8 --val-workers 6) and its global batch of 4,096 on one H100 (--batch-size 4096 --micro-batches 4:
# gradient accumulation, the same update), --norm train (BQ z-scored with the training statistics, as 7fc7396 did).
# MODE=readnorm trains and predicts through ris/tools/readnorm.py (BQ centred on each read's own median); MODE=control
# runs the package unchanged, so the pair differs only in the BQ plane (same code, hardware, seed 0). Predictions: chr1 of
# each set with the off-reference rescue on (--offref-site90 10), as the RIS predictions of 2026-10-12.
# usage: sbatch [resources] -J NAME -o LOG run_joint.sh NAME control|readnorm
set -euo pipefail
NAME=$1; MODE=$2
R=/scratch/jshen/data/pansoma_net_v2_runs; OUT=$R/$NAME; P=/wanglab/jshen/anaconda3/bin/python
D=/scratch/jshen/data/pansoma_v2_tensors; REPO=/scratch/jshen/Github/Pansoma
SETS=($D/COLO829T_ONT/tensors $D/COLO829T_ONT_P75/tensors $D/COLO829T_ONT_P50/tensors $D/COLO829T_ONT_P25/tensors)
cd $REPO/machine_learning; export PYTHONPATH=$REPO/machine_learning
case $MODE in
  control)  T_ENTRY=(-m pansoma_net_v2.train) P_ENTRY=(-m pansoma_net_v2.predict);;
  readnorm) T_ENTRY=($REPO/ris/tools/readnorm.py train) P_ENTRY=($REPO/ris/tools/readnorm.py predict);;
  *) echo "MODE must be control or readnorm" >&2; exit 2;;
esac
mkdir -p $OUT
echo "host $(hostname) job ${SLURM_JOB_ID:-} gpus ${CUDA_VISIBLE_DEVICES:-} mode $MODE repo $(git -C $REPO rev-parse --short HEAD)"
nvidia-smi --query-gpu=index,name,memory.total --format=csv
nvidia-smi --query-gpu=timestamp,index,memory.used,utilization.gpu --format=csv,noheader -l 10 > $OUT/gpu_usage.csv &
echo "== train"
$P -u "${T_ENTRY[@]}" --tensors "${SETS[@]}" --kinds SNV --output $OUT --cache-dir $R/index_cache \
  --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 4096 --micro-batches 4 --lr 1e-4 \
  --class-weights pow:0.75 --non-fraction 1.0 --block v2 \
  --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth \
  --read-threads 8 --val-workers 6 --norm train 2>&1 \
  | grep --line-buffered -v 'UserWarning\|warnings.warn\|x = torch.from_numpy\|isnan'
echo "== test chr1 (rescue on)"
$P -u "${P_ENTRY[@]}" --checkpoint $OUT/best.pth --tensors "${SETS[@]}" --kinds SNV --chroms chr1 \
  --output $OUT/test_chr1 --cache-dir $R/index_cache --num-workers 8 --offref-site90 10
awk -F', ' '{gsub(/ MiB/, "", $3); if ($3 > m[$2]) m[$2] = $3} END {for (g in m) print "GPU " g ": " m[g] " MiB"}' $OUT/gpu_usage.csv
echo "done $(date)"

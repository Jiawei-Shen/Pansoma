#!/bin/bash
# Rerun 2026-11-01 (R1.8): one final INDEL release model retrained with the same recipe and seed by the code copy whose
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/rerun_20261101/indel_select/logs/%x_%j.out
# train.py also logs, each epoch, the selection score without the off-path validation negatives (val.excl_offpath).
# Original: tmp/indel_gap_20261020/final_20261023/train_final.sh (pre-registered INDEL batch).
# Writes best.pth (selection as in the original) and best_excl_offpath.pth (selection without the off-path negatives).
# usage: sbatch --gres gpu:<h100 or 24gb>:1 ... train_select.sh RUN PROTOCOL(release|unseen) PLATFORM SEED MICRO_BATCHES
set -euo pipefail
RUN=$1; PROTO=$2; p=$3; SEED=$4; MICRO=$5
R=/scratch/jshen/data/pansoma_net_v2_runs; D=/scratch/jshen/data/pansoma_v2_tensors
CODE=/scratch/jshen/Github/Pansoma/tmp/rerun_20261101/indel_select/code/machine_learning
P=/wanglab/jshen/anaconda3/bin/python
co=COLO829T_$p; [ $p = PacBio ] && co=COLO829T_fiberseq
tumors="HCC1395T HCC1954T H1437T H2009T"
[ $PROTO = release ] && tumors+=" HG008T COLO829T"
train=""; held=""
for t in $tumors; do
  s=$D/${t}_${p}_d46/tensors; [ $t = COLO829T ] && s=$D/${co}_d46/tensors
  train+=" $s"
  if [ $t = HCC1395T ]; then held+=" --set-test-chroms $s chr1 chr20"; else held+=" --set-test-chroms $s chr1"; fi
done
OUT=$R/selcheck_20261101/$RUN
mkdir -p $OUT
echo "host $(hostname) job ${SLURM_JOB_ID:-} gpus ${CUDA_VISIBLE_DEVICES:-} code $CODE (base $(cat $CODE/../BASE_COMMIT)); $RUN $PROTO $p seed $SEED micro $MICRO"
nvidia-smi --query-gpu=index,name,memory.total --format=csv
nvidia-smi --query-gpu=timestamp,index,memory.used,utilization.gpu --format=csv,noheader -l 10 > $OUT/gpu_usage.csv &
cd $CODE
$P -m pansoma_net_v2.train --tensors $train --kinds INDEL --output $OUT --cache-dir $R/index_cache \
  --test-chroms $held \
  --depths 2 2 6 2 --dims 64 128 256 512 --epochs 12 --batch-size 1024 --lr 1e-4 --weight-decay 0.05 \
  --warmup-epochs 1 --class-weights sqrt --non-fraction 1.0 --block v2 \
  --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth \
  --non-reasons repeat_indel_without_truth --non-reasons-in-bed outside_confident_region \
  --site-features /scratch/jshen/data/indel_gap_20261020/site_features --best-by postpon \
  --seed $SEED --micro-batches $MICRO 2>&1 | grep --line-buffered -v 'UserWarning\|warnings.warn\|x = torch.from_numpy\|isnan'
grep -q "^done; best epoch" $OUT/train.log

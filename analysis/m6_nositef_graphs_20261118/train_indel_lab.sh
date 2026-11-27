#!/bin/bash
# One six-tumor INDEL model without the population site features as a model input on HPRC v1.1 (d9) or v2.1 (d46), on
# the lab cluster (user 2026-11-24; the INDEL tensors of the CASTLE lines are not on RIS storage2). The recipe and code
# of the GRCh38-only model it pairs with (m4m6_graphs_20261110/multicancer6_<P>_INDEL_linear_nositef,
# tmp/m4m6_graphs_20261110/train.sh ... nositef): the final INDEL recipe (code_final, 8.7 M, 12 epochs, batch 1024,
# lr 1e-4, sqrt, wd 0.05, warmup 1, rule (a) negatives, --best-by postpon, --micro-batches 4, seed 0) with
# --pon-features instead of --site-features: the feature files give only the PoN tags for the post-PoN epoch selection
# and the call rule, the five features are not fed to the model and predict reads no feature file.
# Sets = those of the model with the features on the same graph:
#   d9   m4m6_graphs_20261110/multicancer6_<P>_INDEL_d9 (tmp/m4m6_graphs_20261110/sets.sh)
#   d46  indel_gap_20261020/final_20261023 release_s0 (<tumor>_<P>_d46, COLO829T PacBio = COLO829T_fiberseq_d46)
# Held out of training and validation: chr1 of every tumor and HCC1395T chr20. Then predicts chr1 of HCC1395T, HG008T
# and COLO829T and chr20 of HCC1395T with best.pth. A job cut by its time limit continues from last.pth when resubmitted.
# usage: sbatch --gres=gpu:24gb:1 | --gres=gpu:h100:1 [-t ...] train_indel_lab.sh PLAT(PacBio|ONT|Illumina) GRAPH(d9|d46)
#SBATCH -p gpu
#SBATCH -w tequila
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 2-00:00:00
set -euo pipefail
PLAT=$1 G=$2; K=INDEL
R=/scratch/jshen/data/pansoma_net_v2_runs; OUT=$R/m6_nositef_graphs_20261118/multicancer6_${PLAT}_INDEL_${G}_nositef
D=/scratch/jshen/data/pansoma_v2_tensors; P=/wanglab/jshen/anaconda3/bin/python
CODE=/scratch/jshen/Github/Pansoma/tmp/indel_gap_20261020/code_final/machine_learning
source /scratch/jshen/Github/Pansoma/tmp/m4m6_graphs_20261110/sets.sh
case $G in
  d9)  SF=/scratch/jshen/data/m4m6_graphs_20261110/indel_site_features; name() { set_name INDEL $PLAT d9 $1; } ;;
  d46) SF=/scratch/jshen/data/indel_gap_20261020/site_features
       name() { local p=$PLAT; [ $1 = COLO829T ] && [ $PLAT = PacBio ] && p=fiberseq; echo ${1}_${p}_d46; } ;;
  *) echo "GRAPH must be d9 or d46"; exit 1 ;;
esac
T=(); HELD=""
for t in $TUMORS6; do
  s=$(name $t); T+=($D/$s/tensors)
  [ -s $D/$s/tensors/$K/labels.manifest.json ] || { echo "no labels $s"; exit 1; }
  [ -s $SF/$s.$K.npz ] || { echo "no site features $SF/$s.$K.npz"; exit 1; }
  if [ $t = HCC1395T ]; then HELD+=" --set-test-chroms $D/$s/tensors chr1 chr20"; else HELD+=" --set-test-chroms $D/$s/tensors chr1"; fi
done
cd $CODE
mkdir -p $OUT
echo "host $(hostname) job ${SLURM_JOB_ID:-} gpus ${CUDA_VISIBLE_DEVICES:-} code $CODE $K $PLAT $G nositef: ${T[*]} pon features $SF $(date)"
nvidia-smi --query-gpu=index,name,memory.total --format=csv
RES=""; [ -s $OUT/last.pth ] && { RES="--resume $OUT/last.pth"; echo "resuming from $OUT/last.pth"; }
if ! grep -q "^done; best epoch" $OUT/train.log 2>/dev/null; then
  $P -m pansoma_net_v2.train --tensors "${T[@]}" --kinds $K --output $OUT --cache-dir $R/index_cache --test-chroms $HELD \
    --block v2 --lr 1e-4 --batch-size 1024 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth --seed 0 \
    --depths 2 2 6 2 --dims 64 128 256 512 --epochs 12 --weight-decay 0.05 --warmup-epochs 1 --class-weights sqrt \
    --non-fraction 1.0 --non-reasons repeat_indel_without_truth --non-reasons-in-bed outside_confident_region \
    --pon-features $SF --best-by postpon --micro-batches 4 $RES 2>&1 | grep --line-buffered -v 'UserWarning\|warnings.warn\|x = torch.from_numpy\|isnan'
fi
grep -q "^done; best epoch" $OUT/train.log
HCC=$D/$(name HCC1395T)/tensors
$P -m pansoma_net_v2.predict --checkpoint $OUT/best.pth --kinds $K --chroms chr1 --output $OUT/test_chr1 \
  --tensors $HCC $D/$(name HG008T)/tensors $D/$(name COLO829T)/tensors \
  --cache-dir $R/index_cache --num-workers 4 2>&1 | grep --line-buffered -v 'UserWarning\|warnings.warn'
$P -m pansoma_net_v2.predict --checkpoint $OUT/best.pth --kinds $K --chroms chr20 --output $OUT/test_chr20 --tensors $HCC \
  --cache-dir $R/index_cache --num-workers 4 2>&1 | grep --line-buffered -v 'UserWarning\|warnings.warn'
echo "done $(date)"

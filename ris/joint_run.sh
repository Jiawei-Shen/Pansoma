#!/bin/bash
# Train one model on several sets (e.g. a sample's purity dilutions; chr1 left out of every set, 5 % of the chr2-22 node
# blocks validate: the same blocks in every set), then predict chr1 of each set with best.pth into
# <run>/test_chr1_<SET>/. With several GPUs (submit.sh -g N) training runs under torchrun --ddp.
# Usage: joint_run.sh NAME KIND SET[,SET...] [train options]   (joint_pipeline.sh submits it with the recipe)
set -euo pipefail
source "$(dirname "$0")/env.sh"
NAME=$1 KIND=$2 SETS=$3; shift 3
OUT=$PANSOMA_RUNS/$NAME
T=(); for s in ${SETS//,/ }; do T+=("$PANSOMA_DATA/pansoma_v2_tensors/$s/tensors"); done
NGPU=$(python -c "import torch; print(torch.cuda.device_count())")
echo "host $(hostname)  job ${LSB_JOBID:-}  GPUs $NGPU (${CUDA_VISIBLE_DEVICES:-})  $(date)  sets $SETS"
ref=$(sed -n "s/^ref: //p" "$REPO/.git/HEAD")  # no git in the image
echo "repo $ref $(cat "$REPO/.git/$ref" 2>/dev/null || grep " $ref$" "$REPO/.git/packed-refs" | cut -c1-12)"
cd "$REPO/machine_learning"
train=(-m pansoma_net_v2.train --tensors "${T[@]}" --kinds "$KIND" --output "$OUT" --cache-dir "$PANSOMA_RUNS/index_cache" "$@")
if [ "$NGPU" -gt 1 ]; then
  torchrun --standalone --nproc_per_node="$NGPU" "${train[@]}" --ddp
else
  python -u "${train[@]}"
fi
for s in ${SETS//,/ }; do
  python -u -m pansoma_net_v2.predict --checkpoint "$OUT/best.pth" --tensors "$PANSOMA_DATA/pansoma_v2_tensors/$s/tensors" \
    --kinds "$KIND" --chroms chr1 --cache-dir "$PANSOMA_RUNS/index_cache" --output "$OUT/test_chr1_$s" \
    --num-workers "${PREDICT_WORKERS:-4}"
done
echo "done $(date)"

#!/bin/bash
# ris/joint_run.sh with the training ended after epoch STOP: train one model on several sets through
# ris/tools/save_epochs.py --stop-after-epoch STOP (the recipe's --epochs and so its learning-rate schedule unchanged:
# epochs 1-STOP are those of the full run; every epoch's checkpoint is also kept in <run>/epochs/), then predict chr1
# of each set with best.pth (the best of epochs 1-STOP by the run's own selection) into <run>/test_chr1_<SET>/.
# The user (2026-11-12): 10 epochs are enough for the COLO829T dilution models.
# Usage: joint_run_stop.sh STOP NAME KIND SET[,SET...] [train options]. A GPU job:
#   ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 4 -n 16 -m 384 pnv2_NAME \
#       ../analysis/colo829t_dilution_newpipe_20261112/joint_run_stop.sh 10 NAME SNV SET,SET... <recipe>
set -euo pipefail
source "$(dirname "$0")/../../ris/env.sh"
STOP=$1 NAME=$2 KIND=$3 SETS=$4; shift 4
OUT=$PANSOMA_RUNS/$NAME
T=(); for s in ${SETS//,/ }; do T+=("$PANSOMA_DATA/pansoma_v2_tensors/$s/tensors"); done
NGPU=$(python -c "import torch; print(torch.cuda.device_count())")
echo "host $(hostname)  job ${LSB_JOBID:-}  GPUs $NGPU (${CUDA_VISIBLE_DEVICES:-})  $(date)  sets $SETS  stop after epoch $STOP"
ref=$(sed -n "s/^ref: //p" "$REPO/.git/HEAD")  # no git in the image
echo "repo $ref $(cat "$REPO/.git/$ref" 2>/dev/null || grep " $ref$" "$REPO/.git/packed-refs" | cut -c1-12)"
cd "$REPO/machine_learning"
export PYTHONPATH="$REPO/machine_learning${PYTHONPATH:+:$PYTHONPATH}"
train=("$REPO/ris/tools/save_epochs.py" --stop-after-epoch "$STOP" --tensors "${T[@]}" --kinds "$KIND" --output "$OUT"
       --cache-dir "$PANSOMA_RUNS/index_cache" "$@")
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

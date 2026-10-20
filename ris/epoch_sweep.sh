#!/bin/bash
# Train like run.sh (no predict; with more than one GPU under torchrun with --ddp, as run.sh), keeping every epoch's
# checkpoint in <run>/epochs/epoch_NN.pth (tools/save_epochs.py; --stop-after-epoch N ends it once epoch N is saved).
# Usage: epoch_sweep.sh NAME SET KIND [train options]. Each epoch is then predicted and scored
# as its own run: <run>_eNN/ with args.json and best.pth linked to <run>/args.json and <run>/epochs/epoch_NN.pth,
# through predict.sh and vcf.sh (see analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/
# submit_epoch_sweep.sh).
set -euo pipefail
source "$(dirname "$0")/env.sh"
NAME=$1 SET=$2 KIND=$3; shift 3
OUT=$PANSOMA_RUNS/$NAME
NGPU=$(python -c "import torch; print(torch.cuda.device_count())")
echo "host $(hostname)  job ${LSB_JOBID:-}  GPUs $NGPU (${CUDA_VISIBLE_DEVICES:-})  $(date)"
ref=$(sed -n "s/^ref: //p" "$REPO/.git/HEAD")
echo "repo $ref $(cat "$REPO/.git/$ref" 2>/dev/null || grep " $ref$" "$REPO/.git/packed-refs" | cut -c1-12)"
cd "$REPO/machine_learning"
export PYTHONPATH="$REPO/machine_learning${PYTHONPATH:+:$PYTHONPATH}"
train=("$REPO/ris/tools/save_epochs.py" --tensors "$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors" --kinds "$KIND"
       --output "$OUT" --cache-dir "$PANSOMA_RUNS/index_cache" "$@")
if [ "$NGPU" -gt 1 ]; then
  torchrun --standalone --nproc_per_node="$NGPU" "${train[@]}" --ddp
else
  python -u "${train[@]}"
fi
ls -la "$OUT/epochs"
echo "done $(date)"

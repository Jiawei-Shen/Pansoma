#!/bin/bash
# Train on a set (chr1 left out), then predict chr1 with the best checkpoint (predict.sh). Usage: run.sh NAME SET KIND [train options]
# The train options are passed as given (the SNV or INDEL recipe: machine_learning/pansoma_net_v2/README.md,
# "SNV recipe" / "INDEL recipe"; snv_pipeline.sh passes them).
# With more than one GPU (submit.sh -g N) training runs under torchrun with --ddp.
#   ris/submit.sh -n 1 -m 64 pnv2_HG008_Illumina_SNV_test run.sh HG008_Illumina_SNV_test HG008T_Illumina SNV \
#       --train-chroms chr21 chr22 --epochs 2 --read-threads 4 --num-workers 4 --val-workers 4
set -euo pipefail
source "$(dirname "$0")/env.sh"
NAME=$1 SET=$2 KIND=$3; shift 3
T=$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors
OUT=$PANSOMA_RUNS/$NAME
NGPU=$(python -c "import torch; print(torch.cuda.device_count())")
echo "host $(hostname)  job ${LSB_JOBID:-}  GPUs $NGPU (${CUDA_VISIBLE_DEVICES:-})  $(date)"
ref=$(sed -n "s/^ref: //p" "$REPO/.git/HEAD")  # no git in the image
echo "repo $ref $(cat "$REPO/.git/$ref" 2>/dev/null || grep " $ref$" "$REPO/.git/packed-refs" | cut -c1-12)"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader -i "${CUDA_VISIBLE_DEVICES:-0}"
cd "$REPO/machine_learning"
train=(-m pansoma_net_v2.train --tensors "$T" --kinds "$KIND" --output "$OUT" --cache-dir "$PANSOMA_RUNS/index_cache" "$@")
if [ "$NGPU" -gt 1 ]; then
  torchrun --standalone --nproc_per_node="$NGPU" "${train[@]}" --ddp
else
  python -u "${train[@]}"
fi
bash "$REPO/ris/predict.sh" "$NAME" "$SET" "$KIND"
echo "done $(date)"

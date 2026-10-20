#!/bin/bash
# Debugging aid: train one epoch with the code of another checkout (e.g. a git worktree of an older commit), no
# predict, to compare its epoch-1 loss and validation with an earlier run (2026-10-19: COLO829_Illumina_SNV_h90 of
# 9d6bd8d against today's reruns). The first epoch is the learning-rate warm-up, which does not depend on --epochs,
# so epoch 1 of a 1-epoch run equals epoch 1 of a 20-epoch run. A GPU job:
#   ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n 7 -m 192 NAME tools/train_epoch1.sh CODE NAME SET KIND [train options]
# CODE: the checkout whose machine_learning/ is run (this repository's ris/env.sh sets the environment).
set -euo pipefail
source "$(dirname "$0")/../env.sh"
CODE=$1 NAME=$2 SET=$3 KIND=$4; shift 4
OUT=$PANSOMA_RUNS/$NAME
echo "host $(hostname)  job ${LSB_JOBID:-}  code $CODE ($(cat "$CODE/.git" 2>/dev/null | head -1))  $(date)"
cd "$CODE/machine_learning"
export PYTHONPATH="$CODE/machine_learning${PYTHONPATH:+:$PYTHONPATH}"
python -u -m pansoma_net_v2.train --tensors "$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors" --kinds "$KIND" \
    --output "$OUT" "$@" --epochs 1   # last: overrides a recipe's --epochs
echo "done $(date)"

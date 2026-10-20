#!/bin/bash
# Debugging aid: train with the code of another checkout (e.g. a git worktree of an older commit), as run.sh does
# but without predict (predict.sh / vcf.sh of this repository score the result: predict is unchanged). A GPU job:
#   ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n 7 -m 192 NAME tools/train_at.sh CODE NAME SET KIND [train options]
# CODE: the checkout whose machine_learning/ is run (this repository's ris/env.sh sets the environment).
# 2026-10-19: COLO829_Illumina_SNV_h90 (9d6bd8d) is reproduced only by its own code: from f2a0e4e the encoder's
# z-score is written with broadcast statistics, which torch.compile turns into other kernels (last-bit differences).
set -euo pipefail
source "$(dirname "$0")/../env.sh"
CODE=$1 NAME=$2 SET=$3 KIND=$4; shift 4
OUT=$PANSOMA_RUNS/$NAME
echo "host $(hostname)  job ${LSB_JOBID:-}  code $CODE  $(date)"
cd "$CODE/machine_learning"
export PYTHONPATH="$CODE/machine_learning${PYTHONPATH:+:$PYTHONPATH}"
python -u -m pansoma_net_v2.train --tensors "$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors" --kinds "$KIND" \
    --output "$OUT" "$@"
echo "done $(date)"

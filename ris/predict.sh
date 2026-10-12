#!/bin/bash
# Predict chromosomes of a set with a run's best.pth (default chr1). Usage: predict.sh NAME SET KIND [predict options]
# Every DataLoader worker holds a copy of the dataset index (several GB for a 5 M-tensor set): with 12 workers
# chr1 of HG008T_Illumina SNV passed 64 GB, so the default is 4 (PREDICT_WORKERS).
#   ris/submit.sh -n 1 -m 64 pnv2_predict_SNV_test predict.sh HG008_Illumina_SNV_test HG008T_Illumina SNV
set -euo pipefail
source "$(dirname "$0")/env.sh"
NAME=$1 SET=$2 KIND=$3; shift 3
OUT=$PANSOMA_RUNS/$NAME
cd "$REPO/machine_learning"
time python -u -m pansoma_net_v2.predict --checkpoint "$OUT/best.pth" --tensors "$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors" \
    --kinds "$KIND" --chroms chr1 --cache-dir "$PANSOMA_RUNS/index_cache" --output "$OUT/test_chr1" \
    --num-workers "${PREDICT_WORKERS:-4}" "$@"

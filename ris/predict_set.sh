#!/bin/bash
# Predict chr1 of another set with a run's best.pth (e.g. a cross-tumor test) into <run>/test_chr1_<SET>/, the
# directory vcf.sh --set SET scores. Usage: predict_set.sh NAME SET KIND [predict options]. A GPU job:
#   ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -n 4 -m 64 pnv2_xpred predict_set.sh \
#       HG008_Illumina_INDEL_d46_ppr8 COLO829T_Illumina_d46 INDEL
set -euo pipefail
source "$(dirname "$0")/env.sh"
NAME=$1 SET=$2 KIND=$3; shift 3
OUT=$PANSOMA_RUNS/$NAME
echo "host $(hostname)  job ${LSB_JOBID:-}  $(date)  $NAME -> $SET"
cd "$REPO/machine_learning"
time python -u -m pansoma_net_v2.predict --checkpoint "$OUT/best.pth" --tensors "$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors" \
    --kinds "$KIND" --chroms chr1 --cache-dir "$PANSOMA_RUNS/index_cache" --output "$OUT/test_chr1_$SET" \
    --num-workers "${PREDICT_WORKERS:-4}" "$@"

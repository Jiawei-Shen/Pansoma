#!/bin/bash
# Submit one experiment (login node): exp_train.sh on an H100 (epigenome), then exp_calls.sh per set on CPU after it.
#   ris/exp_submit.sh NAME SET[,SET...] [--focal G] [train options]
# Train options default to the INDEL recipe with chr2 held out (--val-chroms chr2); given options are appended.
# SLOTS / MEM (default 5 / 192) for the GPU job.
set -euo pipefail
NAME=$1 SETS=$2; shift 2
HERE=$(cd "$(dirname "$0")" && pwd)
pre=(); if [ "${1:-}" = --focal ]; then pre=(--focal "$2"); shift 2; fi
RECIPE=(--depths 2 2 6 2 --dims 64 128 256 512 --epochs 12 --batch-size 1024 --lr 4e-4 --class-weights sqrt
        --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
        --val-chroms chr2 --read-threads 8 --val-workers 4)
"$HERE/submit.sh" -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -n "${SLOTS:-5}" -m "${MEM:-192}" "$NAME" \
    exp_train.sh "$NAME" INDEL "$SETS" "${pre[@]}" "${RECIPE[@]}" "$@"
for s in ${SETS//,/ }; do
  "$HERE/submit.sh" -g 0 -n 4 -m 32 -w "done($NAME)" "${NAME}_calls_$s" exp_calls.sh "$NAME" INDEL "$s"
done

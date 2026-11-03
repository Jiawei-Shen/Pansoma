#!/bin/bash
# One comparator call set scored on chr1 with Pansoma's protocol (score_comparator_ris.py), in the container. A CPU job:
#   ris/submit.sh -g 0 -n 4 -m 16 NAME ../analysis/indel_ppr8_9919d4a_20261026/comparators/score.sh \
#       CALLER SET KIND OUTPUT VCF [VCF ...]
set -euo pipefail
source "$(dirname "$0")/../../../ris/env.sh"
export PATH=$(dirname "$BCFTOOLS"):$PATH
CALLER=$1 SET=$2 KIND=$3 OUT=$4; shift 4
vcfs=(); for v in "$@"; do vcfs+=(--vcf "$v"); done
python -u "$(dirname "$0")/score_comparator_ris.py" "${vcfs[@]}" --set "$SET" --kind "$KIND" --caller "$CALLER" \
    --output "$OUT" --threads "${LSB_DJOB_NUMPROC:-4}"

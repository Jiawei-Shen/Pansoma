#!/usr/bin/env bash
set -euo pipefail

# FASTQ -> GAM with vg giraffe.
#
# Required variables:
#   GBZ, MIN_INDEX, DIST_INDEX, FASTQ1, OUT_GAM
#
# Optional variables:
#   FASTQ2      paired-end FASTQ
#   READ_TYPE   illumina, hifi, or ont/r10
#   THREADS     default: 12

: "${GBZ:?Set GBZ to the graph .gbz path}"
: "${MIN_INDEX:?Set MIN_INDEX to the .min index path}"
: "${DIST_INDEX:?Set DIST_INDEX to the .dist index path}"
: "${FASTQ1:?Set FASTQ1 to read 1 FASTQ path}"
: "${OUT_GAM:?Set OUT_GAM to output .gam path}"

THREADS="${THREADS:-12}"
READ_TYPE="${READ_TYPE:-illumina}"

mkdir -p "$(dirname "${OUT_GAM}")"

args=(
  vg giraffe
  -Z "${GBZ}"
  -m "${MIN_INDEX}"
  -d "${DIST_INDEX}"
  -f "${FASTQ1}"
  -t "${THREADS}"
  -o gam
)

if [[ -n "${FASTQ2:-}" ]]; then
  args+=( -f "${FASTQ2}" )
fi

case "${READ_TYPE}" in
  illumina)
    ;;
  hifi|pacbio)
    args+=( -b hifi )
    ;;
  ont|r10)
    args+=( -b r10 )
    ;;
  *)
    echo "ERROR: READ_TYPE must be illumina, hifi, pacbio, ont, or r10" >&2
    exit 2
    ;;
esac

"${args[@]}" > "${OUT_GAM}"


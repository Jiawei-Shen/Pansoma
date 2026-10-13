#!/bin/bash
# Index cache and zstd chunk cache of a set's training chromosomes, on CPU, ahead of a GPU run
# (train builds what is missing itself, on the GPU node). Usage: prebuild.sh SET KIND [chunks options]
#   ris/submit.sh -g 0 -n 8 -m 32 pnv2_prebuild prebuild.sh HG008T_Illumina SNV            # chr2-22
#   ris/submit.sh -g 0 -n 8 -m 32 pnv2_prebuild prebuild.sh HG008T_Illumina SNV --chroms chr21 chr22
set -euo pipefail
source "$(dirname "$0")/env.sh"
SET=$1 KIND=$2; shift 2
cd "$REPO/machine_learning"
time python -u -m pansoma_net_v2.chunks --tensors "$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors" --kinds "$KIND" \
    --index-cache "$PANSOMA_RUNS/index_cache" --cache "$PANSOMA_RUNS/index_cache/chunks" \
    --threads "${LSB_DJOB_NUMPROC:-8}" "$@"
du -sh "$PANSOMA_RUNS/index_cache/chunks" 2>/dev/null || true  # other builds may be writing temporaries there

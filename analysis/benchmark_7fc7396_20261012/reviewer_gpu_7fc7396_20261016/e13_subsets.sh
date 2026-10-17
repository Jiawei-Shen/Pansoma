#!/bin/bash
# CPU job (ris/submit.sh -g 0 ... ../analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/e13_subsets.sh):
# the E13 exclude files of HG008T_PacBio_d46 SNV (10 / 25 / 50 / 75 % of the training blocks, 2 replicates) into
# $PANSOMA_RUNS/e13_subsets/, as E21 README section 6.
set -euo pipefail
source "$(dirname "$0")/../../../ris/env.sh"
python -u "$(dirname "$0")/e13_subsets_ris.py" HG008T_PacBio_d46 SNV "$PANSOMA_RUNS/e13_subsets" 0.10 0.25 0.50 0.75 \
    --replicates 2

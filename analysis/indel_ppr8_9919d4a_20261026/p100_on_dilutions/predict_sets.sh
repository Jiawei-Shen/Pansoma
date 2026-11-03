#!/bin/bash
# A run's best.pth on chr1 of several sets (ris/predict_set.sh for each) in one GPU job.
#   usage: predict_sets.sh NAME KIND SET [SET ...]   (submit.sh)
set -euo pipefail
NAME=$1 KIND=$2; shift 2
for SET in "$@"; do bash "$(dirname "$0")/../../../ris/predict_set.sh" "$NAME" "$SET" "$KIND"; done

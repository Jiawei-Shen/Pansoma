#!/bin/bash
# CPU job: seed ensemble of exp_train.sh runs for one set: fuse their chr2 and chr1 predictions (tools/fuse_preds.py)
# into $PANSOMA_RUNS/NAME/pred_<SET>_<chrom>/, ready for exp_calls.sh NAME INDEL SET.
# Usage: ens_fuse.sh NAME SET RUN [RUN ...]
set -euo pipefail
source "$(dirname "$0")/env.sh"
NAME=$1 SET=$2; shift 2
for c in chr2 chr1; do
  m=(); for r in "$@"; do m+=("$PANSOMA_RUNS/$r/pred_${SET}_$c"); done
  python -u "$REPO/ris/tools/fuse_preds.py" "$PANSOMA_RUNS/$NAME/pred_${SET}_$c" "${m[@]}"
done

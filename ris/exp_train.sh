#!/bin/bash
# Experiment GPU job: train one model on one or more sets, then predict the tuning chromosome (chr2) and the test
# chromosome (chr1) of each set. Usage: exp_train.sh NAME KIND SET[,SET...] [train options]
# The train options must hold chr2 out of training (--val-chroms chr2). A first option --focal G trains with
# tools/train_focal.py (focal loss, gamma G).
# Predictions: $PANSOMA_RUNS/NAME/pred_<SET>_<chrom>/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
NAME=$1 KIND=$2 SETS=$3; shift 3
FOCAL=""; if [ "${1:-}" = --focal ]; then FOCAL=$2; shift 2; fi
OUT=$PANSOMA_RUNS/$NAME
T=(); for s in ${SETS//,/ }; do T+=("$PANSOMA_DATA/pansoma_v2_tensors/$s/tensors"); done
echo "host $(hostname)  job ${LSB_JOBID:-}  $(date)  sets ${SETS}  focal ${FOCAL:-no}"
ref=$(sed -n "s/^ref: //p" "$REPO/.git/HEAD"); echo "repo $ref $(cat "$REPO/.git/$ref" 2>/dev/null)"
cd "$REPO/machine_learning"
export PYTHONPATH="$REPO/machine_learning:$REPO/scripts:$REPO/ris/tools${PYTHONPATH:+:$PYTHONPATH}"
if [ -n "$FOCAL" ]; then export FOCAL_GAMMA=$FOCAL; entry=("$REPO/ris/tools/train_focal.py"); else entry=(-m pansoma_net_v2.train); fi
python -u "${entry[@]}" --tensors "${T[@]}" --kinds "$KIND" --output "$OUT" --cache-dir "$PANSOMA_RUNS/index_cache" "$@"
for s in ${SETS//,/ }; do
  for c in chr2 chr1; do
    python -u -m pansoma_net_v2.predict --checkpoint "$OUT/best.pth" --kinds "$KIND" --chroms "$c" \
      --tensors "$PANSOMA_DATA/pansoma_v2_tensors/$s/tensors" --cache-dir "$PANSOMA_RUNS/index_cache" \
      --output "$OUT/pred_${s}_$c" --num-workers "${PREDICT_WORKERS:-4}"
  done
done
echo "done $(date)"

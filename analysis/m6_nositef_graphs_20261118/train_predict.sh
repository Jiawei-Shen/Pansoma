#!/bin/bash
# One six-tumor model: train on SETS with some chromosomes of the test tumors left out of training and validation
# (HELD; the other sets train on every chromosome, 5 % of the node blocks validate), then predict every held
# chromosome of every held set with best.pth into <run>/test_<chrom>_<SET>/ (test_chr1_<SET> is what
# ris/vcf.sh --set SET scores). With several GPUs (submit.sh -g N) training runs under torchrun --ddp.
# Usage: train_predict.sh NAME KIND SET[,SET...] HELD [train options]      HELD = SET:chrom[+chrom...][,SET:...]
# A GPU job (submit.sh of this directory submits it with the recipe):
#   ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -n 8 -m 256 pnv2_NAME \
#       ../analysis/m6_nositef_graphs_20261118/train_predict.sh NAME SNV SET,SET... SET:chr1+chr20,SET:chr1 <recipe>
set -euo pipefail
source "$(dirname "$0")/../../ris/env.sh"
NAME=$1 KIND=$2 SETS=$3 HELD=$4; shift 4
OUT=$PANSOMA_RUNS/$NAME
D=$PANSOMA_DATA/pansoma_v2_tensors
T=(); for s in ${SETS//,/ }; do T+=("$D/$s/tensors"); done
H=(); for h in ${HELD//,/ }; do c=${h#*:}; H+=(--set-test-chroms "$D/${h%%:*}/tensors" ${c//+/ }); done
NGPU=$(python -c "import torch; print(torch.cuda.device_count())")
echo "host $(hostname)  job ${LSB_JOBID:-}  GPUs $NGPU (${CUDA_VISIBLE_DEVICES:-})  $(date)  sets $SETS  held $HELD"
ref=$(sed -n "s/^ref: //p" "$REPO/.git/HEAD")  # no git in the image
echo "repo $ref $(cat "$REPO/.git/$ref" 2>/dev/null || grep " $ref$" "$REPO/.git/packed-refs" | cut -c1-12)"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader -i "${CUDA_VISIBLE_DEVICES:-0}"
mkdir -p "$OUT"
cd "$REPO/machine_learning"
if ! grep -q "^done; best epoch" "$OUT/train.log" 2>/dev/null; then
  train=(-m pansoma_net_v2.train --tensors "${T[@]}" --kinds "$KIND" --output "$OUT" --cache-dir "$PANSOMA_RUNS/index_cache"
         --test-chroms "${H[@]}" "$@")
  if [ "$NGPU" -gt 1 ]; then
    torchrun --standalone --nproc_per_node="$NGPU" "${train[@]}" --ddp
  else
    python -u "${train[@]}"
  fi
fi
for h in ${HELD//,/ }; do
  s=${h%%:*} c=${h#*:}
  for chrom in ${c//+/ }; do
    python -u -m pansoma_net_v2.predict --checkpoint "$OUT/best.pth" --tensors "$D/$s/tensors" --kinds "$KIND" \
      --chroms "$chrom" --cache-dir "$PANSOMA_RUNS/index_cache" --output "$OUT/test_${chrom}_$s" \
      --num-workers "${PREDICT_WORKERS:-4}"
  done
done
echo "done $(date)"

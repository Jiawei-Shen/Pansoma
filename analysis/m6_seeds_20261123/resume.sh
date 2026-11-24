#!/bin/bash
# GPU job: move a running seed model off a throttled GPU of exec-419 (GPUs 0, 1 and 6 run at 86 C and half speed,
# 2026-11-23). Started by LSF on some GPU; on a throttled one (HOT) it leaves the running job alone, submits itself
# again to start when job RETRY ends (a job on a healthy GPU) and exits 3; on a healthy one it kills the running job
# OLD, waits for it to end, and continues its training from <run>/last.pth (train.py --resume: model, optimizer,
# scheduler, epoch and best epoch restored; checkpoints are written atomically, so the kill cannot corrupt them),
# then predicts as train_predict.sh does. LSF commands work inside the container.
# Usage (ris/submit.sh, from the repository root; job name <old name>r):
#   ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n 5 -m 192 pnv2_m6_Illumina_SNV_newpipe_s1r \
#     ../analysis/m6_seeds_20261123/resume.sh OLD HOT RETRY ATTEMPT <train_predict.sh arguments of OLD>
set -euo pipefail
OLD=$1 HOT=$2 RETRY=$3 ATTEMPT=$4; shift 4
source "$(dirname "$0")/../../ris/env.sh"
OUT=$PANSOMA_RUNS/$1
gpu=${CUDA_VISIBLE_DEVICES:-?}
echo "resume.sh attempt $ATTEMPT  job ${LSB_JOBID:-}  ${LSB_HOSTS:-}  GPU $gpu  $(date)  replaces $OLD"
nvidia-smi --query-gpu=index,temperature.gpu,power.draw --format=csv,noheader -i "$gpu" || true
if [[ ",$HOT," == *",$gpu,"* ]]; then
  echo "GPU $gpu is throttled; job $OLD left running"
  if [ "$ATTEMPT" -lt 3 ]; then
    "$REPO/ris/submit.sh" -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n 5 -m 192 -w "ended($RETRY)" \
      "$LSB_JOBNAME" ../analysis/m6_seeds_20261123/resume.sh "$OLD" "$HOT" "$RETRY" $((ATTEMPT + 1)) "$@"
  fi
  exit 3
fi
state=$(bjobs -noheader -o stat "$OLD" 2>/dev/null || true)
if [ "$state" = DONE ]; then echo "job $OLD already done"; exit 0; fi
if [ "$state" = RUN ] || [ "$state" = PEND ]; then
  bkill "$OLD"
  while state=$(bjobs -noheader -o stat "$OLD" 2>/dev/null) && [ "$state" != EXIT ] && [ "$state" != DONE ]; do sleep 10; done
  sleep 30  # container teardown
fi
echo "job $OLD $state; resuming from $OUT/last.pth  $(date)"
exec bash "$REPO/analysis/m6_nositef_graphs_20261118/train_predict.sh" "$@" --resume "$OUT/last.pth"

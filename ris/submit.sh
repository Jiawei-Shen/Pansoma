#!/bin/bash
# Submit a ris/ job script to LSF on RIS compute1, in the pytorch/pytorch:2.8.0-cuda12.6-cudnn9-devel container.
#   ris/submit.sh [-q QUEUE] [-G GROUP] [-g GPUS] [-n SLOTS] [-m MEM_GB] [-s GPU_SELECT] [-w DEPENDENCY] JOBNAME SCRIPT [args]
# Defaults: -q general, -g 1, -n 1, -m 64, -s gmem=40G (A40 / A100 / H100: bf16, no V100). -g 0 is a CPU job.
# -G: compute-epigenome-condo on the epigenome queues, compute-epigenome elsewhere. -w: bsub -w, e.g. "done(JOB)".
# -j: LSF job group (bsub -g); default $PANSOMA_JOB_GROUP, else /<user>/pansoma if it exists. RIS's /<user>/default runs
# at most 5 jobs at once; make a group of your own with `bgadd -L 10 /$USER/pansoma`.
# The label manifests keep the truth files' old-cluster paths; their copies on storage2 (MOUNTS below) are
# bind-mounted there, so the manifests (and their sha256) stay as they are.
set -euo pipefail
QUEUE=general GROUP="" GPUS=1 SLOTS=1 MEM=64 SELECT=gmem=40G WAIT="" JGROUP=${PANSOMA_JOB_GROUP:-}
while getopts "q:G:g:n:m:s:w:j:" o; do
  case $o in q) QUEUE=$OPTARG;; G) GROUP=$OPTARG;; g) GPUS=$OPTARG;; n) SLOTS=$OPTARG;; m) MEM=$OPTARG;;
             s) SELECT=$OPTARG;; w) WAIT=$OPTARG;; j) JGROUP=$OPTARG;; *) exit 2;; esac
done
if [ -z "$JGROUP" ] && bjgroup -s "/$USER/pansoma" 2>/dev/null | grep -q "^/$USER/pansoma "; then JGROUP=/$USER/pansoma; fi
shift $((OPTIND - 1))
NAME=$1 SCRIPT=$2; shift 2
[ -n "$GROUP" ] || case $QUEUE in epigenome*) GROUP=compute-epigenome-condo;; *) GROUP=compute-epigenome;; esac
HERE=$(cd "$(dirname "$0")" && pwd)
RUNS=${PANSOMA_RUNS:-/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs}
DATA=/storage2/fs1/epigenome/Active/jiawei.s/Dataset
MOUNTS=(  # storage2 copy : old-cluster path in the label manifests
  "$DATA/HG008_GIAB/draft_v02_benchmark:/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark"
  "$DATA/HG008_GIAB/dipcall_HG008N_GRCh38:/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38"
  "$DATA/pansoma_v2_tensors/COLO829T_truth:/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth"
  "$DATA/COLO829BL_DSA/dipcall_hg38:/scratch/qfu/COLO829BL_DSA/dipcall_hg38"
)
VOLS="/storage2/fs1/epigenome/Active:/storage2/fs1/epigenome/Active"
for m in "${MOUNTS[@]}"; do
  [ -d "${m%%:*}" ] && VOLS="$VOLS $m"
done
if [ "$GPUS" -gt 0 ]; then
  # gpuhost is an exclusive resource; -Q 99 requeues a GPU job whose container fails to start (exit 99, seen on exec-419)
  gpu=(-gpu "num=$GPUS:$SELECT:j_exclusive=yes" -Q 99); select="gpuhost "
else
  gpu=(); select=""
fi
mkdir -p "$RUNS/logs"
LSF_DOCKER_VOLUMES="$VOLS" LSF_DOCKER_SHM_SIZE=${LSF_DOCKER_SHM_SIZE:-16g} \
  bsub -G "$GROUP" -q "$QUEUE" -J "$NAME" ${JGROUP:+-g "$JGROUP"} "${gpu[@]}" ${WAIT:+-w "$WAIT"} -n "$SLOTS" \
    -R "${select}span[hosts=1] rusage[mem=${MEM}GB]" -M "${MEM}GB" \
    -a 'docker(pytorch/pytorch:2.8.0-cuda12.6-cudnn9-devel)' -oo "$RUNS/logs/$NAME.log" bash "$HERE/$SCRIPT" "$@"

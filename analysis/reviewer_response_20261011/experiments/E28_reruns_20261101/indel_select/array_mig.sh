#!/bin/bash
# One line of mig_tasks.tsv (RUN PROTOCOL PLATFORM SEED MICRO) per array task -> train_select.sh (24gb MIG originals).
set -euo pipefail
W=/scratch/jshen/Github/Pansoma/tmp/rerun_20261101/indel_select
read RUN PROTO p SEED MICRO < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $W/mig_tasks.tsv)
exec bash $W/train_select.sh $RUN $PROTO $p $SEED $MICRO

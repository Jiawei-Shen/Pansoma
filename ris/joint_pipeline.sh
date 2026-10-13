#!/bin/bash
# One joint model's pipeline as LSF jobs: joint_run.sh (train on all SETS with the recipe of KIND, predict chr1 of each
# set) on H100s of the epigenome queue, then one vcf.sh --set job per set once it is done (<run>/vcf_chr1_<SET>/).
# Runs on the login node: KIND=SNV|INDEL ris/joint_pipeline.sh NAME SET[,SET...] [extra train options]
#   KIND=INDEL ris/joint_pipeline.sh COLO829_Illumina_INDEL_joint \
#       COLO829T_Illumina,COLO829T_Illumina_P75,COLO829T_Illumina_P50,COLO829T_Illumina_P25
# GPUS (default 1; DDP with the recipe's --batch-size per GPU), SLOTS / MEM (default 7 / 192), QUEUE / SELECT, WAIT (an
# LSF dependency of the GPU job), VAL_WORKERS (default 6, per rank): as in snv_pipeline.sh. The recipes: recipes.sh.
set -euo pipefail
NAME=$1 SETS=$2; shift 2
HERE=$(cd "$(dirname "$0")" && pwd)
KIND=${KIND:-SNV}; source "$HERE/recipes.sh"
"$HERE/submit.sh" ${WAIT:+-w "$WAIT"} -q "${QUEUE:-epigenome}" -s "${SELECT:-gmodel=NVIDIAH10080GBHBM3}" -g "${GPUS:-1}" \
    -n "${SLOTS:-7}" -m "${MEM:-192}" \
    "pnv2_$NAME" joint_run.sh "$NAME" "$KIND" "$SETS" "${RECIPE[@]}" --read-threads 8 --val-workers "${VAL_WORKERS:-6}" "$@"
for s in ${SETS//,/ }; do
  "$HERE/submit.sh" -g 0 -n 4 -m 32 -w "done(pnv2_$NAME)" "pnv2_vcf_${NAME}_$s" vcf.sh --set "$s" "$NAME" "$KIND"
done

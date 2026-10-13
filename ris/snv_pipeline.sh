#!/bin/bash
# One set's SNV or INDEL pipeline (KIND, default SNV) as two LSF jobs: run.sh (train with the recipe of that kind,
# predict chr1) on one H100 of the epigenome queue (GPUS for more), then vcf.sh (graph VCF, GRCh38 VCF, PoN, rtg
# vcfeval) on CPU once it is done.
# Runs on the login node: ris/snv_pipeline.sh NAME SET [extra train options]
#   ris/snv_pipeline.sh HG008_Illumina_SNV_recipe HG008T_Illumina
# SLOTS / MEM (default 7 / 128) set the GPU job's slots and memory; QUEUE / SELECT its queue and GPU; WAIT an LSF
# dependency of the GPU job (e.g. WAIT='done(pnv2_relabel_HG008T_PacBio)'); KIND=INDEL runs the INDEL recipe.
# GPUS=N (default 1) trains on N GPUs under DDP (run.sh), with the recipe's --batch-size per GPU (global batch N x 1024);
# VAL_WORKERS (default 6) is per rank. NOVCF=1 submits the GPU job only (vcf.sh later, e.g. when a graph index arrives).
# The recipes: recipes.sh.
set -euo pipefail
NAME=$1 SET=$2; shift 2
HERE=$(cd "$(dirname "$0")" && pwd)
KIND=${KIND:-SNV}; source "$HERE/recipes.sh"
"$HERE/submit.sh" ${WAIT:+-w "$WAIT"} -q "${QUEUE:-epigenome}" -s "${SELECT:-gmodel=NVIDIAH10080GBHBM3}" -g "${GPUS:-1}" \
    -n "${SLOTS:-7}" -m "${MEM:-128}" \
    "pnv2_$NAME" run.sh "$NAME" "$SET" "$KIND" "${RECIPE[@]}" --read-threads 8 --val-workers "${VAL_WORKERS:-6}" "$@"
[ "${NOVCF:-0}" = 1 ] || "$HERE/submit.sh" -g 0 -n 4 -m 32 -w "done(pnv2_$NAME)" "pnv2_vcf_$NAME" vcf.sh "$NAME" "$KIND"

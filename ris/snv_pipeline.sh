#!/bin/bash
# One set's SNV pipeline as two LSF jobs: run.sh (train with the SNV recipe, predict chr1) on one H100 of the
# epigenome queue, then vcf.sh (graph VCF, GRCh38 VCF, PoN, rtg vcfeval) on CPU once it is done.
# Runs on the login node: ris/snv_pipeline.sh NAME SET [extra train options]
#   ris/snv_pipeline.sh HG008_Illumina_SNV_recipe HG008T_Illumina
# SLOTS / MEM (default 7 / 128) set the GPU job's slots and memory; QUEUE / SELECT its queue and GPU; WAIT an LSF
# dependency of the GPU job (e.g. WAIT='done(pnv2_relabel_HG008T_PacBio)'); KIND=INDEL runs the INDEL recipe.
# GPUS=N (default 1) trains on N GPUs under DDP (run.sh), with the recipe's --batch-size per GPU (global batch N x 1024);
# VAL_WORKERS (default 6) is per rank.
set -euo pipefail
NAME=$1 SET=$2; shift 2
HERE=$(cd "$(dirname "$0")" && pwd)
# the SNV recipe of machine_learning/pansoma_net_v2/README.md ("SNV recipe")
RECIPE=(--depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75
        --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth)
# the INDEL recipe ("INDEL recipe", 7fc7396: lr 1e-4, was 4e-4 in 45d56ad): the 8.7 M model, no candidate filter
[ "${KIND:-SNV}" = INDEL ] && RECIPE=(--depths 2 2 6 2 --dims 64 128 256 512 --epochs 12 --batch-size 1024 --lr 1e-4
        --class-weights sqrt --non-fraction 1.0 --block v2
        --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth)
GPUS=${GPUS:-1}
"$HERE/submit.sh" ${WAIT:+-w "$WAIT"} -q "${QUEUE:-epigenome}" -s "${SELECT:-gmodel=NVIDIAH10080GBHBM3}" -g "$GPUS" \
    -n "${SLOTS:-7}" -m "${MEM:-128}" \
    "pnv2_$NAME" run.sh "$NAME" "$SET" "${KIND:-SNV}" "${RECIPE[@]}" --read-threads 8 --val-workers "${VAL_WORKERS:-6}" "$@"
"$HERE/submit.sh" -g 0 -n 4 -m 32 -w "done(pnv2_$NAME)" "pnv2_vcf_$NAME" vcf.sh "$NAME" "${KIND:-SNV}"

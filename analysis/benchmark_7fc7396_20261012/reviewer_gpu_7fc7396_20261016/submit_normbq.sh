#!/bin/bash
# --norm dataset as of 2026-10-19 (user: MAPQ comes from the pipeline's own vg alignment, so only BQ is z-scored per
# dataset; the statistics come from tensors of all the dataset's chromosomes, in train and predict): the 8 models of
# group 3 retrained as <S>_<P>_SNV_d46_normbq / <S>_Illumina_INDEL_d46_normbq (SNV on 4 GPUs, INDEL on 1, as group 3),
# then predicted and evaluated by snv_pipeline.sh. Check of the diagnosis: the two ONT group-3 models
# (*_SNV_d46_normds, which predicted with chr1-only statistics) re-predicted with the new predict, which takes a
# training dataset's statistics from the checkpoint, into <run>/test_chr1_<SET> and vcf_chr1_<SET> (vcf.sh --set),
# the earlier test_chr1 / vcf_chr1 kept. Run from the repository root on the login node: bash <this file>
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs
SETS=(HG008:Illumina HG008:PacBio HG008:ONT COLO829:Illumina COLO829:fiberseq COLO829:ONT)
for SP in "${SETS[@]}"; do
  S=${SP%%:*} P=${SP#*:}
  GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 KIND=SNV ris/snv_pipeline.sh ${S}_${P}_SNV_d46_normbq ${S}T_${P}_d46 --norm dataset
done
for S in HG008 COLO829; do
  GPUS=1 SLOTS=5 MEM=192 VAL_WORKERS=6 KIND=INDEL ris/snv_pipeline.sh ${S}_Illumina_INDEL_d46_normbq ${S}T_Illumina_d46 \
      --norm dataset
done
for S in HG008 COLO829; do
  N=${S}_ONT_SNV_d46_normds SET=${S}T_ONT_d46
  ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n 7 -m 128 "pnv2_repredict_$N" predict.sh "$N" "$SET" SNV \
      --output "$RUNS/$N/test_chr1_$SET"
  ris/submit.sh -g 0 -n 4 -m 32 -w "done(pnv2_repredict_$N)" "pnv2_vcf_repredict_$N" vcf.sh --set "$SET" "$N" SNV
done

#!/bin/bash
# The RIS GPU jobs of analysis/reviewer_response_20261011/overview_20261015.md §3.1 (P0) that the data on RIS allow,
# submitted 2026-10-16 from the repository root on the login node: bash <this file>
# Every run: ris/snv_pipeline.sh (train with ris/recipes.sh, predict chr1, then vcf.sh: graph VCF, GRCh38 VCF, PoN, rtg
# vcfeval in and outside the BED), on the epigenome H100s; SNV on 4 GPUs (DDP, --batch-size 1024 per GPU), INDEL on 1.
# --norm train everywhere except group 4: every control (seed 0, with ch6, d46 / linear) was trained before --norm
# existed, i.e. with the training statistics, so the pairs differ only in the factor tested. Submitted in this order
# (LSF dispatches first come, first served):
#   0. repro: HG008_Illumina_INDEL_7fc7396 retrained as it was, to check that the current code with --norm train
#      reproduces a 7fc7396 checkpoint (best epoch, threshold, chr1 calls).
#   1. matched d9 (8): the SNV models of HG008T PacBio / ONT and of COLO829T on 4 GPUs, as the d46 and linear SNV runs
#      and HG008_Illumina_SNV_7fc7396_ddp4_b1024 (the d9 runs they replace were on 1 GPU: HG008_*_SNV_7fc7396,
#      COLO829_*_SNV_h90); the COLO829T INDEL models at the 7fc7396 recipe (lr 1e-4) on 1 GPU (only lr 4e-4 single-set
#      runs or purity-joint runs existed). These are also the COLO829T d9 runs the local analyses lack.
#   2. no ch6 (24): analysis/reviewer_response_20261011/experiments/E1.1_ch6_counterfactual/gpu_nopc_commands.sh
#      (--drop-planes path_count), d9 and d46, SNV and INDEL, with that file's names, GPU counts, slots and memory, plus
#      --norm train; its 3 COLO829T d9 INDEL with-ch6 controls are group 1 (same names). The d9 SNV no-ch6 runs other
#      than HG008T Illumina are on 1 GPU, paired with HG008_*_SNV_7fc7396 / COLO829_*_SNV_h90 (first submitted on 4
#      GPUs as *_SNV_nopc_7fc7396_ddp4_b1024, jobs 133460-133477, killed while pending and resubmitted as in the file).
#   3. --norm dataset (8): the d46 models of E12a re-trained with BQ / MAPQ z-scored per dataset, for the cross-sample
#      tests (predicted on tequila): SNV x 6 sets, Illumina INDEL x 2.
#   4. seeds (36): the main SNV model of each reference arm (d9 = group 1 / HG008_Illumina_SNV_7fc7396_ddp4_b1024,
#      d46, linear) x 6 sets, seeds 1 and 2 (seed 0 = the existing run).
# Not here (no tensor sets on RIS): the HCC1395T / CASTLE in-sample models and leave-one-line-out.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
SNV4=(env GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 KIND=SNV ris/snv_pipeline.sh)
INDEL1=(env GPUS=1 SLOTS=5 MEM=192 VAL_WORKERS=6 KIND=INDEL ris/snv_pipeline.sh)
TRAIN=(--norm train)
NOPC=(--drop-planes path_count)
SETS=(HG008:Illumina HG008:PacBio HG008:ONT COLO829:Illumina COLO829:fiberseq COLO829:ONT)

# 0. repro
"${INDEL1[@]}" HG008_Illumina_INDEL_7fc7396_repro HG008T_Illumina "${TRAIN[@]}"

# 1. matched d9
for SP in HG008:PacBio HG008:ONT COLO829:Illumina COLO829:fiberseq COLO829:ONT; do
  S=${SP%%:*} P=${SP#*:}
  "${SNV4[@]}" ${S}_${P}_SNV_7fc7396_ddp4_b1024 ${S}T_${P} "${TRAIN[@]}"
done
for P in Illumina fiberseq ONT; do
  "${INDEL1[@]}" COLO829_${P}_INDEL_7fc7396 COLO829T_${P} "${TRAIN[@]}"
done

# 2. no ch6
"${SNV4[@]}" HG008_Illumina_SNV_nopc_7fc7396_ddp4_b1024 HG008T_Illumina "${TRAIN[@]}" "${NOPC[@]}"
for SP in HG008:PacBio HG008:ONT COLO829:Illumina COLO829:fiberseq COLO829:ONT; do
  S=${SP%%:*} P=${SP#*:}
  MEM=192 ris/snv_pipeline.sh ${S}_${P}_SNV_nopc_7fc7396 ${S}T_${P} "${TRAIN[@]}" "${NOPC[@]}"
done
for SP in "${SETS[@]}"; do
  S=${SP%%:*} P=${SP#*:}
  "${INDEL1[@]}" ${S}_${P}_INDEL_nopc_7fc7396 ${S}T_${P} "${TRAIN[@]}" "${NOPC[@]}"
done
for SP in "${SETS[@]}"; do
  S=${SP%%:*} P=${SP#*:}
  "${SNV4[@]}" ${S}_${P}_SNV_nopc_d46_7fc7396 ${S}T_${P}_d46 "${TRAIN[@]}" "${NOPC[@]}"
  "${INDEL1[@]}" ${S}_${P}_INDEL_nopc_d46_7fc7396 ${S}T_${P}_d46 "${TRAIN[@]}" "${NOPC[@]}"
done

# 3. --norm dataset (the default of train)
for SP in "${SETS[@]}"; do
  S=${SP%%:*} P=${SP#*:}
  "${SNV4[@]}" ${S}_${P}_SNV_d46_normds ${S}T_${P}_d46 --norm dataset
done
for S in HG008 COLO829; do
  "${INDEL1[@]}" ${S}_Illumina_INDEL_d46_normds ${S}T_Illumina_d46 --norm dataset
done

# 4. seeds
for seed in 1 2; do
  for SP in "${SETS[@]}"; do
    S=${SP%%:*} P=${SP#*:}
    "${SNV4[@]}" ${S}_${P}_SNV_7fc7396_ddp4_b1024_s$seed ${S}T_${P} "${TRAIN[@]}" --seed $seed
    "${SNV4[@]}" ${S}_${P}_SNV_d46_7fc7396_s$seed ${S}T_${P}_d46 "${TRAIN[@]}" --seed $seed
    "${SNV4[@]}" ${S}_${P}_SNV_linear_7fc7396_s$seed ${S}T_${P}_linear "${TRAIN[@]}" --seed $seed
  done
done

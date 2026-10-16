#!/bin/bash
# No-ch6 arm (E3 / E1.1 follow-up): the 7fc7396 recipes retrained with the path-count plane left out
# (train --drop-planes path_count), 12 sample x platform x kind sets on HPRC v1.1 d9 and on HPRC v2.1 d46.
# NOT run here (GPU, RIS compute1). Run on the RIS login node from the repository:
#   cd /storage2/fs1/epigenome/Active/jiawei.s/GitHub/Pansoma && bash <this file>
# ris/snv_pipeline.sh appends extra train options after the recipe (ris/recipes.sh) and submits vcf.sh (graph VCF,
# GRCh38 VCF, PoN, rtg vcfeval in and outside the BED) after each GPU job. The index / chunk caches of every set
# already exist on RIS (the with-ch6 runs built them; --drop-planes does not change them).
# GPU counts follow the with-ch6 run each one is compared with, so the pair differs only in ch6:
#   d9 HG008T SNV: Illumina on 4 GPUs (HG008_Illumina_SNV_7fc7396_ddp4_b1024), PacBio / ONT on 1 (HG008_*_SNV_7fc7396)
#   d9 COLO829T SNV: 1 GPU, as COLO829_<platform>_SNV_h90 (2026-10-03, the same SNV recipe)
#   d9 INDEL: 1 GPU (HG008_*_INDEL_7fc7396). COLO829T d9 has no single-set INDEL run at the 7fc7396 recipe (the
#     2026-10-06 runs used lr 4e-4, the 7fc7396 ones are purity-joint), so its with-ch6 controls are included (3 runs)
#   d46: SNV on 4 GPUs, INDEL on 1, as every *_d46_7fc7396 run
# 27 GPU runs: 24 no-ch6 + 3 COLO829T d9 INDEL controls.
set -euo pipefail
NOPC=(--drop-planes path_count)

# ---- HPRC v1.1 d9 ----
GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 ris/snv_pipeline.sh HG008_Illumina_SNV_nopc_7fc7396_ddp4_b1024 HG008T_Illumina "${NOPC[@]}"
MEM=192 ris/snv_pipeline.sh HG008_PacBio_SNV_nopc_7fc7396 HG008T_PacBio "${NOPC[@]}"
MEM=192 ris/snv_pipeline.sh HG008_ONT_SNV_nopc_7fc7396 HG008T_ONT "${NOPC[@]}"
MEM=192 ris/snv_pipeline.sh COLO829_Illumina_SNV_nopc_7fc7396 COLO829T_Illumina "${NOPC[@]}"
MEM=192 ris/snv_pipeline.sh COLO829_fiberseq_SNV_nopc_7fc7396 COLO829T_fiberseq "${NOPC[@]}"
MEM=192 ris/snv_pipeline.sh COLO829_ONT_SNV_nopc_7fc7396 COLO829T_ONT "${NOPC[@]}"
for P in Illumina PacBio ONT; do
  KIND=INDEL SLOTS=5 MEM=192 ris/snv_pipeline.sh HG008_${P}_INDEL_nopc_7fc7396 HG008T_${P} "${NOPC[@]}"
done
for P in Illumina fiberseq ONT; do
  KIND=INDEL SLOTS=5 MEM=192 ris/snv_pipeline.sh COLO829_${P}_INDEL_nopc_7fc7396 COLO829T_${P} "${NOPC[@]}"
  KIND=INDEL SLOTS=5 MEM=192 ris/snv_pipeline.sh COLO829_${P}_INDEL_7fc7396 COLO829T_${P}     # with-ch6 control
done

# ---- HPRC v2.1 d46 ----
for SP in HG008:Illumina HG008:PacBio HG008:ONT COLO829:Illumina COLO829:fiberseq COLO829:ONT; do
  S=${SP%%:*} P=${SP#*:}
  GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 ris/snv_pipeline.sh ${S}_${P}_SNV_nopc_d46_7fc7396 ${S}T_${P}_d46 "${NOPC[@]}"
  KIND=INDEL GPUS=1 SLOTS=5 MEM=192 VAL_WORKERS=6 ris/snv_pipeline.sh ${S}_${P}_INDEL_nopc_d46_7fc7396 ${S}T_${P}_d46 "${NOPC[@]}"
done

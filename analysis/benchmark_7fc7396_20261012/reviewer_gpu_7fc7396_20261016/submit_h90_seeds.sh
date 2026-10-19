#!/bin/bash
# COLO829T Illumina SNV on v1.1 d9 with the 1-GPU recipe of COLO829_Illumina_SNV_h90 (global batch 1,024; the matched
# d9 runs use 4 GPUs, 4,096; lr 1e-4 in both), seeds 1 and 2, beside seed 0 = COLO829_Illumina_SNV_h90 (PASS F1 0.778
# against 0.566-0.626 for the three 4-GPU seeds). If the two new seeds also land near 0.78 the 4-GPU recipe is what
# is unstable on this set; if they are low, 0.778 was a lucky seed. --norm train as the seed-0 run (trained before
# --norm; note: seed 0 was trained by the code of 9d6bd8d, these by the current code). Submitted 2026-10-19 from the
# repository root on the login node: bash <this file>
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
for seed in 1 2; do
  GPUS=1 MEM=192 KIND=SNV ris/snv_pipeline.sh COLO829_Illumina_SNV_h90_s$seed COLO829T_Illumina --norm train --seed $seed
done

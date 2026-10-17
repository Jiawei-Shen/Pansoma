#!/bin/bash
# E13 (R1.7): the learning curve and the model size of analysis/reviewer_response_20261011/experiments/
# E21_training_composition/README.md section 6, submitted 2026-10-16 after submit.sh. Its commands, with two changes:
#   - the exclude files are made here on RIS (e13_subsets.sh, CPU job; they are SHA-1 draws over the seed-0 split, so
#     e13_subsets_HG008T_PacBio_d46_SNV.tsv here must equal E21's) instead of copied from tequila; the GPU jobs wait for it;
#   - --norm train, as the 100 % / 50.4 M point they extend (HG008_PacBio_SNV_d46_7fc7396, trained before --norm).
# 8 runs: 10 / 25 / 50 / 75 % of the training blocks x 2 replicates; 2 runs: 8.7 M and 199.6 M on all blocks.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
E13=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs/e13_subsets
SNV4=(env GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 KIND=SNV)
ris/submit.sh -g 0 -n 2 -m 32 pnv2_e13_subsets ../analysis/benchmark_7fc7396_20261012/reviewer_gpu_7fc7396_20261016/e13_subsets.sh
for rep in 1 2; do for f in 10 25 50 75; do
  "${SNV4[@]}" WAIT="done(pnv2_e13_subsets)" ris/snv_pipeline.sh HG008_PacBio_SNV_d46_e13_keep${f}_rep${rep} HG008T_PacBio_d46 \
      --norm train --exclude-candidates $E13/HG008T_PacBio_d46.SNV.keep${f}.rep${rep}.exclude.txt
done; done
"${SNV4[@]}" ris/snv_pipeline.sh HG008_PacBio_SNV_d46_e13_8p7M HG008T_PacBio_d46 --norm train \
    --depths 2 2 6 2 --dims 64 128 256 512
"${SNV4[@]}" ris/snv_pipeline.sh HG008_PacBio_SNV_d46_e13_199p6M HG008T_PacBio_d46 --norm train \
    --dims 192 384 768 1536 --micro-batches 2

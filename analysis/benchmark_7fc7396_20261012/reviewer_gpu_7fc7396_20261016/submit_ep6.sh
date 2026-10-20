#!/bin/bash
# Epoch 6 of the COLO829T Illumina d9 --norm dataset SNV seeds (COLO829_Illumina_SNV_d9_normbq, _s1, _s2: PASS F1
# 0.566 / 0.628 / 0.568, all three keep epoch 8), the user, 2026-10-20: "test COLO829T Illumina d9 with epoch 6, first
# in the queue". Each seed is retrained exactly as submit_normbq_all.sh trained it (snv_pipeline.sh: SNV recipe,
# 4 GPUs DDP, --batch-size 1024 per GPU, --read-threads 8 --val-workers 3 --norm dataset [--seed K]) with every epoch's
# checkpoint kept and the run ended once epoch 6 is saved (ris/epoch_sweep.sh, ris/tools/save_epochs.py
# --stop-after-epoch 6; --epochs 20 and so the learning-rate schedule unchanged); metrics.jsonl epochs 1-6 are then
# compared with the original run's. Epoch 6 is predicted (chr1) and scored (vcf.sh) as <run>_sweep_e06, its PASS
# threshold that epoch's validation recall-0.9 one, as the epoch sweep of submit_epoch_sweep.sh.
# Submitted 2026-10-20 from the repository root on the login node: bash <this file>; then btop (first in the queue).
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs
SET=COLO829T_Illumina
KIND=SNV; source ris/recipes.sh
for seed in 0 1 2; do
  RUN=COLO829_Illumina_SNV_d9_normbq$([ $seed = 0 ] || echo _s$seed)
  NAME=${RUN}_sweep E=${RUN}_sweep_e06
  ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 4 -n 16 -m 192 "pnv2_$NAME" epoch_sweep.sh "$NAME" "$SET" SNV \
      "${RECIPE[@]}" --read-threads 8 --val-workers 3 --norm dataset $([ $seed = 0 ] || echo --seed $seed) \
      --stop-after-epoch 6
  mkdir -p "$RUNS/$E"
  ln -sfn "../$NAME/args.json" "$RUNS/$E/args.json"
  ln -sfn "../$NAME/epochs/epoch_06.pth" "$RUNS/$E/best.pth"
  ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n 7 -m 128 -w "done(pnv2_$NAME)" "pnv2_predict_$E" \
      predict.sh "$E" "$SET" SNV
  ris/submit.sh -g 0 -n 4 -m 32 -w "done(pnv2_predict_$E)" "pnv2_vcf_$E" vcf.sh "$E" SNV
done
# Then (2026-10-20 15:05, the user: "seed 1 at epoch 5?") seed 1's epoch 5, the same way:
#   E=COLO829_Illumina_SNV_d9_normbq_s1_sweep_e05: args.json -> ../COLO829_Illumina_SNV_d9_normbq_s1_sweep/args.json,
#   best.pth -> ../COLO829_Illumina_SNV_d9_normbq_s1_sweep/epochs/epoch_05.pth; predict.sh "$E" COLO829T_Illumina SNV
#   (job 204610), vcf.sh "$E" SNV (204611). The user's choice: seeds 0 / 1 / 2 at epochs 6 / 5 / 6 (collect.py CHOSEN).

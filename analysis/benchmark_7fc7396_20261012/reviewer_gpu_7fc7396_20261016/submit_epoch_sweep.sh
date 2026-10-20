#!/bin/bash
# Epoch sweep of COLO829T Illumina SNV on v1.1 d9 with the h90 recipe (1 GPU, seed 0, --norm train): is the PASS F1
# 0.778 (h90, best epoch 6) against 0.624 (the rerun h90_repro, best epoch 8) the choice between two near-tied
# epochs, or something else? The run trains as snv_pipeline.sh would (ris/recipes.sh SNV recipe, --read-threads 8
# --val-workers 6) but keeps every epoch's checkpoint (ris/epoch_sweep.sh, ris/tools/save_epochs.py); each epoch is
# then predicted (chr1) and scored (vcf.sh) as its own run <name>_eNN. Submitted 2026-10-19 from the repository root
# on the login node: bash <this file>
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs
NAME=COLO829_Illumina_SNV_h90_sweep SET=COLO829T_Illumina
KIND=SNV; source ris/recipes.sh
ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n 7 -m 192 "pnv2_$NAME" epoch_sweep.sh "$NAME" "$SET" SNV \
    "${RECIPE[@]}" --read-threads 8 --val-workers 6 --norm train --seed 0
for e in $(seq -w 1 20); do
  E=${NAME}_e$e
  mkdir -p "$RUNS/$E"
  ln -sfn "../$NAME/args.json" "$RUNS/$E/args.json"
  ln -sfn "../$NAME/epochs/epoch_$e.pth" "$RUNS/$E/best.pth"
  ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n 7 -m 128 -w "done(pnv2_$NAME)" "pnv2_predict_$E" \
      predict.sh "$E" "$SET" SNV
  ris/submit.sh -g 0 -n 4 -m 32 -w "done(pnv2_predict_$E)" "pnv2_vcf_$E" vcf.sh "$E" SNV
done

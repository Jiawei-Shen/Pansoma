#!/bin/bash
# The COLO829T purity dilutions (P100 + P75 + P50 + P25) rerun with the current SNV pipeline (ris/recipes.sh; f0171af:
# one SNV recipe for every platform), submitted 2026-11-12. The user: the original dilution sets (HPRC v1.1 d9, SNV build AF 0.08:
# COLO829T_<P>, _P75, _P50, _P25), not AF 0.05 rebuilds; the pipeline is the latest one. The design of the earlier
# dilution runs (COLO829_<P>_SNV_joint_7fc7396_gpu4) is kept: one model per platform on the four purities, chr1 of
# every set left out, 5 % of the chr2-22 node blocks validate, chr1 of each purity tested with one threshold per model;
# 4 H100s per model under DDP (--batch-size 1024 per GPU), the slots and memory of those runs.
#   S. SNV site features of the 8 Illumina / fiberseq sets (CPU, general) into $PANSOMA_RUNS/site_features     8 jobs
#   A. SNV joint models, COLO829_<P>_SNV_joint_newpipe:                                                        2 runs
#        Illumina, fiberseq: site features, non tensors with AF <= 0.1 sampled 25 % per epoch at loss weight 1,
#          best epoch by the post-PoN validation F1, PASS at the post-PoN validation recall 0.95 of the callable
#          positives (on AF 0.08 sets here, so only the recipe differs from the earlier runs, not the tensors);
#        ONT, first submission (repository 0c32edf, recipes.sh then by platform): the original SNV recipe (no site
#          features, every negative, PASS at validation recall 0.9). Killed in epoch 2 after f0171af made the ONT SNV
#          recipe the Illumina / fiberseq one; step A now submits Illumina and fiberseq only.
#      Each: joint_run.sh (train, predict chr1 of the four sets), then vcf.sh --set per purity.
#   B. ONT with the SNV recipe of f0171af (user 2026-11-12: ONT too; one recipe for every platform): SNV site features
#      of the three ONT dilution sets (COLO829T_ONT.SNV.npz exists), then COLO829_ONT_SNV_joint_newpipe, same GPUs /
#      slots / memory. The ONT sets are built at SNV AF 0.08 (manifest min_af 0.08: no tensor below it), so no
#      negative with AF < 0.08 enters training (the user); the non tensors at AF 0.08-0.1 are sampled 25 % per epoch.
#   C. the best of the first 5 epochs, tested while the runs train on (user 2026-11-12): snapshot_best5.sh keeps each
#      run's best.pth once epoch 5 is done (COLO829_<P>_SNV_joint_newpipe_e5/; fiberseq's was copied by hand), then
#      predict_set.sh on chr1 of the four purities (one H100 each, epigenome) and vcf.sh --set.  3 x (4 + 4) jobs
#   10 epochs (user 2026-11-12: enough for these dilution models; the recipes keep --epochs 20): the Illumina and
#      fiberseq trainings of step A were killed after epochs 10 and 15; their best of epochs 1-10 is the step-C
#      snapshot (epochs 4 and 3, best.pth identical), so the step-C tests are their results. The 20-epoch ONT job of
#      step B was killed while pending; step B now trains ONT through joint_run_stop.sh 10 (the 20-epoch
#      learning-rate schedule, ended after epoch 10: the same first 10 epochs as the other two platforms).
#   D. final PASS thresholds (user 2026-11-12, after the tests at the call rule's threshold - 0.10 / - 0.15, vcf.sh
#      --set SET RUN SNV THRESHOLD _tm010 / _tm015): Illumina and fiberseq call at the validation threshold - 0.10
#      (0.816548 / 0.503872), ONT keeps the recipe's threshold. set_call_threshold.sh writes the threshold into
#      COLO829_<P>_SNV_joint_newpipe_e5/best.pth (the unchanged checkpoint stays as COLO829_<P>_SNV_joint_newpipe/
#      best.pth); the calls at the recipe's threshold move to vcf_chr1_<SET>_cr095 and vcf.sh --set without a
#      threshold fills vcf_chr1_<SET> again.                                                        2 + 8 jobs
# INDEL is not rerun: COLO829_<P>_INDEL_joint_ppr8 (analysis/indel_ppr8_9919d4a_20261026, step C) already used the
# current INDEL recipe on these sets.
# Run from the repository root on the login node: bash analysis/colo829t_dilution_newpipe_20261112/submit.sh [S|A|B|C|D ...]
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
STEPS=${*:-S A B C D}
has() { [[ " $STEPS " == *" $1 "* ]]; }
sets_of() { echo "COLO829T_$1,COLO829T_$1_P75,COLO829T_$1_P50,COLO829T_$1_P25"; }
sitef() { echo "pnv2_sitefsnv_$1"; }

# S. SNV site features (the INDEL files of these sets exist; the SNV rule also reads dbSNP and 1000G)
if has S; then
  for P in Illumina fiberseq; do for s in $(sets_of $P | tr ',' ' '); do
    KIND=SNV ris/submit.sh -g 0 -n 8 -m 32 "$(sitef $s)" site_features.sh "$s"
  done; done
fi
# A. the joint SNV models of Illumina and fiberseq
if has A; then
  for P in Illumina fiberseq; do
    case $P in Illumina) MEM=320;; fiberseq) MEM=256;; esac
    sets=$(sets_of $P)
    w=$(for s in ${sets//,/ }; do printf 'done(%s) && ' "$(sitef $s)"; done); w=${w% && }
    env GPUS=4 SLOTS=16 MEM=$MEM VAL_WORKERS=3 WAIT="$w" ris/joint_pipeline.sh COLO829_${P}_SNV_joint_newpipe "$sets"
  done
fi
# B. ONT: site features of the three dilution sets, then the joint SNV model
if has B; then
  sets=$(sets_of ONT) w=""
  for s in COLO829T_ONT_P75 COLO829T_ONT_P50 COLO829T_ONT_P25; do
    [ -s "${PANSOMA_RUNS:-/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs}/site_features/$s.SNV.npz" ] && continue
    KIND=SNV ris/submit.sh -g 0 -n 8 -m 48 "$(sitef $s)" site_features.sh "$s"
    w="${w}done($(sitef $s)) && "
  done
  name=COLO829_ONT_SNV_joint_newpipe
  KIND=SNV; source ris/recipes.sh
  ris/submit.sh ${w:+-w "${w% && }"} -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 4 -n 16 -m 384 pnv2_$name \
      ../analysis/colo829t_dilution_newpipe_20261112/joint_run_stop.sh "${STOP:-10}" $name SNV "$sets" "${RECIPE[@]}" \
      --read-threads 8 --val-workers 3
  for s in ${sets//,/ }; do
    ris/submit.sh -g 0 -n 4 -m 32 -w "done(pnv2_$name)" pnv2_vcf_${name}_$s vcf.sh --set $s $name SNV
  done
fi
# C. best of epochs 1-5: snapshot, predict chr1 of each purity, calls + vcfeval
if has C; then
  for P in ${E5_PLATFORMS:-fiberseq Illumina ONT}; do
    run=COLO829_${P}_SNV_joint_newpipe; w=""
    if [ ! -s "${PANSOMA_RUNS:-/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs}/${run}_e5/best.pth" ]; then
      ris/submit.sh -g 0 -n 1 -m 4 pnv2_snap5_$run ../analysis/colo829t_dilution_newpipe_20261112/snapshot_best5.sh $run
      w="done(pnv2_snap5_$run)"
    fi
    for s in $(sets_of $P | tr ',' ' '); do
      ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -n 4 -m 64 ${w:+-w "$w"} pnv2_pred_${run}_e5_$s \
          predict_set.sh ${run}_e5 $s SNV
      ris/submit.sh -g 0 -n 4 -m 32 -w "done(pnv2_pred_${run}_e5_$s)" pnv2_vcf_${run}_e5_$s vcf.sh --set $s ${run}_e5 SNV
    done
  done
fi
# D. Illumina and fiberseq: the checkpoint's PASS threshold - 0.10, then the calls of each purity at it
if has D; then
  R=${PANSOMA_RUNS:-/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs}
  for P in Illumina fiberseq; do
    run=COLO829_${P}_SNV_joint_newpipe_e5
    ris/submit.sh -g 0 -n 2 -m 16 pnv2_setthr_$run ../analysis/colo829t_dilution_newpipe_20261112/set_call_threshold.sh \
        $run "${OFFSET:--0.10}" COLO829_${P}_SNV_joint_newpipe
    for s in $(sets_of $P | tr ',' ' '); do
      if [ -d "$R/$run/vcf_chr1_$s" ] && [ ! -e "$R/$run/vcf_chr1_${s}_cr095" ]; then
        mv "$R/$run/vcf_chr1_$s" "$R/$run/vcf_chr1_${s}_cr095"
      fi
      ris/submit.sh -g 0 -n 4 -m 32 -w "done(pnv2_setthr_$run)" pnv2_vcf_${run}_$s vcf.sh --set $s $run SNV
    done
  done
fi
bjobs -w

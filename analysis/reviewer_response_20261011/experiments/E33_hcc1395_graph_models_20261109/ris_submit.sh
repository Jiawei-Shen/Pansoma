#!/bin/bash
# Step 2 of 3 (RIS login node, after to_ris.sh and git pull): submit the five HCC1395T in-sample SNV models on the
# epigenome H100s. Each: run.sh (train on chr2-22 without chr20, chr1 + chr20 held out of training and validation;
# then predict chr1) and, after it, predict.sh for chr20. The same train options as the tequila runs
# (tmp/hcc1395_linear_d9_20261104/train_af05/train.sh = lnf25w1, train_normfix.sh = ONT normfix), seed 0; only
# --read-threads / --val-workers differ (speed). The code must be 1497adb's pansoma_net_v2 (= tequila's code_ml3).
# usage: bash analysis/reviewer_response_20261011/experiments/E33_hcc1395_graph_models_20261109/ris_submit.sh
set -euo pipefail
cd "$(dirname "$0")/../../../.."
c=$(git log -1 --format=%H -- 'machine_learning/pansoma_net_v2/*.py')
[[ $c == 1497adb* ]] || { echo "machine_learning/pansoma_net_v2 code last changed in ${c:0:7}, not 1497adb (the tequila code)"; exit 1; }
grep -q CASTLE/HCC1395 ris/submit.sh || { echo "ris/submit.sh has no HCC1395T truth mount: git pull"; exit 1; }
A=/storage2/fs1/epigenome/Active/jiawei.s
T=$A/Dataset/pansoma_v2_tensors RUNS=$A/pansoma_net_v2_runs
[ -s $A/Dataset/CASTLE/HCC1395/truth/SEQC2_High-Confidence_Regions_v1.2.bed ] || { echo "no HCC1395T truth on storage2: run to_ris.sh"; exit 1; }
for s in HCC1395T_ONT_linear HCC1395T_ONT_d9 HCC1395T_ONT_d46 HCC1395T_Illumina_d9_snvaf05 HCC1395T_Illumina_d46_snvaf05; do
  [ -s $T/$s/tensors/SNV/labels.manifest.json ] || { echo "no $T/$s/tensors/SNV: run to_ris.sh"; exit 1; }
done
for s in HCC1395T_Illumina_d9_snvaf05 HCC1395T_Illumina_d46_snvaf05; do
  [ -s $RUNS/site_features/$s.SNV.npz ] || { echo "no site features $RUNS/site_features/$s.SNV.npz: run to_ris.sh"; exit 1; }
done

BASE=(--depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 --class-weights pow:0.75
      --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth
      --seed 0 --read-threads 8 --val-workers 6)
# PacBio HiFi / Illumina (lnf25w1): AF 0.05 sets, SNV site features, non with AF <= 0.1 sampled 25% at weight 1,
# best epoch by post-PoN validation F1, PASS at validation recall 0.9 (rule 1; rule 2 is applied in compare_graphs.py)
LNF=(--site-features $RUNS/site_features --best-by postpon --call-recall 0.9
     --non-fraction 0.25 --keep-non-af 0.100001 --no-importance-weights)
sub() {  # NAME SET [options]
  local N=$1 S=$2; shift 2
  ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -n 7 -m 128 pnv2_$N run.sh $N $S SNV \
      --test-chroms --set-test-chroms $T/$S/tensors chr1 chr20 "${BASE[@]}" "$@"
  ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -n 4 -m 64 -w "done(pnv2_$N)" pnv2_pred20_$N \
      predict.sh $N $S SNV --chroms chr20 --output $RUNS/$N/test_chr20
}
# ONT: the normfix recipe (AF 0.08 sets, no site features, every negative, best epoch by validation F1, PASS at recall 0.9)
sub HCC1395T_ONT_SNV_linear_normfix    HCC1395T_ONT_linear
sub HCC1395T_Illumina_SNV_d9_lnf25w1   HCC1395T_Illumina_d9_snvaf05  "${LNF[@]}"
sub HCC1395T_Illumina_SNV_d46_lnf25w1  HCC1395T_Illumina_d46_snvaf05 "${LNF[@]}"
sub HCC1395T_ONT_SNV_d46_normfix       HCC1395T_ONT_d46
sub HCC1395T_ONT_SNV_d9_normfix        HCC1395T_ONT_d9
bjobs -w
echo "logs: $RUNS/logs/pnv2_<name>.log, training: $RUNS/<name>/train.log; when every pnv2_pred20_* is done: from_ris.sh on wangcluster"

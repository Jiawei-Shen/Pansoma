#!/bin/bash
# Seeds 1 and 2 of the six-tumor SNV models of the current pipeline on HPRC v2.1 (d46), Illumina / PacBio HiFi
# (COLO829T: Fiber-seq) / ONT: 6 models (user 2026-11-23: multi-seed six-tumor v2.1 models of the latest pipeline,
# tested on the three test tumors). Seed 0 is snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe (old cluster, one
# H100), the six-tumor model of Supplementary Fig. R1.6 / Table R1.7b; these two seeds repeat it with only --seed
# changed: ris/recipes.sh SNV recipe with the site features as a model input (NOSITEF=0), the six d46 sets (Illumina /
# PacBio HiFi on the AF 0.05 builds of the four CASTLE lines and COLO829T, HG008T and ONT at 0.08), HCC1395T chr1 +
# chr20, HG008T chr1 and COLO829T chr1 left out of training and validation, one H100, batch 1024 (not DDP, as seed 0).
# Each model: ../m6_nositef_graphs_20261118/train_predict.sh (train, then predict the held chromosomes of the three
# test tumors), then ris/vcf.sh --set per test tumor (chr1, PASS = rule 2: post-PoN validation recall 0.95 of the
# callable positives; population filter; rtg vcfeval in the benchmark BED).
# Runs: $PANSOMA_RUNS/m6_seeds_20261123/multicancer6_<P>_SNV_newpipe_s<SEED>. A run whose directory exists, or whose
# job is in LSF, is skipped; a set that is not complete on storage2 (ready.py) stops the script.
# Run from the repository root on the login node:
#   bash analysis/m6_seeds_20261123/submit.sh            [SEEDS="1 2"] [PLATFORMS="Illumina PacBio ONT"] [SLOTS=5]
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
A=/storage2/fs1/epigenome/Active/jiawei.s
D=$A/Dataset/pansoma_v2_tensors RUNS=${PANSOMA_RUNS:-$A/pansoma_net_v2_runs}
B=m6_seeds_20261123
TP=../analysis/m6_nositef_graphs_20261118/train_predict.sh   # relative to ris/ (submit.sh runs ris/<script>)

sets() {  # PLATFORM -> HCC1395T HCC1954T H1437T H2009T HG008T COLO829T, the sets of multicancer6_<P>_SNV_newpipe
  local P=$1 cp=$1 af=_snvaf05 l
  [ "$P" = PacBio ] && cp=fiberseq
  [ "$P" = ONT ] && af=""
  for l in HCC1395T HCC1954T H1437T H2009T; do printf '%s ' "${l}_${P}_d46$af"; done
  echo "HG008T_${P}_d46 COLO829T_${cp}_d46$af"
}

for P in ${PLATFORMS:-Illumina PacBio ONT}; do
  read -r -a S <<< "$(sets "$P")"
  for s in "${S[@]}"; do
    r=$(python3 analysis/m6_nositef_graphs_20261118/ready.py "$D/$s/tensors" SNV) || { echo "$s/tensors/SNV: $r"; exit 1; }
    [ -s "$RUNS/site_features/$s.SNV.npz" ] && [ ! -L "$RUNS/site_features/$s.SNV.npz" ] || { echo "no site features of $s"; exit 1; }
  done
  for seed in ${SEEDS:-1 2}; do
    name=multicancer6_${P}_SNV_newpipe_s$seed job=pnv2_m6_${P}_SNV_newpipe_s$seed
    if [ -e "$RUNS/$B/$name" ]; then echo "$name: $RUNS/$B/$name exists, skipped"; continue; fi
    if bjobs -J "$job" 2>/dev/null | grep -q "$job"; then echo "$name: job $job is in LSF, skipped"; continue; fi
    case $P in Illumina) mem=192;; *) mem=128;; esac  # LSF memory counts page cache; need ~40-75 GB (2026-11-23)
    KIND=SNV NOSITEF=0 PLATFORM=$P; source ris/recipes.sh
    held="${S[0]}:chr1+chr20,${S[4]}:chr1,${S[5]}:chr1"
    ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n "${SLOTS:-5}" -m $mem "$job" \
        "$TP" "$B/$name" SNV "$(IFS=,; echo "${S[*]}")" "$held" "${RECIPE[@]}" --seed "$seed" \
        --read-threads 8 --val-workers 6
    for s in "${S[0]}" "${S[4]}" "${S[5]}"; do
      ris/submit.sh -g 0 -n 4 -m 32 -w "done($job)" "pnv2_vcf_m6_${P}_SNV_newpipe_s${seed}_$s" vcf.sh --set "$s" "$B/$name" SNV
    done
  done
done
echo "logs: $RUNS/logs/pnv2_m6_<P>_SNV_newpipe_s<SEED>.log, training: $RUNS/$B/<name>/train.log"

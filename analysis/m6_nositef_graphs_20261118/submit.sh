#!/bin/bash
# The six-tumor models of the current pipeline WITHOUT the population site features as a model input, on HPRC v1.1 (d9)
# and HPRC v2.1 (d46), SNV and INDEL, Illumina / PacBio HiFi (COLO829T: Fiber-seq) / ONT: 12 models (user 2026-11-18).
# They pair with the GRCh38-only six-tumor models, which never had the features as an input
# (m4m6_graphs_20261110/multicancer6_<P>_<KIND>_linear_nositef on the old cluster), so the three references differ in
# the graph only. Each model: ris/recipes.sh with NOSITEF=1 (the recipe + --site-features-pon-only: the feature files
# give only the PoN tags for the post-PoN epoch selection and the call rule), the design of the other six-tumor models
# (HCC1395T chr1 + chr20, HG008T chr1 and COLO829T chr1 left out of training and validation; HCC1954T, H1437T and
# H2009T train on every chromosome; seed 0; one H100, batch 1024), then chr1 of the three test tumors predicted
# (HCC1395T also chr20) and called: vcf.sh --set (graph VCF, GRCh38 VCF, population filter, rtg vcfeval in the
# benchmark BED) at the checkpoint's call rule (SNV rule 2; INDEL ppr8, ONT pp).
# Sets = those of the six-tumor models with the features: SNV Illumina / PacBio on the AF 0.05 builds of the four CASTLE
# lines and COLO829T (HG008T at 0.08), SNV ONT and every INDEL model on the AF 0.08 sets.
# Runs: $PANSOMA_RUNS/m6_nositef_graphs_20261118/multicancer6_<P>_<KIND>_{d9,d46}_nositef
# A model whose tensors are not all on storage2 is skipped with the list of what is missing (to_ris.sh copies them from
# the old cluster); a model whose run directory exists is skipped too, so the script can be run again as data arrive.
# Run from the repository root on the login node:
#   bash analysis/m6_nositef_graphs_20261118/submit.sh            # every model that is ready
#   KINDS=SNV GRAPHS=d9 PLATFORMS="Illumina ONT" bash analysis/m6_nositef_graphs_20261118/submit.sh
# SLOTS (default 8): slots of each GPU job; exec-419 has 48, so 8 jobs at once need SLOTS=5 or fewer.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
A=/storage2/fs1/epigenome/Active/jiawei.s
D=$A/Dataset/pansoma_v2_tensors RUNS=${PANSOMA_RUNS:-$A/pansoma_net_v2_runs}
B=m6_nositef_graphs_20261118
HERE=../analysis/$B

sets() {  # PLATFORM KIND GRAPH -> HCC1395T HCC1954T H1437T H2009T HG008T COLO829T sets
  local P=$1 K=$2 G=$3 cp=$1 af="" hg co l
  [ "$P" = PacBio ] && cp=fiberseq
  [ "$K" = SNV ] && [ "$P" != ONT ] && af=_snvaf05          # AF 0.05 builds (CASTLE lines, COLO829T)
  if [ "$G" = d9 ]; then hg=HG008T_$P co=COLO829T_$cp$af; else hg=HG008T_${P}_d46 co=COLO829T_${cp}_d46$af; fi
  for l in HCC1395T HCC1954T H1437T H2009T; do printf '%s ' "${l}_${P}_$G$af"; done
  echo "$hg $co"
}

for K in ${KINDS:-SNV INDEL}; do for G in ${GRAPHS:-d9 d46}; do for P in ${PLATFORMS:-Illumina PacBio ONT}; do
  name=multicancer6_${P}_${K}_${G}_nositef job=pnv2_m6_${P}_${K}_${G}_nositef
  read -r -a S <<< "$(sets "$P" "$K" "$G")"
  if [ -e "$RUNS/$B/$name" ]; then echo "$name: $RUNS/$B/$name exists, skipped"; continue; fi
  if bjobs -J "$job" 2>/dev/null | grep -q "$job"; then echo "$name: job $job is in LSF, skipped"; continue; fi
  missing=() build=() wait=""
  for s in "${S[@]}"; do
    # every shard of the manifest at its full size (ready.py), so a set that is still being copied does not train
    if ! r=$(python3 analysis/$B/ready.py "$D/$s/tensors" "$K"); then missing+=("$s/tensors/$K ($r)"); continue; fi
    if [ ! -s "$RUNS/site_features/$s.$K.npz" ]; then       # PoN tags of the set (a link to the old cluster counts as missing)
      if [ -L "$RUNS/site_features/$s.$K.npz" ]; then missing+=("site_features/$s.$K.npz"); else build+=("$s"); fi
    fi
  done
  if [ ${#missing[@]} -gt 0 ]; then echo "$name: not on storage2, skipped: ${missing[*]}"; continue; fi
  for s in ${build[@]+"${build[@]}"}; do
    KIND=$K ris/submit.sh -g 0 -n 8 -m 48 "pnv2_sitef${K}_$s" site_features.sh "$s"
    wait="${wait}done(pnv2_sitef${K}_$s) && "
  done
  case $P in Illumina) mem=384;; *) mem=256;; esac
  KIND=$K NOSITEF=1 PLATFORM=$P; source ris/recipes.sh
  held="${S[0]}:chr1+chr20,${S[4]}:chr1,${S[5]}:chr1"
  ris/submit.sh ${wait:+-w "${wait% && }"} -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -g 1 -n "${SLOTS:-8}" -m $mem "$job" \
      "$HERE/train_predict.sh" "$B/$name" "$K" "$(IFS=,; echo "${S[*]}")" "$held" "${RECIPE[@]}" --seed 0 \
      --read-threads 8 --val-workers 6
  for s in "${S[0]}" "${S[4]}" "${S[5]}"; do
    ris/submit.sh -g 0 -n 4 -m 32 -w "done($job)" "pnv2_vcf_m6_${P}_${K}_${G}_nositef_$s" vcf.sh --set "$s" "$B/$name" "$K"
  done
done; done; done
bjobs -w 2>/dev/null | grep -c nositef | sed 's/^/jobs in LSF with "nositef" in the name: /'
echo "logs: $RUNS/logs/pnv2_m6_<P>_<KIND>_<GRAPH>_nositef.log, training: $RUNS/$B/<name>/train.log"

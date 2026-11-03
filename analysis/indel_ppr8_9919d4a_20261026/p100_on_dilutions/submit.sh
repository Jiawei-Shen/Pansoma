#!/bin/bash
# COLO829T P100 single-set models (d9, seeds 0-2) on chr1 of the P75 / P50 / P25 in-silico mixtures: one GPU job per
# model (predict_sets.sh -> <run>/test_chr1_<SET>/), then vcf.sh --set per mixture (INDEL also at the post-PoN best-F1
# threshold, vcf_chr1_<SET>_postpon/).
#   INDEL: COLO829_<P>_INDEL_d9_ppr8{,_s1,_s2} (the ppr8 recipe)
#   SNV: COLO829_<P>_SNV_d9_normbq{,_s1,_s2}; Illumina at the epochs the seed table reports (6 / 5 / 6, the user
#        2026-10-20): COLO829_Illumina_SNV_d9_normbq{_sweep_e06,_s1_sweep_e05,_s2_sweep_e06}
# Run from the repo root: bash analysis/indel_ppr8_9919d4a_20261026/p100_on_dilutions/submit.sh
set -euo pipefail
S=../analysis/indel_ppr8_9919d4a_20261026/p100_on_dilutions/predict_sets.sh
for P in Illumina fiberseq ONT; do
  sets=(COLO829T_${P}_P75 COLO829T_${P}_P50 COLO829T_${P}_P25)
  if [ $P = Illumina ]; then snv=(COLO829_Illumina_SNV_d9_normbq_sweep_e06 COLO829_Illumina_SNV_d9_normbq_s1_sweep_e05
                                  COLO829_Illumina_SNV_d9_normbq_s2_sweep_e06)
  else snv=(COLO829_${P}_SNV_d9_normbq COLO829_${P}_SNV_d9_normbq_s1 COLO829_${P}_SNV_d9_normbq_s2); fi
  indel=(COLO829_${P}_INDEL_d9_ppr8 COLO829_${P}_INDEL_d9_ppr8_s1 COLO829_${P}_INDEL_d9_ppr8_s2)
  for K in INDEL SNV; do
    if [ $K = INDEL ]; then runs=("${indel[@]}"); else runs=("${snv[@]}"); fi
    for name in "${runs[@]}"; do
      job=pnv2_p100pred_$name
      ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -n 4 -m 64 "$job" "$S" "$name" $K "${sets[@]}"
      for set in "${sets[@]}"; do
        ris/submit.sh -g 0 -n 4 -m 16 -w "done($job)" pnv2_vcf_${name}_$set vcf.sh --set "$set" "$name" $K
        if [ $K = INDEL ]; then
          ris/submit.sh -g 0 -n 4 -m 16 -w "done($job)" pnv2_vcfpp_${name}_$set vcf.sh --set "$set" "$name" INDEL postpon _postpon
        fi
      done
    done
  done
done

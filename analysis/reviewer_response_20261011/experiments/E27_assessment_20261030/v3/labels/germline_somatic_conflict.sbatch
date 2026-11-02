#!/bin/bash
#SBATCH -J assess_gs_conflict
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 2:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/v3/labels/logs/germline_somatic_conflict-%j.out
# R1.8: how often is a tensor's representative allele an exact germline truth allele (label 2, step 3) while another
# allele of its site is a somatic truth allele (labels.ndjson "somatic" list non-empty)? Under classify() the exact
# germline match wins (step 3 before the partial somatic steps 4 / 9). chr1-22, six HPRC v2.1 d46 sets.
set -euo pipefail
T=/scratch/jshen/data/pansoma_v2_tensors
printf "set\tkind\texact_germline\texact_germline_with_somatic_site_mate\texact_somatic_with_germline_site_mate\n"
for s in HG008T_Illumina HG008T_PacBio HG008T_ONT COLO829T_Illumina COLO829T_fiberseq COLO829T_ONT; do
  for k in SNV INDEL; do
    d=$T/${s}_d46/tensors/$k
    g=$(cat $d/chr*_labels.ndjson | grep -F '"reason": "representative_allele_is_germline_truth' | tee >(grep -cF '"somatic": [{' > /tmp/gsc.$$.$s.$k) | wc -l)
    sleep 1; gs=$(cat /tmp/gsc.$$.$s.$k); rm -f /tmp/gsc.$$.$s.$k
    sg=$(cat $d/chr*_labels.ndjson | grep -F '"reason": "representative_allele_is_somatic_truth"' | grep -vcF '"germline": []' || true)
    printf "%s\t%s\t%s\t%s\t%s\n" $s $k $g $gs $sg
  done
done

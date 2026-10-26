#!/bin/bash
# Share of chr1 truth positions with another truth position within 50 bp (all records, not restricted to the BED):
# HapMap positives come from donor haplotypes, so nearby positives travel on the same reads; somatic truth rarely does.
# usage: bash neighbours.sh
ml load bcftools 2>/dev/null
for v in /scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/HapMap_somatic_snv_indel.vcf.gz \
         /scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz; do
  bcftools query -r chr1 -f '%POS\n' $v | sort -n -u | awk -v name=$(basename $v) '{p[NR]=$1} END {n=NR; c=0;
    for (i=1;i<=n;i++) if ((i>1 && p[i]-p[i-1]<=50) || (i<n && p[i+1]-p[i]<=50)) c++;
    printf "%s: %d positions on chr1, %.3f with another truth position within 50 bp\n", name, n, c/n}'
done

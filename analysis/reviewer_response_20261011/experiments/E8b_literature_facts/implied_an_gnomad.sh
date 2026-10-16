#!/bin/bash
# Implied allele number (AN = AC / AF) of the GATK af-only-gnomad.hg38 resource, on chr1:1-60,000,000.
# The file carries only AC and AF; max implied AN separates gnomAD v2.0.x (max AN 277,264) from v2.1 (282,912).
module load bcftools/1.18
F=/scratch/jshen/data/Pansoma/panel_of_normal_VCFs/af-only-gnomad.hg38.vcf.gz
bcftools query -r chr1:1-60000000 -f '%INFO/AC\t%INFO/AF\n' "$F" \
 | awk -F'\t' '$1>=2000 && $2>0 {an=$1/$2; n++; if(an>m)m=an; b=int(an/1000)*1000; h[b]++} END {printf "records_AC_ge_2000\t%d\nmax_implied_AN\t%.0f\n", n, m; for (k in h) print "bin_"k"\t"h[k]}' \
 | sort -k1,1V

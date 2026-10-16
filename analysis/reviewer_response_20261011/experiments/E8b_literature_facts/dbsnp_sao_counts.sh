#!/bin/bash
# Count dbSNP 138 (GATK hg38) records by SAO (variant allele origin): 0 unspecified, 1 germline, 2 somatic, 3 both.
module load bcftools/1.18
F=/scratch/jshen/data/Pansoma/panel_of_normal_VCFs/Homo_sapiens_assembly38.dbsnp138.vcf.gz
bcftools query -f '%INFO/SAO\n' "$F" | sort | uniq -c | awk '{print "SAO="$2"\t"$1}'

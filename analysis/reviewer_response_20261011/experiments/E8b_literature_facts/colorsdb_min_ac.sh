#!/bin/bash
# CoLoRSdb v1.1.0 on chr20: number of records, records with AC = 1, minimum AC (-> colorsdb_chr20_min_ac.txt).
module load bcftools/1.18
F=/scratch/jshen/data/Pansoma/panel_of_normal_VCFs/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz
bcftools query -r chr20 -f '%INFO/AC\n' "$F" \
 | awk '{n++; if($1+0==1)a1++; if(min==""||$1+0<min)min=$1+0} END{print "chr20_records",n,"AC_eq_1",a1+0,"min_AC",min}'

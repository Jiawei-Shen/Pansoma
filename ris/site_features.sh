#!/bin/bash
# The site features of a set's INDEL tensors (machine_learning/pansoma_net_v2/site_features.py build) into
# $SITE_FEATURES/<SET>.INDEL.npz, the directory the INDEL recipe's --site-features reads (recipes.sh). bcftools is not
# in the image: $BCFTOOLS (env.sh) is bcftools 1.21 from bioconda in $PANSOMA_RUNS/env/bcftools.
# Usage: site_features.sh SET. A CPU job (one chromosome per slot):
#   ris/submit.sh -g 0 -n 8 -m 48 pnv2_sitef_HG008T_Illumina_d46 site_features.sh HG008T_Illumina_d46
set -euo pipefail
source "$(dirname "$0")/env.sh"
SET=$1
PON=$PANSOMA_DATA/panel_of_normal_VCFs
cd "$REPO/machine_learning"
time python -u -m pansoma_net_v2.site_features build --tensors "$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors" \
    --fasta "$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" \
    --gnomad "$PON/af-only-gnomad.hg38.vcf.gz" --colors "$PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz" \
    --out "$SITE_FEATURES" --cache-dir "$PANSOMA_RUNS/index_cache" --workers "${LSB_DJOB_NUMPROC:-8}" \
    --bcftools "$BCFTOOLS"

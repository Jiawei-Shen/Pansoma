#!/bin/bash
# The PoN + rtg vcfeval steps of vcf.sh only, on a run's existing GRCh38 calls (<run>/SRC/*.KIND.linear.vcf.gz), into
# <run>/OUT/eval_bed and eval_nobed: re-scores calls made earlier (an earlier tensor build, an earlier PoN rule) with
# this checkout's PoN script. As upstream's jobs/vcfeval_only.sh. Usage: vcfeval_only.sh RUN KIND SRC OUT. A CPU job:
#   ris/submit.sh -g 0 -n 4 -m 32 pnv2_rescore_INDEL_small vcfeval_only.sh HG008_Illumina_INDEL_small INDEL vcf_chr1 vcf_chr1_pon7fc7396
set -euo pipefail
source "$(dirname "$0")/env.sh"
RUN=$1 KIND=$2 SRC=$3 OUT=$4
D=$PANSOMA_RUNS/$RUN
SDF=$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf
PON=$PANSOMA_DATA/panel_of_normal_VCFs
T=$(python -c "import json, sys; print(json.load(open(sys.argv[1]))['tensors'][0])" "$D/args.json")
L=$(ls "$D/$SRC"/*."$KIND".linear.vcf.gz)
read -r TRUTH BED < <(python -c "import json, sys; t = json.load(open(sys.argv[1]))['truth']['somatic']; print(t['vcf'], t['bed'])" \
    "$T/$KIND/labels.manifest.json")
O=$D/$OUT; rm -rf "$O/eval_bed" "$O/eval_nobed"; mkdir -p "$O"
cd "$REPO/machine_learning"
python -u -m pansoma_net_v2.vcfeval --calls "$L" \
    --pon "$PON/af-only-gnomad.hg38.vcf.gz" "$PON/Homo_sapiens_assembly38.dbsnp138.vcf.gz" "$PON/1000g_pon.hg38.vcf.gz" \
          "$PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz" \
    --truth "$TRUTH" --bed "$BED" --sdf "$SDF" --rtg "$RTG" --output "$O/eval_bed"
python -u -m pansoma_net_v2.vcfeval --calls "$L" --pon-vcf "$O/eval_bed/$(basename "$L" .vcf.gz).pon.vcf.gz" \
    --truth "$TRUTH" --sdf "$SDF" --rtg "$RTG" --output "$O/eval_nobed"
cat "$O/eval_bed/report.txt" "$O/eval_nobed/report.txt"

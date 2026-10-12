#!/bin/bash
# A run's chr1 predictions -> graph VCF -> GRCh38 VCF -> PoN + rtg vcfeval against the truth VCF of the
# set's labels, inside the benchmark BED (eval_bed) and on the whole chromosome (eval_nobed), into <run>/vcf_chr1/.
# PASS is graph_vcf's default (the checkpoint's validation recall-0.9 threshold); THRESHOLD (a number, or
# "predict" for predict's own call) sets it, as needed for checkpoints trained before val.somatic_at_recall.
# Usage: vcf.sh RUN KIND [THRESHOLD]. A CPU job:
#   ris/submit.sh -g 0 -n 4 -m 32 pnv2_vcf vcf.sh HG008_Illumina_SNV_recipe SNV
#   ris/submit.sh -g 0 -n 4 -m 32 pnv2_vcf vcf.sh HG008_Illumina_SNV_test SNV predict
set -euo pipefail
source "$(dirname "$0")/env.sh"
RUN=$1 KIND=$2 THRESHOLD=${3:-}
D=$PANSOMA_RUNS/$RUN
G=$PANSOMA_DATA/pansoma_v2_tensors/graph_index
FASTA=$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
SDF=$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf
PON=$PANSOMA_DATA/panel_of_normal_VCFs
T=$(python -c "import json, sys; print(json.load(open(sys.argv[1]))['tensors'][0])" "$D/args.json")
P=$(ls "$D"/test_chr1/*."$KIND".predictions.ndjson.gz)
name=$(basename "$P" ."$KIND".predictions.ndjson.gz)
# the truth VCF and BED the labels were made with (old-cluster paths, bind-mounted by submit.sh)
read -r TRUTH BED < <(python -c "import json, sys; t = json.load(open(sys.argv[1]))['truth']['somatic']; print(t['vcf'], t['bed'])" \
    "$T/$KIND/labels.manifest.json")
O=$D/vcf_chr1; mkdir -p "$O"
cd "$REPO/machine_learning"
python -u -m pansoma_net_v2.graph_vcf --predictions "$P" --tensors "$T" --kind "$KIND" \
    --graph-index "$G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite" --output "$O/$name.$KIND.graph.vcf.gz" \
    ${THRESHOLD:+--threshold "$THRESHOLD"}
python -u -m pansoma_net_v2.linear_vcf --graph-vcf "$O/$name.$KIND.graph.vcf.gz" \
    --reference-path "$G/hprc-v1.1-mc-grch38.d9.grch38_path" --fasta "$FASTA" --output "$O/$name.$KIND.linear.vcf.gz"
rm -rf "$O/eval_bed" "$O/eval_nobed"
# the PoN for both kinds: filter_panel_of_normals.py applies its INDEL rule (gnomAD / CoLoRSdb AF >= 0.01) to INDELs
pon=(--pon "$PON/af-only-gnomad.hg38.vcf.gz" "$PON/Homo_sapiens_assembly38.dbsnp138.vcf.gz" "$PON/1000g_pon.hg38.vcf.gz"
     "$PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")
pon_vcf=(--pon-vcf "$O/eval_bed/$name.$KIND.linear.pon.vcf.gz")
python -u -m pansoma_net_v2.vcfeval --calls "$O/$name.$KIND.linear.vcf.gz" "${pon[@]}" \
    --truth "$TRUTH" --bed "$BED" --sdf "$SDF" --rtg "$RTG" --output "$O/eval_bed"
python -u -m pansoma_net_v2.vcfeval --calls "$O/$name.$KIND.linear.vcf.gz" "${pon_vcf[@]}" \
    --truth "$TRUTH" --sdf "$SDF" --rtg "$RTG" --output "$O/eval_nobed"
cat "$O/eval_bed/report.txt" "$O/eval_nobed/report.txt"

#!/bin/bash
# A run's chr1 predictions -> graph VCF -> GRCh38 VCF -> PoN + rtg vcfeval against the truth VCF of the
# set's labels, inside the benchmark BED (eval_bed) and on the whole chromosome (eval_nobed), into <run>/vcf_chr1/.
# PASS is graph_vcf's default (the checkpoint's validation recall-0.9 threshold); THRESHOLD (a number, or
# "predict" for predict's own call) overrides it.
# --set SET (a joint_run.sh run): the predictions <run>/test_chr1_<SET>/ of that set, into <run>/vcf_chr1_<SET>/.
# The graph index and GRCh38 reference path are those the set was built and labelled with (manifest.json graph_index,
# labels.manifest.json reference_path), by name in Dataset/pansoma_v2_tensors/graph_index/ (HPRC v1.1 d9, HPRC v2.1 d46
# or the GRCh38-only graph).
# Usage: vcf.sh [--set SET] RUN KIND [THRESHOLD]. A CPU job:
#   ris/submit.sh -g 0 -n 4 -m 32 pnv2_vcf vcf.sh HG008_Illumina_SNV_recipe SNV
#   ris/submit.sh -g 0 -n 4 -m 32 pnv2_vcf vcf.sh HG008_Illumina_SNV_test SNV predict
#   ris/submit.sh -g 0 -n 4 -m 32 pnv2_vcf vcf.sh --set COLO829T_Illumina_P50 COLO829_Illumina_SNV_joint SNV
set -euo pipefail
source "$(dirname "$0")/env.sh"
SET=""; if [ "${1:-}" = --set ]; then SET=$2; shift 2; fi
RUN=$1 KIND=$2 THRESHOLD=${3:-}
D=$PANSOMA_RUNS/$RUN
G=$PANSOMA_DATA/pansoma_v2_tensors/graph_index
FASTA=$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
SDF=$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf
PON=$PANSOMA_DATA/panel_of_normal_VCFs
if [ -n "$SET" ]; then
  T=$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors PRED=$D/test_chr1_$SET O=$D/vcf_chr1_$SET
else
  T=$(python -c "import json, sys; print(json.load(open(sys.argv[1]))['tensors'][0])" "$D/args.json") PRED=$D/test_chr1 O=$D/vcf_chr1
fi
P=$(ls "$PRED"/*."$KIND".predictions.ndjson.gz)
read -r GI RP < <(python -c "import json, os, sys; print(os.path.basename(json.load(open(sys.argv[1]))['graph_index']['path']),
os.path.basename(json.load(open(sys.argv[2]))['reference_path']))" "$T/$KIND/manifest.json" "$T/$KIND/labels.manifest.json")
[ -f "$G/$GI" ] && [ -d "$G/$RP" ] || { echo "missing $G/$GI or $G/$RP (copy the set's graph index there)"; exit 1; }
name=$(basename "$P" ."$KIND".predictions.ndjson.gz)
# the truth VCF and BED the labels were made with (old-cluster paths, bind-mounted by submit.sh)
read -r TRUTH BED < <(python -c "import json, sys; t = json.load(open(sys.argv[1]))['truth']['somatic']; print(t['vcf'], t['bed'])" \
    "$T/$KIND/labels.manifest.json")
mkdir -p "$O"
cd "$REPO/machine_learning"
python -u -m pansoma_net_v2.graph_vcf --predictions "$P" --tensors "$T" --kind "$KIND" \
    --graph-index "$G/$GI" --output "$O/$name.$KIND.graph.vcf.gz" \
    ${THRESHOLD:+--threshold "$THRESHOLD"}
python -u -m pansoma_net_v2.linear_vcf --graph-vcf "$O/$name.$KIND.graph.vcf.gz" \
    --reference-path "$G/$RP" --fasta "$FASTA" --output "$O/$name.$KIND.linear.vcf.gz"
rm -rf "$O/eval_bed" "$O/eval_nobed"
# the PoN for both kinds: filter_panel_of_normals.py applies its INDEL rule to INDELs (gnomAD / CoLoRSdb exact allele,
# AF >= 0.01 in homopolymers >= 8 bp, else >= 0.0001; 71cb51a)
pon=(--pon "$PON/af-only-gnomad.hg38.vcf.gz" "$PON/Homo_sapiens_assembly38.dbsnp138.vcf.gz" "$PON/1000g_pon.hg38.vcf.gz"
     "$PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")
pon_vcf=(--pon-vcf "$O/eval_bed/$name.$KIND.linear.pon.vcf.gz")
python -u -m pansoma_net_v2.vcfeval --calls "$O/$name.$KIND.linear.vcf.gz" "${pon[@]}" \
    --truth "$TRUTH" --bed "$BED" --sdf "$SDF" --rtg "$RTG" --output "$O/eval_bed"
python -u -m pansoma_net_v2.vcfeval --calls "$O/$name.$KIND.linear.vcf.gz" "${pon_vcf[@]}" \
    --truth "$TRUTH" --sdf "$SDF" --rtg "$RTG" --output "$O/eval_nobed"
cat "$O/eval_bed/report.txt" "$O/eval_nobed/report.txt"

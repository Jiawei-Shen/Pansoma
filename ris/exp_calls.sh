#!/bin/bash
# Experiment CPU job for one set of an exp_train.sh run: calls of chr2 and chr1 with the default decision (the
# checkpoint's validation recall-0.9 threshold, the PoN script's rules), then two decisions chosen on chr2 (the
# training data's held-out chromosome) and applied once to chr1:
#  - general: the threshold only (best call-level F1 on chr2, tools/tune_indel.py with one flat 0.01 INDEL PoN floor),
#    then chr1 with the PoN script's own rules; no read-AF floor.
#    The same rule for every set: what a model would carry to an unlabelled sample.
#  - tuned: threshold, a flat INDEL PoN AF floor (tools/pon_tag.py) and read-AF floor (tools/tune_indel.py grid):
#    per-set, for reference only.
# Usage: exp_calls.sh NAME KIND SET
# Output: $PANSOMA_RUNS/NAME/vcf_<SET>_<chrom>[_general|_tuned]/, tune_<SET>[.general].json.
set -euo pipefail
shopt -s inherit_errexit  # a failing command inside $(...) stops the job
source "$(dirname "$0")/env.sh"
NAME=$1 KIND=$2 SET=$3
D=$PANSOMA_RUNS/$NAME T=$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors
G=$PANSOMA_DATA/pansoma_v2_tensors/graph_index
FASTA=$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
SDF=$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf
P=$PANSOMA_DATA/panel_of_normal_VCFs
PONS=("$P/af-only-gnomad.hg38.vcf.gz" "$P/Homo_sapiens_assembly38.dbsnp138.vcf.gz" "$P/1000g_pon.hg38.vcf.gz"
      "$P/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")
read -r TRUTH BED < <(python -c "import json, sys; t = json.load(open(sys.argv[1]))['truth']['somatic']; print(t['vcf'], t['bed'])" \
    "$T/$KIND/labels.manifest.json")
cd "$REPO/machine_learning"
export PYTHONPATH="$REPO/machine_learning:$REPO/scripts:$REPO/ris/tools${PYTHONPATH:+:$PYTHONPATH}"
calls() {  # PRED_DIR OUT_DIR [graph_vcf options]: graph VCF -> GRCh38 VCF
  local pred=$1 o=$2; shift 2
  local pf; pf=$(ls "$pred"/*."$KIND".predictions.ndjson.gz); mkdir -p "$o"
  local n; n=$(basename "$pf" ."$KIND".predictions.ndjson.gz)
  python -u -m pansoma_net_v2.graph_vcf --predictions "$pf" --tensors "$T" --kind "$KIND" \
    --graph-index "$G/hprc-v1.1-mc-grch38.d9.graph_index.sqlite" --output "$o/$n.$KIND.graph.vcf.gz" "$@" >&2
  python -u -m pansoma_net_v2.linear_vcf --graph-vcf "$o/$n.$KIND.graph.vcf.gz" \
    --reference-path "$G/hprc-v1.1-mc-grch38.d9.grch38_path" --fasta "$FASTA" --output "$o/$n.$KIND.linear.vcf.gz" >&2
  echo "$o/$n.$KIND.linear.vcf.gz"
}
evals() {  # LINEAR_VCF OUT_DIR PON_ARGS...: vcfeval in the BED, then without it on the same PoN tags
  local lin=$1 o=$2; shift 2
  rm -rf "$o/eval_bed" "$o/eval_nobed"
  python -u -m pansoma_net_v2.vcfeval --calls "$lin" "$@" --truth "$TRUTH" --bed "$BED" --sdf "$SDF" --rtg "$RTG" --output "$o/eval_bed"
  python -u -m pansoma_net_v2.vcfeval --calls "$lin" --pon-vcf "$o/eval_bed/$(basename "$lin" .vcf.gz).pon.vcf.gz" \
    --truth "$TRUTH" --sdf "$SDF" --rtg "$RTG" --output "$o/eval_nobed"
}
for c in chr2 chr1; do  # the default decision
  lin=$(calls "$D/pred_${SET}_$c" "$D/vcf_${SET}_$c" | tail -1)
  evals "$lin" "$D/vcf_${SET}_$c" --pon "${PONS[@]}"
done
# general: threshold only, chosen on chr2
python -u "$REPO/ris/tools/tune_indel.py" "$D/vcf_${SET}_chr2/eval_bed" "$D/tune_${SET}.general.json" --floors 0.01 --read-af 0
GTHR=$(python -c "import json, sys; print(json.load(open(sys.argv[1]))['threshold'])" "$D/tune_${SET}.general.json")
lin=$(calls "$D/pred_${SET}_chr1" "$D/vcf_${SET}_chr1_general" --threshold "$GTHR" | tail -1)
evals "$lin" "$D/vcf_${SET}_chr1_general" --pon "${PONS[@]}"
# tuned per set (reference): threshold, PoN floor, read-AF floor
python -u "$REPO/ris/tools/tune_indel.py" "$D/vcf_${SET}_chr2/eval_bed" "$D/tune_${SET}.json"
read -r THR FLOOR RAF < <(python -c "import json, sys; j = json.load(open(sys.argv[1])); print(j['threshold'], j['pon_floor'], j['read_af_min'])" "$D/tune_${SET}.json")
o=$D/vcf_${SET}_chr1_tuned
lin=$(calls "$D/pred_${SET}_chr1" "$o" --threshold "$THR" --min-af "$RAF" | tail -1)
PONVCF=$o/$(basename "$lin" .vcf.gz).floor$FLOOR.pon.vcf.gz
python -u "$REPO/ris/tools/pon_tag.py" "$FLOOR" "$lin" "$PONVCF" "${PONS[@]}"
rm -rf "$o/eval_bed" "$o/eval_nobed"
python -u -m pansoma_net_v2.vcfeval --calls "$lin" --pon-vcf "$PONVCF" --truth "$TRUTH" --bed "$BED" --sdf "$SDF" --rtg "$RTG" --output "$o/eval_bed"
python -u -m pansoma_net_v2.vcfeval --calls "$lin" --pon-vcf "$PONVCF" --truth "$TRUTH" --sdf "$SDF" --rtg "$RTG" --output "$o/eval_nobed"
echo "== $NAME $SET: general = threshold $GTHR from chr2; tuned = threshold $THR, PoN floor $FLOOR, read AF >= $RAF"
for x in vcf_${SET}_chr2 vcf_${SET}_chr1 vcf_${SET}_chr1_general vcf_${SET}_chr1_tuned; do echo "-- $x"; grep -E "^(raw|pon) " "$D/$x/eval_bed/report.txt"; done

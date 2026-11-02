#!/usr/bin/env bash
# Pansoma minimal example: sorted GAM -> target nodes -> tensors -> SNV + INDEL models -> GRCh38 VCF -> PoN filter
# (-> rtg vcfeval against the truth, when rtg is on PATH).
#
# Input (DATA, ~300 MB): HG008-T Illumina WGS reads (GIAB, BCM 2024-03-13) that vg giraffe aligned to the
# HPRC v1.1 d9 graph, the records visiting GRCh38 chr1:242,498,000-243,002,000 (vg 1.77 gamsort -i); the graph
# index, GRCh38 reference path and chromosome table of that graph cut to the visited nodes; GRCh38 chr1 (no-alt
# analysis set, N outside chr1:242.4-243.1 Mb); the four PoN VCFs and the HG008-T somatic benchmark v0.2 cut to
# the window; the HG008T Illumina SNV and INDEL models (trained on chr2-22; chr1 is their held-out test chromosome).
#
# Usage (from the repository root, in the environment of conda/pansoma-cpu.yml):
#   bash example/run_example.sh DATA OUT [THREADS]        # THREADS default 4; OUT must not exist
# Expected output: example/expected/ (check: python example/check_example.py OUT example/expected).
set -euo pipefail
DATA=$(cd "$1" && pwd); OUT=$2; THREADS=${3:-4}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
[ -e "$OUT" ] && { echo "$OUT exists" >&2; exit 1; }
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
GI=$DATA/graph/hprc-v1.1-mc-grch38.d9.chr1_242.5Mb.graph_index.sqlite
RP=$DATA/graph/hprc-v1.1-mc-grch38.d9.chr1_242.5Mb.grch38_path
CI=$DATA/graph/hprc-v1.1-mc-grch38.d9.chr_node_ranges.tsv
GAM=$DATA/gam/HG008T_Illumina.chr1_242.5-243.0Mb.sorted.gam
FA=$DATA/reference/GRCh38_no_alt.chr1_masked_outside_242.4-243.1Mb.fa.gz
PON=$DATA/pon; TRUTH=$DATA/truth
PONS="$PON/af-only-gnomad.hg38.chr1_242.4-243.1Mb.vcf.gz $PON/Homo_sapiens_assembly38.dbsnp138.chr1_242.4-243.1Mb.vcf.gz
      $PON/1000g_pon.hg38.chr1_242.4-243.1Mb.vcf.gz $PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.chr1_242.4-243.1Mb.vcf.gz"
S=$OUT/HG008T_Illumina            # sample directory: discovery/ run/ tensors/
T=$S/tensors
cd "$REPO"
t() { local t0=$SECONDS; echo "== $1"; shift; "$@"; echo "   ($(( SECONDS - t0 )) s)"; }

# 0. the native GAM decoder (C++17 + pybind11), once per checkout and Python version; the builder's output is identical
#    with the Python decoder (PANSOMA_DECODER=python), which is ~28x slower.
#    `run discover` needs it in any case.
if ! python -c "import sys; from indexed_gam_pipeline_v4 import native; sys.exit(native.load()[0] is None)"; then
  t "compile the native decoder" python -m indexed_gam_pipeline_v4.native compile
fi

# 1. target nodes (the genome run's rule: > 5 % of the MAPQ > 5 mappings carry an edit)
t "discover" python -m indexed_gam_pipeline_v4.run discover --gam "$GAM" --output "$S/discovery" \
  --graph-index "$GI" --processes "$THREADS"

# 2. tensors with the settings of the HG008T Illumina genome run, merged per chromosome (no labels)
t "prepare" python -m indexed_gam_pipeline_v4.orchestrate prepare --root "$S/run" --tensors "$T" --gam "$GAM" \
  --nodes "$S/discovery/target_nodes.txt" --node-stats "$S/discovery/node_stats.json" --graph-index "$GI" \
  --haplotypes 90 --tasks "$THREADS" --processes "$THREADS" --gam-cache-mb 1024 \
  --snv-min-af 0.08 --indel-min-af 0.08 --min-variants 3 --max-node-reads 800 \
  --chromosomes autosome --chr-index "$CI" --merge-shard-size 32768 --reference-path "$RP"
t "build + merge" bash "$S/run/run.sh"

# 3. the models (CPU or one GPU), the graph VCF and the GRCh38 VCF
cd "$REPO/machine_learning"
for KIND in SNV INDEL; do
  NAME=HG008T_Illumina.tensors.$KIND
  t "predict $KIND" python -m pansoma_net_v2.predict --checkpoint "$DATA/models/HG008_Illumina_${KIND}_d9_normbq.pth" \
    --tensors "$T" --kinds "$KIND" --output "$OUT/predict_$KIND" --num-workers "$THREADS"
  t "graph_vcf $KIND" python -m pansoma_net_v2.graph_vcf --predictions "$OUT/predict_$KIND/$NAME.predictions.ndjson.gz" \
    --tensors "$T" --kind "$KIND" --graph-index "$GI" --output "$OUT/vcf/$NAME.graph.vcf.gz"
  t "linear_vcf $KIND" python -m pansoma_net_v2.linear_vcf --graph-vcf "$OUT/vcf/$NAME.graph.vcf.gz" \
    --reference-path "$RP" --fasta "$FA" --output "$OUT/vcf/$NAME.linear.vcf.gz"
  # 4. panel-of-normals filter (README section 5): matched records get FILTER=PanelOfNormals
  t "PoN $KIND" python -u "$REPO/scripts/filter_panel_of_normals.py" "$OUT/vcf/$NAME.linear.vcf.gz" \
    "$OUT/vcf/$NAME.pon.vcf.gz" --pon $PONS
done

# 5. optional: rtg vcfeval against the HG008-T somatic benchmark v0.2 inside its BED (PASS after the PoN)
if command -v rtg >/dev/null; then
  t "rtg format" rtg format -o "$OUT/reference.sdf" "$FA" >/dev/null
  for KIND in SNV INDEL; do
    NAME=HG008T_Illumina.tensors.$KIND
    t "vcfeval $KIND" python -m pansoma_net_v2.vcfeval --calls "$OUT/vcf/$NAME.linear.vcf.gz" \
      --pon $PONS \
      --truth "$TRUTH/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.chr1_242.4-243.1Mb.vcf.gz" \
      --bed "$TRUTH/HG008-T_somatic_smvar_benchmark_v0.2_all.chr1_242.5-243.0Mb.bed" --sdf "$OUT/reference.sdf" \
      --rtg rtg --threads "$THREADS" --output "$OUT/eval_$KIND"
  done
fi
echo "done: $OUT/vcf/*.pon.vcf.gz"

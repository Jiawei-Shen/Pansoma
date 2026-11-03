#!/bin/bash
# tmp/snv_sitefeat_20261103 copy of machine_learning/pansoma_net_v2/jobs/vcf.sh that runs the code/ copy (graph_vcf takes the checkpoint call rule).
# One run's chr1 predictions -> graph VCF -> GRCh38 VCF -> PoN -> rtg vcfeval, inside the benchmark BED and without it.
# usage: [THRESHOLD=t] [SET_NAME=set] [TEST=test_chr1] [CHM13=1] sbatch -J NAME -o LOG vcf.sh RUN KIND [OUT]   (RUN under
# pansoma_net_v2_runs; OUT default <run>/vcf_chr1). SET_NAME picks one set's predictions when the run predicted several
# sets (e.g. COLO829T_fiberseq_P25); TEST the predictions directory of the run (e.g. test_chr20; the chromosomes are
# the predicted ones).
# The set: pansoma_v2_tensors/<set>/tensors, <set> the predictions file's prefix (any graph: d9, d46, linear). Truth VCF,
# BED, FASTA and reference path: the set's labels (labels.manifest.json); graph index: the set's manifest.json.
# PASS: p_somatic >= the checkpoint's validation threshold at recall 0.9 (graph_vcf's default); THRESHOLD=t instead
# (checkpoints from before val.somatic_at_recall). AF floor: the checkpoint's training labels.
# PoN (scripts/filter_panel_of_normals.py: allele match; SNV gnomAD/CoLoRSdb AF >= 0.0001, dbSNP non-somatic, 1000G;
# INDEL gnomAD/CoLoRSdb only, AF >= 0.01 in homopolymers >= 8 bp, else >= 0.0001) for both kinds. SNV off-reference
# tensors: offref_call when predict ran the rescue (--offref-site90 10; off by default). eval_*/report.txt and report.json record the PASS result
# and the chr1 curve's best point (p_somatic, P, R, TP, FP, FN; test-tuned, reference only); summary prints both.
# Optional, CHM13=1: also the CHM13 VCF <name>.chm13.vcf.gz (+ .chm13.unplaced.vcf.gz), linear_vcf on the graph's CHM13
# reference path (<graph>.chm13_path next to the GRCh38 one; graph_prep.sh with CHM13_FASTA) and the FASTA it was checked
# with. Not evaluated and without the PoN: the truth sets and the PoN VCFs are on GRCh38.
#SBATCH -p general
#SBATCH -c 4
#SBATCH --mem=3G
#SBATCH -t 2:00:00
set -euo pipefail
RUN=$1; KIND=$2
R=/scratch/jshen/data/pansoma_net_v2_runs; OUT=${3:-$R/$RUN/vcf_chr1}
D=/scratch/jshen/data/pansoma_v2_tensors
P=/wanglab/jshen/anaconda3/bin/python; RTG=/scratch/jshen/bin/rtg-tools-3.12.1/rtg
PON=/scratch/jshen/data/Pansoma/panel_of_normal_VCFs
PONS="$PON/af-only-gnomad.hg38.vcf.gz $PON/Homo_sapiens_assembly38.dbsnp138.vcf.gz $PON/1000g_pon.hg38.vcf.gz $PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"
PRED=$(ls $R/$RUN/${TEST:-test_chr1}/${SET_NAME:-*}.*.$KIND.predictions.ndjson.gz)
[ $(echo $PRED | wc -w) = 1 ] || { echo "need exactly one predictions file, got: $PRED (set SET_NAME)"; exit 1; }
# the set: <set>.tensors.<KIND>.predictions.ndjson.gz
B=$(basename $PRED); SET=$D/${B%%.*}/tensors
[ -d $SET/$KIND ] || { echo "no tensor set $SET/$KIND for $PRED"; exit 1; }
read TRUTH BED FASTA REFPATH GRAPH < <($P -c "
import json; m = json.load(open('$SET/$KIND/labels.manifest.json'))
print(m['truth']['somatic']['vcf'], m['truth']['somatic']['bed'], m['fasta'], m['reference_path'],
      json.load(open('$SET/$KIND/manifest.json'))['graph_index']['path'])")
if ! command -v java > /dev/null; then module load java; fi
echo "host $(hostname) job ${SLURM_JOB_ID:-}; $PRED; set $SET; graph $GRAPH; truth $TRUTH; java $(java -version 2>&1 | head -1)"
cd /scratch/jshen/Github/Pansoma/tmp/snv_sitefeat_20261103/code/machine_learning
mkdir -p $OUT
NAME=$(basename $(dirname $SET)).$KIND
T="/usr/bin/time -f %es_%MKB"
$T $P -m pansoma_net_v2.graph_vcf --predictions $PRED --tensors $SET --kind $KIND \
  --graph-index $GRAPH --output $OUT/$NAME.graph.vcf.gz ${THRESHOLD:+--threshold $THRESHOLD}
$T $P -m pansoma_net_v2.linear_vcf --graph-vcf $OUT/$NAME.graph.vcf.gz --reference-path $REFPATH --fasta $FASTA \
  --output $OUT/$NAME.linear.vcf.gz
if [ -n "${CHM13:-}" ]; then
  CP=${REFPATH%.grch38_path}.chm13_path
  CF=$($P -c "import json; print(json.load(open('$CP/meta.json'))['checks']['graph_index_and_fasta']['fasta'])")
  $T $P -m pansoma_net_v2.linear_vcf --graph-vcf $OUT/$NAME.graph.vcf.gz --reference-path $CP --fasta $CF \
    --output $OUT/$NAME.chm13.vcf.gz
fi
EVAL="$P -m pansoma_net_v2.vcfeval --calls $OUT/$NAME.linear.vcf.gz --truth $TRUTH --sdf /scratch/jshen/REF/GRCh38.primary.sdf --rtg $RTG --threads 4"
PON_BED="--pon $PONS"; PON_NOBED="--pon-vcf $OUT/eval_bed/$NAME.linear.pon.vcf.gz"
echo "== vcfeval inside the BED"
$T $EVAL $PON_BED --bed $BED --output $OUT/eval_bed
echo "== vcfeval on the whole chromosomes (the PoN-tagged calls of eval_bed)"
$T $EVAL $PON_NOBED --output $OUT/eval_nobed

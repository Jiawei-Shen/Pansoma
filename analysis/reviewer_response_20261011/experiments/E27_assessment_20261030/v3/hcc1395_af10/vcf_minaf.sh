#!/bin/bash
# machine_learning/pansoma_net_v2/jobs/vcf.sh with graph_vcf --min-af $MIN_AF (calls whose A1 allele fraction is below it
# are LowAF, not PASS), inside the benchmark BED only. Same predictions, threshold, PoN and rtg as vcf.sh.
# usage: MIN_AF=0.1 TEST=test_chr1_normfix SET_NAME=HCC1395T_Illumina_d46 sbatch vcf_minaf.sh RUN KIND OUT
#SBATCH -p general
#SBATCH -x maotai,guinness,meow,miao
#SBATCH -c 4
#SBATCH --mem=2G
#SBATCH -t 2:00:00
set -euo pipefail
RUN=$1; KIND=$2; OUT=$3
R=/scratch/jshen/data/pansoma_net_v2_runs; D=/scratch/jshen/data/pansoma_v2_tensors
P=/wanglab/jshen/anaconda3/bin/python; RTG=/scratch/jshen/bin/rtg-tools-3.12.1/rtg
PON=/scratch/jshen/data/Pansoma/panel_of_normal_VCFs
PONS="$PON/af-only-gnomad.hg38.vcf.gz $PON/Homo_sapiens_assembly38.dbsnp138.vcf.gz $PON/1000g_pon.hg38.vcf.gz $PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"
PRED=$(ls $R/$RUN/${TEST:-test_chr1}/${SET_NAME:-*}.*.$KIND.predictions.ndjson.gz)
[ $(echo $PRED | wc -w) = 1 ] || { echo "need exactly one predictions file, got: $PRED"; exit 1; }
B=$(basename $PRED); SET=$D/${B%%.*}/tensors
read TRUTH BED FASTA REFPATH GRAPH < <($P -c "
import json; m = json.load(open('$SET/$KIND/labels.manifest.json'))
print(m['truth']['somatic']['vcf'], m['truth']['somatic']['bed'], m['fasta'], m['reference_path'],
      json.load(open('$SET/$KIND/manifest.json'))['graph_index']['path'])")
command -v java > /dev/null || module load java
echo "host $(hostname) job ${SLURM_JOB_ID:-}; $PRED; min AF $MIN_AF; truth $TRUTH"
cd /scratch/jshen/Github/Pansoma/machine_learning
mkdir -p $OUT
NAME=$(basename $(dirname $SET)).$KIND
$P -m pansoma_net_v2.graph_vcf --predictions $PRED --tensors $SET --kind $KIND --graph-index $GRAPH \
  --output $OUT/$NAME.graph.vcf.gz --min-af $MIN_AF
$P -m pansoma_net_v2.linear_vcf --graph-vcf $OUT/$NAME.graph.vcf.gz --reference-path $REFPATH --fasta $FASTA \
  --output $OUT/$NAME.linear.vcf.gz
$P -m pansoma_net_v2.vcfeval --calls $OUT/$NAME.linear.vcf.gz --truth $TRUTH --sdf /scratch/jshen/REF/GRCh38.primary.sdf \
  --rtg $RTG --threads 4 --pon $PONS --bed $BED --output $OUT/eval_bed

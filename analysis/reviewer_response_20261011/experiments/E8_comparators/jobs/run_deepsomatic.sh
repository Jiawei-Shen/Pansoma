#!/bin/bash
# E8: DeepSomatic 1.10.0 (google/deepsomatic:1.10.0, CPU, Singularity) tumor-only on chr1, then the Pansoma-protocol
# scoring (score_comparator.py) for SNV and INDEL.
# usage: sbatch run_deepsomatic.sh DATASET BAM MODEL_TYPE SET [extra run_deepsomatic options]
#SBATCH -p general
#SBATCH -c 16
#SBATCH --mem=80G
#SBATCH -t 3-00:00:00
set -euo pipefail
DATASET=$1; BAM=$2; MODEL=$3; SET=$4; shift 4
module load singularity/4.1.0 bcftools/1.18 samtools/1.18
if ! command -v java > /dev/null; then module load java; fi
SIF=/scratch/jshen/bin/deepsomatic_1.10.0.sif
REF=/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
E=/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators
X=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E8_comparators
P=/wanglab/jshen/anaconda3/bin/python
OUT=$E/calls/deepsomatic/$DATASET
mkdir -p $OUT/logs $OUT/intermediate
echo "host $(hostname) job ${SLURM_JOB_ID:-}; $(date); image $SIF ($(sha256sum $SIF | cut -c1-16)); bam $BAM; model $MODEL; set $SET; extra: $*"
HELP=$(singularity exec -B /scratch $SIF run_deepsomatic --helpfull 2>&1 || true)
singularity exec -B /scratch $SIF run_deepsomatic --version 2>&1 | head -3 || true
grep -q -- "$MODEL" <<< "$HELP" || { echo "model type $MODEL is not in run_deepsomatic --helpfull"; grep -A8 model_type <<< "$HELP"; exit 2; }
# the tumor-only recommendation of the DeepSomatic docs: the default panel of normals
PON=""; grep -q use_default_pon_filtering <<< "$HELP" && PON="--use_default_pon_filtering=true"
CMD="run_deepsomatic --model_type=$MODEL --ref=$REF --reads_tumor=$BAM --output_vcf=$OUT/output.vcf.gz --sample_name_tumor=$DATASET --num_shards=16 --logging_dir=$OUT/logs --intermediate_results_dir=$OUT/intermediate --regions=chr1 $PON $*"
echo "singularity exec -B /scratch $SIF $CMD" | tee $OUT/command.txt
/usr/bin/time -f "deepsomatic %es %MKB" singularity exec -B /scratch $SIF $CMD
ls -la $OUT
for KIND in SNV INDEL; do
  S=$E/score/deepsomatic/$DATASET/$KIND
  if [ -s $S/report.json ]; then echo "$S done"; continue; fi
  /usr/bin/time -f "score $KIND %es %MKB" $P $X/score_comparator.py --vcf $OUT/output.vcf.gz --set $SET --kind $KIND --caller deepsomatic_1.10.0 --output $S
done
date

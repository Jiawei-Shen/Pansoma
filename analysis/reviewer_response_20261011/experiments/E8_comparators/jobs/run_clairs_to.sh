#!/bin/bash
# E8: ClairS-TO v0.5.1 (hkubal/clairs-to:v0.5.1, Singularity) tumor-only on chr1, then the Pansoma-protocol scoring
# (score_comparator.py) for SNV and INDEL. usage: sbatch run_clairs_to.sh DATASET BAM PLATFORM SET [extra ClairS-TO options]
# DATASET: output name; SET: tensor set whose labels manifest gives truth/BED/FASTA (e.g. HG008T_PacBio).
#SBATCH -p general
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 2-00:00:00
set -euo pipefail
DATASET=$1; BAM=$2; PLATFORM=$3; SET=$4; shift 4
module load singularity/4.1.0 bcftools/1.18 samtools/1.18
if ! command -v java > /dev/null; then module load java; fi
SIF=/scratch/jshen/bin/clairs-to_v0.5.1.sif
REF=/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
E=/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators
X=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E8_comparators
P=/wanglab/jshen/anaconda3/bin/python
OUT=$E/calls/clairs_to/$DATASET
mkdir -p $OUT
echo "host $(hostname) job ${SLURM_JOB_ID:-}; $(date); image $SIF ($(sha256sum $SIF | cut -c1-16)); bam $BAM; platform $PLATFORM; set $SET; extra: $*"
singularity exec -B /scratch $SIF /opt/bin/run_clairs_to --version || true
CMD="/opt/bin/run_clairs_to --tumor_bam_fn $BAM --ref_fn $REF --threads 16 --platform $PLATFORM --output_dir $OUT --ctg_name chr1 --conda_prefix /opt/micromamba/envs/clairs-to $*"
echo "singularity exec -B /scratch $SIF $CMD" | tee $OUT/command.txt
/usr/bin/time -f "clairs-to %es %MKB" singularity exec -B /scratch $SIF $CMD
ls -la $OUT
VCFS=""
for v in snv.vcf.gz indel.vcf.gz; do [ -s $OUT/$v ] && VCFS="$VCFS --vcf $OUT/$v"; done
for KIND in SNV INDEL; do
  S=$E/score/clairs_to/$DATASET/$KIND
  if [ -s $S/report.json ]; then echo "$S done"; continue; fi
  /usr/bin/time -f "score $KIND %es %MKB" $P $X/score_comparator.py $VCFS --set $SET --kind $KIND --caller clairs_to_v0.5.1 --output $S
done
date

#!/bin/bash
#SBATCH -p general
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/clairs_to_nopon_rerun_20261102/logs/%x_%j.out
# ClairS-TO rerun without its panel-of-normals step (--disable_nonsomatic_tagging), everything else as the original E8
# run of the dataset (same version, platform, contigs, reference and BAM; threads do not change the 50 chunks per contig).
# A FILTER rewrite (NonSomatic -> PASS) is not equivalent: the tagging (step 3) runs before realignment, postfilter and
# haplotype filtering, which only process FILTER == PASS records. Then E8 score_comparator.py per kind and contig.
# usage: sbatch -J cto_nopon_<DS> run_nopon.sh DS VERSION BAM REF PLATFORM CONTIGS(comma) SET
set -euo pipefail
DS=$1; VER=$2; BAM=$3; REF=$4; PLAT=$5; CTGS=$6; SET=$7
module load singularity/4.1.0 bcftools/1.18
command -v java > /dev/null || module load java
SIF=/scratch/jshen/bin/clairs-to_v$VER.sif
OUT=/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/calls_nopon/clairs_to/$DS
SCORE=/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score_nopon_rerun/clairs_to
X=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E8_comparators
echo "host $(hostname) job $SLURM_JOB_ID; $(date)"
if [ ! -s $OUT/snv.vcf.gz ] || [ ! -s $OUT/indel.vcf.gz ]; then
  mkdir -p $OUT
  CMD="/opt/bin/run_clairs_to --tumor_bam_fn $BAM --ref_fn $REF --threads $SLURM_CPUS_PER_TASK --platform $PLAT --output_dir $OUT --ctg_name $CTGS --conda_prefix /opt/micromamba/envs/clairs-to --disable_nonsomatic_tagging"
  echo "ClairS-TO $VER (singularity exec --cleanenv -B /scratch $SIF) $CMD" > $OUT/command.txt; cat $OUT/command.txt
  /usr/bin/time -f "clairs-to %es %MKB" singularity exec --cleanenv -B /scratch $SIF $CMD
fi
for c in ${CTGS//,/ }; do
  d=$DS; [ "$c" = chr1 ] || d=${DS}_$c
  for k in SNV INDEL; do
    [ -s $SCORE/$d/$k/report.json ] && continue
    /wanglab/jshen/anaconda3/bin/python $X/score_comparator.py --vcf $OUT/snv.vcf.gz --vcf $OUT/indel.vcf.gz --set $SET --kind $k \
        --caller clairs_to --chrom $c --threads $SLURM_CPUS_PER_TASK --output $SCORE/$d/$k
  done
done
date

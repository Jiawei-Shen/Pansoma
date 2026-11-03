#!/bin/bash
#SBATCH -p general
#SBATCH -x maotai,guinness,meow,miao
#SBATCH -c 4
#SBATCH --mem=2G
#SBATCH -t 2:00:00
#SBATCH --array=0-5
# Comparator calls with the allele-fraction floor (rewrite_minaf.py), scored by E8 score_comparator.py exactly as the E8 rows.
set -euo pipefail
module load bcftools/1.18
command -v java > /dev/null || module load java
W=/scratch/jshen/Github/Pansoma/tmp/af10_hcc1395_20261103
X=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E8_comparators
D=/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score_minaf${MIN_AF}
IFS=$'\t' read -r CALLER DS KIND SET CHR VCFS < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $W/tasks.tsv)
echo "$CALLER $DS $KIND $SET $CHR min AF $MIN_AF host $(hostname) job ${SLURM_JOB_ID:-}"
O=$D/$CALLER/$DS/$KIND; I=$D/inputs/$CALLER/$DS/$KIND
mkdir -p $I
ARGS=()
for v in ${VCFS//,/ }; do
  b=$(basename $v .gz)
  /wanglab/jshen/anaconda3/bin/python $W/rewrite_minaf.py $CALLER $MIN_AF $v $I/$b.tmp && mv $I/$b.tmp $I/$b
  ARGS+=(--vcf $I/$b)
done
/wanglab/jshen/anaconda3/bin/python $X/score_comparator.py "${ARGS[@]}" --set $SET --kind $KIND --caller $CALLER \
    --chrom $CHR --threads 4 --output $O

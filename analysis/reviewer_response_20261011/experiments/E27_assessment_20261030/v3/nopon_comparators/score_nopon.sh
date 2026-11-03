#!/bin/bash
#SBATCH -p general
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH --array=0-35
# Comparator calls without the tool's own population filtering (rewrite_nopon.py), scored by E8 score_comparator.py
# exactly as the E8 / E12b comparator rows. Output: /scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score_nopon/
set -euo pipefail
module load bcftools/1.18
command -v java > /dev/null || module load java
W=/scratch/jshen/Github/Pansoma/tmp/nopon_comparators_20261102
X=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E8_comparators
D=/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score_nopon
IFS=$'\t' read -r CALLER DS KIND SET CHR VCFS < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $W/tasks.tsv)
echo "$CALLER $DS $KIND $SET $CHR host $(hostname) job ${SLURM_JOB_ID:-}"
O=$D/$CALLER/$DS/$KIND; I=$D/inputs/$CALLER/$DS/$KIND
[ -s $O/report.json ] && { echo "$O done"; exit 0; }
mkdir -p $I
ARGS=()
for v in ${VCFS//,/ }; do
  b=$(basename $v .gz)
  /wanglab/jshen/anaconda3/bin/python $W/rewrite_nopon.py $CALLER $v $I/$b.tmp && mv $I/$b.tmp $I/$b
  ARGS+=(--vcf $I/$b)
done
/usr/bin/time -f "score %es %MKB" /wanglab/jshen/anaconda3/bin/python $X/score_comparator.py "${ARGS[@]}" \
    --set $SET --kind $KIND --caller $CALLER --chrom $CHR --threads 4 --output $O

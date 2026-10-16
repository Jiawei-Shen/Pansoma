#!/bin/bash
# One E9a/E9e task: line SLURM_ARRAY_TASK_ID + 1 of the task list $1 ("RUN VCFDIR" under reviewer_cpu_20261013):
# projection (project_unplaced.py) then complete-output evaluation (evaluate_complete.py), then the run's coordination
# TSV <scratch>/<RUN>.projected.tsv (all its vcf dirs; under flock).
#SBATCH -p general
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
set -euo pipefail
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $1)
read RUN VDIR <<< "$LINE"
R=${E9A_RUNS:-/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013}
S=/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E9a_unplaced_projection
P=/wanglab/jshen/anaconda3/bin/python
O=${E9A_OUT:-$S/runs}
OUT=$O/$RUN/$VDIR
echo "task ${SLURM_ARRAY_TASK_ID}: $RUN $VDIR host $(hostname) job ${SLURM_JOB_ID:-}"
if ! command -v java > /dev/null; then module load java; fi
mkdir -p $OUT
[ -s $OUT/projection.json ] || /usr/bin/time -f "project %es %MKB" $P $E/project_unplaced.py --vcf-dir $R/$RUN/$VDIR --run $RUN --output $OUT
[ -s $OUT/metrics.json ] || /usr/bin/time -f "evaluate %es %MKB" $P $E/evaluate_complete.py --vcf-dir $R/$RUN/$VDIR --projection $OUT --threads 4
(
  flock 9
  TMP=$S/$RUN.projected.tsv.tmp.$$
  first=1
  for f in $O/$RUN/*/projected.tsv; do
    if [ $first = 1 ]; then cat $f > $TMP; first=0; else tail -n +2 $f >> $TMP; fi
  done
  if [ "$O" = "$S/runs" ]; then mv $TMP $S/$RUN.projected.tsv; else mv $TMP $O/$RUN.projected.tsv; fi
) 9> $O/$RUN/.tsv.lock
echo "task done"

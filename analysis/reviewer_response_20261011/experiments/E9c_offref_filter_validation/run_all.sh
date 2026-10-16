#!/bin/bash
# E9c / E9d driver: submit every vcf dir that is ready and not yet processed, the gold-set calibration (once), the
# assembly checks (after the runs) and collect.py (after everything). Re-run it whenever new RIS re-predictions
# (reviewer_cpu_20261013/<run>/vcf_chr1*) have finished; finished runs are skipped.
#
# usage: bash run_all.sh [--no-asm]
# Job scripts and logs: /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E9c_offref_filter_validation
set -euo pipefail
E=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E9c_offref_filter_validation
J=/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E9c_offref_filter_validation
D=/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation
P=/wanglab/jshen/anaconda3/bin/python
mkdir -p $J/logs
if [ -n "$(squeue -u $USER -h -n e9c_runs,e9c_gold,e9c_gold_d46,e9d_asm_HG008T,e9d_asm_COLO829T,e9c_collect -o %i)" ]; then
  echo "E9c jobs are still queued or running (squeue -n e9c_runs,...): run this again when they are done"; exit 0
fi
STAMP=$(date +%Y%m%d_%H%M%S)
$P $E/discover.py > $J/targets_$STAMP.tsv
N=$(wc -l < $J/targets_$STAMP.tsv)
DEPS=""
if [ "$N" -gt 0 ]; then
  A=$(sbatch --parsable -J e9c_runs --array=0-$((N - 1))%4 -o $J/logs/run_%A_%a.out $J/e9c_task.sh $J/targets_$STAMP.tsv)
  echo "e9c runs: $N tasks, array job $A ($J/targets_$STAMP.tsv)"
  DEPS="$DEPS:$A"
else
  echo "e9c runs: nothing new"
fi
if [ ! -s $D/gold/HG008T_ONT/calls.tsv ]; then
  G=$(sbatch --parsable -J e9c_gold -o $J/logs/gold_%j.out $J/gold.sh)
  echo "gold: job $G"
  DEPS="$DEPS:$G"
fi
if [ "${1:-}" != "--no-asm" ]; then
  M=""
  for S in HG008T COLO829T; do
    X=$(sbatch --parsable ${DEPS:+--dependency=afterany$DEPS} -J e9d_asm_$S -o $J/logs/asm_${S}_%j.out $J/asm.sh $S)
    echo "assembly check $S: job $X"
    M="$M:$X"
  done
  C=$(sbatch --parsable --dependency=afterany$M -J e9c_collect -o $J/logs/collect_%j.out $J/collect.sh)
  echo "collect: job $C"
fi

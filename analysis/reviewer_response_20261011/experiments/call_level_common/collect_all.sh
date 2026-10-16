#!/bin/bash
# Rerun all four call-level analyses (E4b, E19, E10, E2.4) on every finished call set: two Slurm chains,
# E4b -> E10 and E19 -> E2.4 (-p general, 4 CPUs each, so at most two jobs run at once). Each job rewrites its
# folder's tables.md and summary.json; logs in tmp/reviewer_cpu_20261013/<ID>/logs/. README numbers are not updated.
set -euo pipefail
J=/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013
for chain in "E4b_fp_decomposition E10_region_strata" "E19_germline_density E2.4_vaf_bins"; do
  dep=""
  for id in $chain; do
    jid=$(sbatch --parsable ${dep:+--dependency=afterany:$dep} -J $id -o $J/$id/logs/%x-%j.out $J/$id/job.sh)
    echo "$id $jid"
    dep=$jid
  done
done

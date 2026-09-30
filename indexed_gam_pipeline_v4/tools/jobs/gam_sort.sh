#!/bin/bash
# One-time sort of one GAM with this package's tools.gam_prep: vg gamsort -i into OUTPUT and OUTPUT.gai, then the
# check the pipeline needs (GAI readable, first records in node order); both are published only after the check.
# With VERIFY=1 it also compares the vg stats -a counts of OUTPUT with the input (reads both GAMs once more).
# Usage: sbatch -J NAME -o LOG gam_sort.sh GAM [OUTPUT]   (OUTPUT default: GAM without .gam + .sorted.gam)
# Environment: PANSOMA_VG (vg; default vg on PATH; scripts/use_vg.sh sets it), GAMSORT_TMP (vg gamsort's temporary
# files; default /scratch/jshen/tmp_gamsort), VERIFY.
# Runs the package it lives in (<package>/tools/jobs/): the checkout or pipeline_code/.
#SBATCH -p general
#SBATCH -c 10
#SBATCH --mem=64G
#SBATCH -t 6-00:00:00
set -euo pipefail
S=${BASH_SOURCE[0]}  # sbatch runs a copy of this file; the submitted path is in the job record (always absolute)
[ -f "$(dirname "$S")/../../__init__.py" ] ||
  S=$(scontrol show job "${SLURM_JOB_ID:?this copy is not in <package>/tools/jobs/ and not a Slurm job}" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p')
PKG=$(cd "$(dirname "$S")/../.." && pwd)
[ -f "$S" ] && [ -f "$PKG/__init__.py" ] || { echo "Error: $S is not in <package>/tools/jobs/"; exit 1; }
echo "package: $PKG"
[ $# -ge 1 ] && [ $# -le 2 ] || { sed -n 5p "$S"; exit 2; }
abs() { case $1 in /*) echo "$1" ;; *) echo "$PWD/$1" ;; esac; }  # the job runs from the package's parent
GAM=$(abs "$1"); OUTPUT=$(abs "${2:-${GAM%.gam}.sorted.gam}")
TMP=${GAMSORT_TMP:-/scratch/jshen/tmp_gamsort}
T=$(( ${SLURM_CPUS_PER_TASK:-4} > 2 ? ${SLURM_CPUS_PER_TASK:-4} - 2 : 1 ))  # gamsort workers; the rest for its I/O
PY=/wanglab/jshen/anaconda3/bin/python
echo "GAM $GAM -> $OUTPUT; tmp $TMP; vg $(command -v "${PANSOMA_VG:-vg}"): $("${PANSOMA_VG:-vg}" version | head -1)"
cd "$(dirname "$PKG")"
M=$(basename "$PKG").tools.gam_prep
/usr/bin/time -v $PY -m $M sort --gam "$GAM" --output "$OUTPUT" --threads "$T" --tmp-dir "$TMP"
if [ "${VERIFY:-0}" = 1 ]; then
  /usr/bin/time -v $PY -m $M check --gam "$OUTPUT" --input "$GAM" --threads "$(( T > 1 ? T / 2 : 1 ))"
fi
echo GAM_SORT_DONE

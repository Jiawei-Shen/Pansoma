#!/bin/bash
# One-time record index <GAM>.gri of one sorted long-read GAM with this package's gam_record_index: build (one
# parallel pass over the GAM, refused with nothing written unless every block, record and GAI run checks out), then
# check (trailer, stamps, structure, GAI map, re-read records and groups). Never replaces an existing <GAM>.gri.
# Usage: sbatch -J NAME -o LOG build_record_index.sh GAM [INDEX]   (INDEX default: GAM + .gai)
# Environment: FORCE=1 indexes a short-read GAM too. Runs the package it lives in (<package>/tools/jobs/).
#SBATCH -p general,gpu
#SBATCH --exclude=tsingtao
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 1-00:00:00
set -euo pipefail
S=${BASH_SOURCE[0]}  # sbatch runs a copy of this file; the submitted path is in the job record (always absolute)
[ -f "$(dirname "$S")/../../__init__.py" ] ||
  S=$(scontrol show job "${SLURM_JOB_ID:?this copy is not in <package>/tools/jobs/ and not a Slurm job}" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p')
PKG=$(cd "$(dirname "$S")/../.." && pwd)
[ -f "$S" ] && [ -f "$PKG/__init__.py" ] || { echo "Error: $S is not in <package>/tools/jobs/"; exit 1; }
echo "package: $PKG"
[ $# -ge 1 ] && [ $# -le 2 ] || { sed -n 5p "$S"; exit 2; }
abs() { case $1 in /*) echo "$1" ;; *) echo "$PWD/$1" ;; esac; }  # the job runs from the package's parent
GAM=$(abs "$1"); INDEX=$(abs "${2:-$GAM.gai}")
PY=/wanglab/jshen/anaconda3/bin/python
cd "$(dirname "$PKG")"
M=$(basename "$PKG").gam_record_index
/usr/bin/time -v $PY -m $M build --gam "$GAM" --index "$INDEX" --processes "${SLURM_CPUS_PER_TASK:-8}" \
  ${FORCE:+--force}
/usr/bin/time -v $PY -m $M check --gam "$GAM" --index "$INDEX"
echo BUILD_RECORD_INDEX_DONE

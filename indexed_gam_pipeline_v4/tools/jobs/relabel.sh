#!/bin/bash
# Labels one merged tensor set with this package's truth-label rules (tensor_postprocessing label), the graph's
# reference-path directory and the GRCh38 FASTA. No backup is kept: label replaces a kind's label files only when all
# of them are written, so a failed run leaves the current labels in place.
# Usage: sbatch -J NAME -o LOG relabel.sh TENSORS SOMATIC_VCF SOMATIC_BED GERMLINE_VCF GERMLINE_BED TRUTH_DIR [SNV_MIN_AF [INDEL_MIN_AF]]
# SNV_MIN_AF (short-read sets: 0.07), INDEL_MIN_AF: lower-AF SNV / INDEL tensors become -1; '' skips one (e.g. '' 0.10).
# REFERENCE_PATH (environment): the ref-path-scan directory of the set's graph; default: the one its current labels
# used (reference_path in TENSORS/SNV/labels.manifest.json).
# Runs the package it lives in (<package>/tools/jobs/): the checkout, or
# /scratch/jshen/data/pansoma_v2_tensors/pipeline_code/ (git archive, commit in pipeline_code/git_head.txt).
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=19G
#SBATCH -t 6:00:00
set -euo pipefail
S=${BASH_SOURCE[0]}  # sbatch runs a copy of this file; the submitted path is in the job record (always absolute)
[ -f "$(dirname "$S")/../../__init__.py" ] ||
  S=$(scontrol show job "${SLURM_JOB_ID:?this copy is not in <package>/tools/jobs/ and not a Slurm job}" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p')
PKG=$(cd "$(dirname "$S")/../.." && pwd)
[ -f "$S" ] && [ -f "$PKG/__init__.py" ] || { echo "Error: $S is not in <package>/tools/jobs/"; exit 1; }
echo "package: $PKG"
[ $# -ge 6 ] && [ $# -le 8 ] || { sed -n 5p "$S"; exit 2; }
abs() { case $1 in /*) echo "$1" ;; *) echo "$PWD/$1" ;; esac; }  # the job runs from the package's parent
T=$(abs "$1"); SV=$(abs "$2"); SB=$(abs "$3"); GV=$(abs "$4"); GB=$(abs "$5"); TD=$(abs "$6"); AF=${7:-}; IAF=${8:-}
PY=/wanglab/jshen/anaconda3/bin/python
grep -qE '"layout": "chromosome-shards(-v1)?"' "$T/SNV/manifest.json" 2>/dev/null || { echo "$T is not merged"; exit 1; }
RP=${REFERENCE_PATH:-}
if [ -z "$RP" ] && [ -f "$T/SNV/labels.manifest.json" ]; then
  RP=$($PY -c 'import json, sys; print(json.load(open(sys.argv[1])).get("reference_path") or "")' "$T/SNV/labels.manifest.json")
fi
[ -n "$RP" ] || { echo "$T has no labels.manifest.json with a reference_path: set REFERENCE_PATH"; exit 1; }
RP=$(abs "$RP"); echo "reference path: $RP"
cd "$(dirname "$PKG")"
/usr/bin/time -v $PY -m "$(basename "$PKG")".tensor_postprocessing label --tensors "$T" \
  --reference-path "$RP" \
  --fasta /scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --somatic-vcf "$SV" --somatic-bed "$SB" --germline-vcf "$GV" --germline-bed "$GB" --truth-dir "$TD" \
  ${AF:+--snv-min-af $AF} ${IAF:+--indel-min-af $IAF}

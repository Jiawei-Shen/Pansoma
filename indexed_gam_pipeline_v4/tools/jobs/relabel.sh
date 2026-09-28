#!/bin/bash
# Labels one merged tensor set with this package's truth-label rules (tensor_postprocessing label), the HPRC v1.1 d9
# reference-path directory and the GRCh38 FASTA. The current label manifests and recall files are copied first to
# <sample dir>/labels_backup_<tensors dir>_<time>_<job>/.
# Usage: sbatch -J NAME -o LOG relabel.sh TENSORS SOMATIC_VCF SOMATIC_BED GERMLINE_VCF GERMLINE_BED TRUTH_DIR [SNV_MIN_AF]
# SNV_MIN_AF: short-read sets only (0.07): lower-AF SNV tensors become -1.
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
[ $# -ge 6 ] && [ $# -le 7 ] || { sed -n 5p "$S"; exit 2; }
abs() { case $1 in /*) echo "$1" ;; *) echo "$PWD/$1" ;; esac; }  # the job runs from the package's parent
T=$(abs "$1"); SV=$(abs "$2"); SB=$(abs "$3"); GV=$(abs "$4"); GB=$(abs "$5"); TD=$(abs "$6"); AF=${7:-}
G=/scratch/jshen/data/pansoma_v2_tensors/graph_index; PY=/wanglab/jshen/anaconda3/bin/python
grep -qE '"layout": "chromosome-shards(-v1)?"' "$T/SNV/manifest.json" 2>/dev/null || { echo "$T is not merged"; exit 1; }
if [ -f "$T/SNV/labels.manifest.json" ]; then
  B=$(dirname "$T")/labels_backup_$(basename "$T")_$(date +%Y%m%dT%H%M%S)_${SLURM_JOB_ID:-local}; mkdir "$B"
  for k in SNV INDEL; do cp -p "$T/$k/labels.manifest.json" "$B/$k.labels.manifest.json"; done
  cp -p "$T"/*.recall.tsv "$T/truth_recall.json" "$B/" 2>/dev/null || true
  echo "previous labels backed up to $B"
fi
cd "$(dirname "$PKG")"
/usr/bin/time -v $PY -m "$(basename "$PKG")".tensor_postprocessing label --tensors "$T" \
  --reference-path $G/hprc-v1.1-mc-grch38.d9.grch38_path \
  --fasta /scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --somatic-vcf "$SV" --somatic-bed "$SB" --germline-vcf "$GV" --germline-bed "$GB" --truth-dir "$TD" ${AF:+--snv-min-af $AF}

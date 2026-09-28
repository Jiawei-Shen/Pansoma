#!/bin/bash
# One-time preparation of one graph with this package's tools.graph_prep: ref-path-scan, ref-path-check, chr-index.
# Writes OUTDIR/<name>.grch38_path/ and OUTDIR/<name>.chr_node_ranges.{tsv,json} (<name> = GFA file name without
# .gfa); stops before the scan when any of them exists. Defaults: HPRC v1.1 d9 (the files under
# /scratch/jshen/data/pansoma_v2_tensors/graph_index).
# Usage: sbatch -J graph_prep -o LOG graph_prep.sh OUTDIR [GFA [GRAPH_INDEX [FASTA [COMPONENTS_DIR]]]]
# COMPONENTS_DIR holds chrN/chrN.component.nodes.raw.txt (vg chunk -C) for chr1-22.
# Runs the package it lives in (<package>/tools/jobs/): the checkout or pipeline_code/.
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=5G
#SBATCH -t 1:00:00
set -euo pipefail
S=${BASH_SOURCE[0]}  # sbatch runs a copy of this file; the submitted path is in the job record (always absolute)
[ -f "$(dirname "$S")/../../__init__.py" ] ||
  S=$(scontrol show job "${SLURM_JOB_ID:?this copy is not in <package>/tools/jobs/ and not a Slurm job}" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p')
PKG=$(cd "$(dirname "$S")/../.." && pwd)
[ -f "$S" ] && [ -f "$PKG/__init__.py" ] || { echo "Error: $S is not in <package>/tools/jobs/"; exit 1; }
echo "package: $PKG"
[ $# -ge 1 ] && [ $# -le 5 ] || { sed -n 6p "$S"; exit 2; }
abs() { case $1 in /*) echo "$1" ;; *) echo "$PWD/$1" ;; esac; }  # the job runs from the package's parent
OUT=$(abs "$1")
GFA=$(abs "${2:-/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gfa}")
DB=$(abs "${3:-/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.graph_index.sqlite}")
FASTA=$(abs "${4:-/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta}")
COMPONENTS=$(abs "${5:-/scratch/jshen/data/AF-Filtered_VG_Indexes/chr_component_vs_GRCh38_summary}")
PY=/wanglab/jshen/anaconda3/bin/python
RP=$OUT/$(basename "$GFA" .gfa).grch38_path; CI=$OUT/$(basename "$GFA" .gfa).chr_node_ranges
for f in "$RP" "$RP.tmp" "$CI.tsv" "$CI.json"; do [ ! -e "$f" ] || { echo "Error: $f exists"; exit 1; }; done
cd "$(dirname "$PKG")"
M=$(basename "$PKG").tools.graph_prep
/usr/bin/time -v $PY -m $M ref-path-scan --gfa "$GFA" --output "$RP"
/usr/bin/time -v $PY -m $M ref-path-check --path "$RP" --graph-index "$DB" --fasta "$FASTA"
/usr/bin/time -v $PY -m $M chr-index --components-dir "$COMPONENTS" --reference-path "$RP" --output "$CI" --graph-index "$DB"
echo GRAPH_PREP_DONE

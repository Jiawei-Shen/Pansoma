#!/bin/bash
# One-time preparation of one graph with this package's tools; a step whose output exists is skipped, so a
# resubmission resumes. With <name> = the GBZ file name without .gbz, it writes:
#   <GBZ dir>/<name>.gfa                         graph_prep gfa (segment names = GBZ node IDs)
#   OUTDIR/gbz_graph_index (+ .build.json)       graph_index_build compile (once per OUTDIR)
#   OUTDIR/<name>.graph_index.sqlite, OUTDIR/graph_index_report.json   graph_index_build build
#   OUTDIR/<name>.grch38_path/                   graph_prep ref-path-scan, then ref-path-check (graph index, FASTA)
#   OUTDIR/<name>.components/                    graph_prep components (vg chunk -C per autosome)
#   OUTDIR/<name>.chr_node_ranges.{tsv,json}     graph_prep chr-index
#   OUTDIR/graph_audit.json                      graph_prep audit (graph index vs GFA)
#   OUTDIR/<name>.chm13_path/                    with CHM13_FASTA only: ref-path-scan --reference-sample CHM13, then
#                                                ref-path-check against CHM13_FASTA (pansoma_net_v2 linear_vcf: CHM13 VCF)
# The GFA and the graph index are built side by side (both read only the GBZ).
# Usage: sbatch -J graph_prep -o LOG graph_prep.sh GBZ OUTDIR [FASTA]
# Environment: PANSOMA_VG (vg for gfa and components; default vg on PATH; scripts/use_vg.sh sets it), GFA (default
# next to the GBZ), GBZ_DEPS (gbwtgraph prefix for the builder; default /scratch/jshen/Github/gbz-tool/dependency),
# REFERENCE_SAMPLE (ref-path-scan --reference-sample; default GRCh38; _gbwt_ref for a GBZ whose reference is generic
# paths, e.g. a GRCh38-only vg autoindex graph), CHM13_FASTA (optional: the T2T-CHM13 FASTA of the graph's CHM13
# reference path, e.g. chm13v2.0.fa; adds the CHM13 reference path), CHM13_EXCLUDE_CONTIGS (CHM13 contigs left off
# that path, space-separated; chrY for the HPRC v1.1 graphs, whose CHM13 chrY is GRCh38's chrY).
# Runs the package it lives in (<package>/tools/jobs/): the checkout or pipeline_code/.
#SBATCH -p general
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 2-00:00:00
set -euo pipefail
S=${BASH_SOURCE[0]}  # sbatch runs a copy of this file; the submitted path is in the job record (always absolute)
[ -f "$(dirname "$S")/../../__init__.py" ] ||
  S=$(scontrol show job "${SLURM_JOB_ID:?this copy is not in <package>/tools/jobs/ and not a Slurm job}" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p')
PKG=$(cd "$(dirname "$S")/../.." && pwd)
[ -f "$S" ] && [ -f "$PKG/__init__.py" ] || { echo "Error: $S is not in <package>/tools/jobs/"; exit 1; }
echo "package: $PKG"
[ $# -ge 2 ] && [ $# -le 3 ] || { sed -n 14p "$S"; exit 2; }
abs() { case $1 in /*) echo "$1" ;; *) echo "$PWD/$1" ;; esac; }  # the job runs from the package's parent
GBZ=$(abs "$1"); OUT=$(abs "$2")
FASTA=$(abs "${3:-/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta}")
NAME=$(basename "$GBZ" .gbz)
CHM13_FASTA=${CHM13_FASTA:+$(abs "$CHM13_FASTA")}
GFA=$(abs "${GFA:-$(dirname "$GBZ")/$NAME.gfa}")
DEPS=${GBZ_DEPS:-/scratch/jshen/Github/gbz-tool/dependency}
T=${SLURM_CPUS_PER_TASK:-4}
PY=/wanglab/jshen/anaconda3/bin/python
DB=$OUT/$NAME.graph_index.sqlite; RP=$OUT/$NAME.grch38_path; CO=$OUT/$NAME.components; CI=$OUT/$NAME.chr_node_ranges
CP=$OUT/$NAME.chm13_path
for f in "$RP.tmp" "$CO.tmp" "$CP.tmp"; do [ ! -e "$f" ] || { echo "Error: $f is left from a failed step; remove it first"; exit 1; }; done
echo "GBZ $GBZ; OUTDIR $OUT; GFA $GFA; FASTA $FASTA; vg $(command -v "${PANSOMA_VG:-vg}"): $("${PANSOMA_VG:-vg}" version | head -1)"
mkdir -p "$OUT"
cd "$(dirname "$PKG")"
M=$(basename "$PKG").tools
TIME="/usr/bin/time -v"
done_or() { if [ -e "$2" ]; then echo "== $1: $2 exists, skipped"; return 1; fi; echo "== $1"; }
if done_or builder "$OUT/gbz_graph_index"; then
  $PY -m $M.graph_index_build compile --deps "$DEPS" --output "$OUT/gbz_graph_index"
fi
pids=()
if done_or gfa "$GFA"; then
  $TIME $PY -m $M.graph_prep gfa --gbz "$GBZ" --output "$GFA" --threads "$T" & pids+=($!)
fi
if done_or graph-index "$DB"; then
  ( $TIME $PY -m $M.graph_index_build build --gbz "$GBZ" --builder "$OUT/gbz_graph_index" --output "$DB" \
      > "$OUT/graph_index_report.json.tmp" && mv "$OUT/graph_index_report.json.tmp" "$OUT/graph_index_report.json" ) &
  pids+=($!)
fi
for p in "${pids[@]}"; do wait "$p"; done
if done_or ref-path-scan "$RP"; then
  $TIME $PY -m $M.graph_prep ref-path-scan --gfa "$GFA" --output "$RP" --reference-sample "${REFERENCE_SAMPLE:-GRCh38}"
fi
echo "== ref-path-check"
$TIME $PY -m $M.graph_prep ref-path-check --path "$RP" --graph-index "$DB" --fasta "$FASTA"
if done_or components "$CO"; then
  $TIME $PY -m $M.graph_prep components --gbz "$GBZ" --reference-path "$RP" --output "$CO" --threads "$T"
fi
if done_or chr-index "$CI.json"; then
  $TIME $PY -m $M.graph_prep chr-index --components-dir "$CO" --reference-path "$RP" --output "$CI" --graph-index "$DB"
fi
if done_or audit "$OUT/graph_audit.json"; then
  $TIME $PY -m $M.graph_prep audit --graph-index "$DB" --gfa "$GFA" --chr-index "$CI.tsv" --output "$OUT/graph_audit.json" \
    --processes "$T"
fi
if [ -n "$CHM13_FASTA" ]; then
  if done_or chm13-path-scan "$CP"; then
    $TIME $PY -m $M.graph_prep ref-path-scan --gfa "$GFA" --output "$CP" --reference-sample CHM13 \
      ${CHM13_EXCLUDE_CONTIGS:+--exclude-contigs $CHM13_EXCLUDE_CONTIGS}
  fi
  echo "== chm13-path-check"
  $TIME $PY -m $M.graph_prep ref-path-check --path "$CP" --graph-index "$DB" --fasta "$CHM13_FASTA"
fi
echo GRAPH_PREP_DONE

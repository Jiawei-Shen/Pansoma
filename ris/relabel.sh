#!/bin/bash
# Relabel one merged set on RIS, as indexed_gam_pipeline_v4/tools/jobs/relabel.sh does on the old cluster:
# tensor_postprocessing label of this checkout, with the truth files of the set's current labels (their old-cluster
# paths, bind-mounted by submit.sh), the graph's GRCh38 path and the GRCh38 FASTA, and a new SNV AF floor.
# label keeps no backup, so the current label files (SNV and INDEL), recall reports and truth tables are copied
# first to <set>/labels_backup_<date the current labels were made>/ (once; an existing backup is kept).
# Usage: relabel.sh SET SNV_MIN_AF ('' for none). A CPU job, 9 slots (label's recall scan reads the filtered-candidate
# streams in 8 processes besides the main one):
#   ris/submit.sh -g 0 -n 9 -m 32 pnv2_relabel_HG008T_PacBio relabel.sh HG008T_PacBio ''
set -euo pipefail
source "$(dirname "$0")/env.sh"
SET=$1 AF=$2
S=$PANSOMA_DATA/pansoma_v2_tensors/$SET T=$PANSOMA_DATA/pansoma_v2_tensors/$SET/tensors
read -r SV SB GV GB CREATED < <(python -c "import json, sys; m = json.load(open(sys.argv[1])); t = m['truth']
print(t['somatic']['vcf'], t['somatic']['bed'], t['germline']['vcf'], t['germline']['bed'], m['created'][:10])" "$T/SNV/labels.manifest.json")
B=$S/labels_backup_$CREATED
if [ ! -e "$B" ]; then
  rm -rf "$B.tmp"; mkdir -p "$B.tmp"
  for k in SNV INDEL; do
    mkdir -p "$B.tmp/$k"; cp -p "$T/$k"/*_labels.npy "$T/$k"/*_labels.ndjson "$T/$k/labels.manifest.json" "$B.tmp/$k/"
  done
  cp -p "$T"/*.recall.tsv "$T/truth_recall.json" "$B.tmp/"; cp -rp "$S/truth" "$B.tmp/truth"
  mv "$B.tmp" "$B"; echo "backup of the $CREATED labels: $B ($(du -sh "$B" | cut -f1))"
else
  echo "backup exists, kept: $B"
fi
for f in "$SV" "$SB" "$GV" "$GB"; do [ -f "$f" ] || { echo "missing truth file $f (submit.sh mounts)"; exit 1; }; done
cd "$REPO"
time python -u -m indexed_gam_pipeline_v4.tensor_postprocessing label --tensors "$T" \
  --reference-path "$PANSOMA_DATA/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path" \
  --fasta "$PANSOMA_DATA/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" \
  --somatic-vcf "$SV" --somatic-bed "$SB" --germline-vcf "$GV" --germline-bed "$GB" --truth-dir "$S/truth" \
  ${AF:+--snv-min-af "$AF"}
python -c "import json, sys; m = json.load(open(sys.argv[1])); print('new labels', m['created'], m['rules_sha256'][:12],
'snv_min_af', m['snv_min_af'], m['totals'])" "$T/SNV/labels.manifest.json"

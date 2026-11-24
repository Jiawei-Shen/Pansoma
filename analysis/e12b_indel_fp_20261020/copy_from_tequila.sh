#!/bin/bash
# Copy the multi-cancer INDEL data from tequila to RIS (run ON TEQUILA, inside tmux/screen; written 2026-10-20).
# INDEL only: no SNV/, no filtered_candidates.ndjson, no labels_backup_*. The R13-free labels are derived on RIS
# from the copied labels (dropping R13 only turns repeat_indel_without_truth -1 into 0; truth_labels.classify, non()),
# so do NOT relabel on tequila before copying.
#   bash copy_from_tequila.sh phase1     # truth, E12b checkpoints, Illumina + PacBio of the four lines (~1.55 TB)
#   bash copy_from_tequila.sh ont        # ONT of the four lines (~2.53 TB): only after space is freed on storage2
set -euo pipefail
RIS=${RIS:-jiawei.s@compute1-client-1.ris.wustl.edu}   # the RIS host you copied through this morning
DATA=/storage2/fs1/epigenome/Active/jiawei.s/Dataset
RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs
SOCK=$HOME/.ssh/ris_copy.sock
LINES="HCC1395T HCC1954T H1437T H2009T"

# one authenticated ssh connection shared by every rsync (one password / Duo prompt)
[ -S "$SOCK" ] || ssh -fNM -o ControlPersist=24h -S "$SOCK" "$RIS"
SSH="ssh -S $SOCK"

tensor_set() {  # SET: tensors/INDEL + the recall reports, same layout as the RIS sets
  local SET=$1
  rsync -a --partial --info=stats1 -e "$SSH" \
    --include='/INDEL/' --include='/INDEL/manifest.json' --include='/INDEL/labels.manifest.json' \
    --include='/INDEL/*_data.npy' --include='/INDEL/*_labels.npy' --include='/INDEL/*_labels.ndjson' \
    --include='/INDEL/*_variant_summary.ndjson' \
    --include='/somatic.recall.tsv' --include='/germline.recall.tsv' --include='/truth_recall.json' \
    --exclude='*' \
    --rsync-path="mkdir -p $DATA/pansoma_v2_tensors/$SET/tensors && rsync" \
    "/scratch/jshen/data/pansoma_v2_tensors/$SET/tensors/" "$RIS:$DATA/pansoma_v2_tensors/$SET/tensors/" \
    > "copy_$SET.log" 2>&1 && echo "done $SET" || echo "FAILED $SET (see copy_$SET.log; rerun resumes)"
}
export -f tensor_set; export RIS DATA SSH

case ${1:-} in
  phase1)
    # 1. truth VCFs (+ .tbi) and BEDs of the four lines (label manifests' paths; mounted at /scratch/jshen/data/CASTLE)
    for L in HCC1395 HCC1954 H1437 H2009; do
      rsync -a -e "$SSH" --include='*.vcf.gz' --include='*.vcf.gz.tbi' --include='*.bed' --exclude='*' \
        --rsync-path="mkdir -p $DATA/CASTLE/$L/truth && rsync" \
        "/scratch/jshen/data/CASTLE/$L/truth/" "$RIS:$DATA/CASTLE/$L/truth/"
    done
    # 2. the six E12b INDEL checkpoints (next to the small files copied this morning)
    for R in multicancer4 multicancer6; do for P in Illumina PacBio ONT; do
      rsync -a -e "$SSH" --include='best.pth' --include='last.pth' --exclude='*' \
        "/scratch/jshen/data/pansoma_net_v2_runs/${R}_${P}_INDEL/" "$RIS:$RUNS/from_tequila_20261020/runs/${R}_${P}_INDEL/"
    done; done
    # 3. Illumina and PacBio INDEL sets, 4 at a time
    for L in $LINES; do echo "${L}_Illumina_d46"; done | xargs -P 4 -I{} bash -c 'tensor_set {}'
    for L in $LINES; do echo "${L}_PacBio_d46"; done | xargs -P 4 -I{} bash -c 'tensor_set {}'
    ;;
  ont)
    for L in $LINES; do echo "${L}_ONT_d46"; done | xargs -P 4 -I{} bash -c 'tensor_set {}'
    ;;
  *) sed -n 2,7p "$0"; exit 2;;
esac

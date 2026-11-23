#!/bin/bash
# (wangcluster, in tmux) Copy to storage2 what the nine six-tumor nositef models RIS cannot train yet read; submit.sh
# on RIS then submits every model whose sets have arrived. Written on RIS without access to /scratch: a source
# directory that does not exist under these names is reported and skipped (rename it here), nothing else is touched.
#   snv        SNV v2.1 (3 models): SNV tensors of {HCC1954T,H1437T,H2009T}_{Illumina,PacBio}_d46_snvaf05 and
#              {HCC1954T,H1437T,H2009T}_ONT_d46 (H1437T_ONT_d46 is partly on storage2 from 2026-11-12; rsync resumes),
#              and the eight SNV site-feature files that are links to /scratch on storage2 (copied as files)
#   indel_d46  INDEL v2.1 (3 models): INDEL tensors of the four CASTLE lines, <line>_<P>_d46
#   indel_d9   INDEL v1.1 (3 models): INDEL tensors of the four CASTLE lines, <line>_<P>_d9
# The INDEL site-feature files of the copied sets are built on RIS (submit.sh submits ris/site_features.sh for a set
# without one). Not copied: labels_backup_*, filtered_candidates.ndjson.
# usage: [RIS=user@host] [JOBS=4] bash analysis/m6_nositef_graphs_20261118/to_ris.sh [snv] [indel_d46] [indel_d9]
set -euo pipefail
RIS=${RIS:-jiawei.s@compute1-client-1.ris.wustl.edu}
A=/storage2/fs1/epigenome/Active/jiawei.s
D=/scratch/jshen/data
LINES="HCC1954T H1437T H2009T"
LOGS=$HOME/ris_rsync_20261118; mkdir -p "$LOGS"
SOCK=$HOME/.ssh/ris_mux_$$
ssh -fNM -S "$SOCK" "$RIS"                     # one login (password / Duo) for every transfer below
trap 'ssh -S "$SOCK" -O exit "$RIS" 2>/dev/null || true' EXIT
E="ssh -S $SOCK"
rv=$($E "$RIS" rsync --version)
if [[ $rv == *zstd* ]]; then Z=(--compress --compress-choice=zstd); else Z=(-z); fi

todo=()  # SET/tensors/KIND
for step in ${*:-snv indel_d46 indel_d9}; do
  case $step in
    snv)
      for l in $LINES; do todo+=("${l}_Illumina_d46_snvaf05/tensors/SNV" "${l}_PacBio_d46_snvaf05/tensors/SNV" "${l}_ONT_d46/tensors/SNV"); done
      f=(); for s in HG008T_Illumina_d46 HG008T_PacBio_d46 HG008T_ONT_d46 COLO829T_ONT_d46 HCC1395T_ONT_d46 HCC1954T_ONT_d46 H1437T_ONT_d46 H2009T_ONT_d46; do
        f+=("$D/snv_sitefeat_20261103/site_features/$s.SNV.npz"); done
      rsync -aL -e "$E" "${f[@]}" "$RIS:$A/pansoma_net_v2_runs/site_features/";;
    indel_d46|indel_d9)
      for l in HCC1395T $LINES; do for p in Illumina PacBio ONT; do todo+=("${l}_${p}_${step#indel_}/tensors/INDEL"); done; done;;
    *) echo "steps: snv indel_d46 indel_d9"; exit 1;;
  esac
done
fail=0
for t in "${todo[@]}"; do
  if [ ! -s "$D/pansoma_v2_tensors/$t/labels.manifest.json" ]; then echo "NOT FOUND, skipped: $D/pansoma_v2_tensors/$t"; fail=1; continue; fi
  while [ "$(jobs -rp | wc -l)" -ge "${JOBS:-4}" ]; do wait -n || fail=1; done
  ( rsync -aR "${Z[@]}" --exclude 'labels_backup_*' --exclude filtered_candidates.ndjson -e "$E" \
      "$D/./pansoma_v2_tensors/$t" "$RIS:$A/Dataset/" > "$LOGS/${t//\//_}.log" 2>&1 \
      && echo "copied $t" || { echo "FAILED $t (see $LOGS; rerun resumes)"; exit 1; } ) &
done
while [ -n "$(jobs -rp)" ]; do wait -n || fail=1; done
[ $fail = 0 ] && echo "all copied; next, on RIS: git pull (if needed) and bash analysis/m6_nositef_graphs_20261118/submit.sh" \
  || echo "some sets were not copied (see above); submit.sh on RIS skips the models that miss a set"

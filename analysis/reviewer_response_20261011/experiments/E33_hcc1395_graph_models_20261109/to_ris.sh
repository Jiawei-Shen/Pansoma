#!/bin/bash
# Step 1 of 3 (wangcluster, in tmux): copy what the five RIS trainings read to storage2 -- the SNV tensors of the five
# sets (4.4 TB apparent, 613 GB on /scratch's disk; rsync with zstd when RIS's rsync has it), the HCC1395T truth (the
# label manifests' /scratch/jshen/data/CASTLE/HCC1395/truth, bind-mounted there by ris/submit.sh) and the SNV site
# features of the two Illumina sets. Ends with a file-count and byte check of every set.
# usage: [RIS=user@host] bash analysis/reviewer_response_20261011/experiments/E33_hcc1395_graph_models_20261109/to_ris.sh
set -euo pipefail
RIS=${RIS:-jiawei.s@compute1-client-1.ris.wustl.edu}
A=/storage2/fs1/epigenome/Active/jiawei.s
D=/scratch/jshen/data
SETS=(HCC1395T_ONT_linear HCC1395T_Illumina_d9_snvaf05 HCC1395T_Illumina_d46_snvaf05 HCC1395T_ONT_d46 HCC1395T_ONT_d9)
LOGS=$HOME/ris_rsync_20261109; mkdir -p "$LOGS"
SOCK=$HOME/.ssh/ris_mux_$$
ssh -fNM -S "$SOCK" "$RIS"                     # one login (password / Duo) for every transfer below
trap 'ssh -S "$SOCK" -O exit "$RIS" 2>/dev/null || true' EXIT
E="ssh -S $SOCK"
rv=$($E "$RIS" rsync --version)
if [[ $rv == *zstd* ]]; then Z=(--compress --compress-choice=zstd); else Z=(-z); fi
echo "RIS rsync: $(head -1 <<< "$rv"); compression ${Z[*]}"

rsync -aR -e "$E" $D/./CASTLE/HCC1395/truth "$RIS:$A/Dataset/"
rsync -a -e "$E" $D/snv_af05_20261103/site_features/HCC1395T_Illumina_{d9,d46}_snvaf05.SNV.npz "$RIS:$A/pansoma_net_v2_runs/site_features/"
pids=()
for s in "${SETS[@]}"; do
  rsync -aR "${Z[@]}" -e "$E" $D/./pansoma_v2_tensors/$s/tensors/SNV "$RIS:$A/Dataset/" > "$LOGS/$s.log" 2>&1 &
  pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
[ $fail = 0 ] || { echo "an rsync failed: see $LOGS/*.log (rerun this script; rsync resumes)"; exit 1; }

count='find "$0" -type f -printf "%s\n" | awk "{s += \$1} END {print NR, s}"'
ok=1
for s in "${SETS[@]}"; do
  here=$(bash -c "$count" "$D/pansoma_v2_tensors/$s/tensors/SNV")
  there=$($E "$RIS" bash -c "'$count'" "$A/Dataset/pansoma_v2_tensors/$s/tensors/SNV")
  [ "$here" = "$there" ] && st=OK || { st=DIFFERENT; ok=0; }
  echo "$s files bytes: here $here, RIS $there $st"
done
[ $ok = 1 ] && echo "all five sets copied; next: on RIS, git pull and bash .../E33_hcc1395_graph_models_20261109/ris_submit.sh"

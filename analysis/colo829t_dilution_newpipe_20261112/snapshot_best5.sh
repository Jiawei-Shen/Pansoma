#!/bin/bash
# Keep the best checkpoint of a run's first N epochs (default 5) while the run trains on: waits for the "epoch N/"
# line of <run>/train.log, then copies best.pth (the best epoch so far by the run's own selection, the lines marked
# "| best") with args.json and the statistics into $PANSOMA_RUNS/<RUN>_e<N>/. Stops with an error when a later epoch
# has already replaced best.pth (best epoch > N), or when best.pth changes during the copy. The user (2026-11-12):
# test the best epoch once 5 epochs are done. predict_set.sh / vcf.sh --set then run on <RUN>_e<N>.
# Usage: snapshot_best5.sh RUN [N]. A CPU job (it polls until the run reaches epoch N):
#   ris/submit.sh -g 0 -n 1 -m 4 pnv2_snap5_RUN ../analysis/colo829t_dilution_newpipe_20261112/snapshot_best5.sh RUN
set -euo pipefail
RUN=$1 N=${2:-5}
R=${PANSOMA_RUNS:-/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs}
S=$R/$RUN D=$R/${RUN}_e$N
[ -s "$D/best.pth" ] && { echo "exists: $D/best.pth"; cat "$D/snapshot.txt"; exit 0; }
until grep -q "^epoch $N/" "$S/train.log" 2>/dev/null; do sleep 30; done
sleep 60  # the checkpoint of epoch N is written by now; epoch N + 1 takes minutes
best=$(grep -E "^epoch [0-9]+/.*\| best$" "$S/train.log" | tail -1 | sed -E 's/^epoch ([0-9]+).*/\1/')
[ "$best" -le "$N" ] || { echo "best.pth is already epoch $best (> $N): the best of epochs 1-$N is gone"; exit 1; }
m0=$(stat -c %Y "$S/best.pth")
mkdir -p "$D"
cp -p "$S/best.pth" "$D/best.pth.tmp"
[ "$(stat -c %Y "$S/best.pth")" = "$m0" ] || { rm -f "$D/best.pth.tmp"; echo "best.pth changed during the copy"; exit 1; }
cp -p "$S/args.json" "$S/stats.json" "$S/dataset_stats.json" "$D/"
grep -E "^epoch ([1-9]|[1-9][0-9])/" "$S/train.log" | awk -F'[ /]' -v n="$N" '$2 <= n' > "$D/train_epochs1-$N.log"
echo "$RUN: best of epochs 1-$N = epoch $best; best.pth of $(date -d @"$m0" '+%F %T') copied $(date '+%F %T')" | tee "$D/snapshot.txt"
mv "$D/best.pth.tmp" "$D/best.pth"

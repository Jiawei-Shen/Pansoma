#!/bin/bash
# Move the PASS threshold recorded in $PANSOMA_RUNS/<RUN>/best.pth by OFFSET (set_call_threshold.py; the user,
# 2026-11-12: Illumina and fiberseq call at the validation threshold - 0.10), then print the threshold graph_vcf now
# takes by default. ORIGINAL_RUN holds the unchanged, byte-identical best.pth, which stays. A line is added to
# <RUN>/snapshot.txt. vcf.sh without a THRESHOLD then calls at the new threshold.
# Usage: set_call_threshold.sh RUN OFFSET ORIGINAL_RUN. A CPU job:
#   ris/submit.sh -g 0 -n 2 -m 16 pnv2_setthr_RUN ../analysis/colo829t_dilution_newpipe_20261112/set_call_threshold.sh \
#       COLO829_Illumina_SNV_joint_newpipe_e5 -0.10 COLO829_Illumina_SNV_joint_newpipe
set -euo pipefail
source "$(dirname "$0")/../../ris/env.sh"
RUN=$1 OFFSET=$2 ORIGINAL=$3
D=$PANSOMA_RUNS/$RUN
python -u "$(dirname "$0")/set_call_threshold.py" "$D/best.pth" "$OFFSET" --original "$PANSOMA_RUNS/$ORIGINAL/best.pth" \
    --note "user 2026-11-12: PASS at the validation call threshold $OFFSET" | tee -a "$D/snapshot.txt"
cd "$REPO/machine_learning"
python - "$D/best.pth" <<'EOF'
import sys
from pansoma_net_v2 import graph_vcf
c = graph_vcf.read_checkpoint(sys.argv[1])
print("graph_vcf default:", c["call_rule"], c["call_recall"],
      graph_vcf.postpon_recall_threshold(c, sys.argv[1], c["call_recall"]))
EOF

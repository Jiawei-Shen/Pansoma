#!/bin/bash
# The CPU unit tests of pansoma_net_v2 in the container (GPU and real-data tests skip themselves), with the full
# failure output. Usage: unit_tests.sh [unittest discover options]. A CPU job:
#   ris/submit.sh -g 0 -n 4 -m 16 pnv2_unit_tests unit_tests.sh
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$REPO/machine_learning"
PANSOMA_TEST_TENSORS=${PANSOMA_TEST_TENSORS:-/nonexistent} python -m unittest discover -s pansoma_net_v2/tests -t . "$@" 2>&1 \
    | grep -vE "^\s*$" | tail -n 80

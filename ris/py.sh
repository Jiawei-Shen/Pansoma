#!/bin/bash
# Run a Python script of ris/tools (or any path) in the container with the ris environment, from machine_learning/
# so that pansoma_net_v2 imports. Usage: py.sh SCRIPT [args]   e.g. ris/submit.sh -g 0 -n 4 -m 16 NAME py.sh tools/x.py
set -euo pipefail
source "$(dirname "$0")/env.sh"
S=$1; shift
case $S in /*) ;; *) S=$REPO/ris/$S ;; esac
cd "$REPO/machine_learning"
PYTHONPATH="$REPO/machine_learning:$REPO/scripts:$REPO/ris/tools${PYTHONPATH:+:$PYTHONPATH}" python -u "$S" "$@"

#!/usr/bin/env bash
set -euo pipefail

unset PYTHONHOME
export PYTHONNOUSERSITE=1

PANSOMA_PYTHON="${PANSOMA_PYTHON:-/opt/venv/bin/python}"
PANSOMA_WORKFLOW="/opt/pansoma/scripts/pansoma_workflow.py"

if [[ $# -eq 0 ]]; then
    exec "${PANSOMA_PYTHON}" "${PANSOMA_WORKFLOW}" --help
fi

case "$1" in
    doctor|build-node-filters|train|infer|-h|--help)
        exec "${PANSOMA_PYTHON}" "${PANSOMA_WORKFLOW}" "$@"
        ;;
    *)
        exec "$@"
        ;;
esac

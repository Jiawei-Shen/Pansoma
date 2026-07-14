#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

PYTHON_BIN="${PANSOMA_PYTHON:-python}"

echo "Building fast_writer with: $(${PYTHON_BIN} -c 'import sys; print(sys.executable)')"

if ! "${PYTHON_BIN}" -c 'import setuptools, pybind11' >/dev/null 2>&1; then
    cat >&2 <<'EOF'
ERROR: setuptools and pybind11 must be installed in the active Python environment.

Create and activate the supported environment first:
  conda env create -f environment.yml
  conda activate pangenome-ml-data-generation

If the environment already exists:
  conda env update -n pangenome-ml-data-generation -f environment.yml
EOF
    exit 1
fi

"${PYTHON_BIN}" setup.py build_ext --inplace

PYTHONPATH="${PWD}/src:${PWD}${PYTHONPATH:+:${PYTHONPATH}}" \
    "${PYTHON_BIN}" -c 'import fast_writer; print("fast_writer ready:", fast_writer.__file__)'

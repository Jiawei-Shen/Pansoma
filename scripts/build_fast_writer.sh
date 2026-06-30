#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
python setup.py build_ext --inplace


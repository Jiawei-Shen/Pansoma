#!/bin/bash
# E22 check 1, step 2: the per-allele join and summary (after subset_resources.sh). 1 core, a few GB.
set -euo pipefail
HERE=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E22_truth_and_version_checks
export PYTHONDONTWRITEBYTECODE=1
cd $HERE
/usr/bin/time -f "truth_germline_check: %es, max RSS %MKB" /wanglab/jshen/anaconda3/bin/python truth_germline_check.py
[ -f make_tables.py ] && /wanglab/jshen/anaconda3/bin/python make_tables.py || true

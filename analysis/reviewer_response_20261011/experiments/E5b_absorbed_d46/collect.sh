#!/bin/bash
# E5b collect: finish the tables from the job outputs (login node, a few minutes).
#   needs: spell_d9.tsv.gz, spell_d46.tsv.gz (Slurm s4 jobs); for linear_recovery.md the CPU re-predictions' rtg outputs
set -euo pipefail
cd /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E5b_absorbed_d46
PY=/wanglab/jshen/anaconda3/bin/python
D=/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46
for g in d9 d46; do [ -s $D/spell_$g.tsv.gz ] || { echo "missing $D/spell_$g.tsv.gz (s4 job)"; exit 1; }; done
$PY t2_compare_d9.py | tee t2_compare_d9.txt
$PY s5_tables.py > /dev/null
$PY s6_linear_recovery.py > /dev/null
echo "written: tables.md summary.json lost_alleles.tsv.gz linear_recovery.md t2_compare_d9.txt"

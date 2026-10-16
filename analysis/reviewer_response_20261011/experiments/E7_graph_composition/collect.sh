#!/bin/bash
# E7: submit the final step (s3 curves, s4 calibration, s6 relatedness) once s1 / s2 / s5 are complete, then build the
# markdown tables (s7_tables.py). Run from anywhere; job logs in /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E7_graph_composition.
set -euo pipefail
W=/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition
n1=$(ls $W/hprc_alleles/chr{1..22}.tsv.gz 2>/dev/null | wc -l)
n2=$(ls $W/sample_alleles/{HG008,COLO829,HCC1395}.{somatic,germline}.tsv.gz 2>/dev/null | wc -l)
n5=$(ls $W/d9_alleles/chr{1..22}.tsv.gz 2>/dev/null | wc -l)
echo "s1 $n1/22  s2 $n2/6  s5 $n5/22"
if [ ! -s /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E7_graph_composition/calibration.tsv ]; then
  [ $n1 = 22 ] && [ $n2 = 6 ] && [ $n5 = 22 ] || { echo 'inputs incomplete'; exit 1; }
  sbatch --wait /scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E7_graph_composition/s3_s4_s6.sh
fi
/wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E7_graph_composition/s7_tables.py

#!/bin/bash
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 02:00:00
#SBATCH -J q9_colo_ill_d9
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q9/logs/%x-%j.out
set -euo pipefail
P=/wanglab/jshen/anaconda3/bin/python
Q=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q9
cd $Q
echo "host $(hostname) start $(date)"
/usr/bin/time -f "e3b_e8 %es %MKB" $P $Q/e3b_e8_q9.py
/usr/bin/time -f "e16 %es %MKB" $P $Q/e16_pairs_q9.py
/usr/bin/time -f "e19 ep8 %es %MKB" $P $Q/e19_colo_ill_d9_q9.py ep8
/usr/bin/time -f "e19 ep6 %es %MKB" $P $Q/e19_colo_ill_d9_q9.py ep6
echo "end $(date)"

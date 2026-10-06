#!/bin/bash
#SBATCH -J s7_parse
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 1:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s7_parse-%j.out
set -euo pipefail
cd /scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001
/usr/bin/time -f '%e s %M KB' /wanglab/jshen/anaconda3/bin/python s7_assembly_paths.py parse
/wanglab/jshen/anaconda3/bin/python s7_assembly_paths.py crosscheck > /scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/asm_crosscheck.log
cat /scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/asm_crosscheck.log

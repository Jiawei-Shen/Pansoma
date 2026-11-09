#!/bin/bash
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 6:00:00
set -euo pipefail
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/tmp/af005_candidates_20261103/count.py $1

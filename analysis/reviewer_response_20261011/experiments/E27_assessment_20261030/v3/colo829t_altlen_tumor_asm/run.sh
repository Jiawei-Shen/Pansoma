#!/bin/bash
#SBATCH -J altlen_asm22
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=28G
#SBATCH -t 2:00:00
#SBATCH -w miao
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/colo829t_altlen_d46_20261103/logs/%x_%j.out
set -euo pipefail
W=/scratch/jshen/Github/Pansoma/tmp/colo829t_altlen_d46_20261103
/wanglab/jshen/anaconda3/bin/python $W/altlen_tumor_asm.py $W/altlen_tumor_asm_22.tsv $W/ids_22.txt

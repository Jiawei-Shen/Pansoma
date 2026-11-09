#!/bin/bash
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 6:00:00
/usr/bin/time -f "%es %MKB" /wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/tmp/snv_af05_20261103/check_tensors.py $1 chr22 chr1

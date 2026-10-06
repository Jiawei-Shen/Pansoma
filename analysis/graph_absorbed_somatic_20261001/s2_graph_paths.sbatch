#!/bin/bash
#SBATCH --job-name=s2_graph_paths
#SBATCH --partition=general
#SBATCH --cpus-per-task=32
#SBATCH --mem=3G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s2_graph_paths.%j.log
cd /scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001
/wanglab/jshen/anaconda3/bin/python s2_graph_paths.py 32

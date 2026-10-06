#!/bin/bash
# Step 4 (analysis/graph_absorbed_somatic_20261001/s4_hprc_walks.py): every HPRC haplotype's local walk through each
# locus of loci.tsv from one parallel pass over the d9 GFA (24 processes over 96 byte ranges).
# Usage: sbatch tmp/graph_absorbed_somatic_20261001/s4_hprc_walks.sh
#SBATCH -p general
#SBATCH -J s4_hprc_walks
#SBATCH -c 24
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s4_hprc_walks-%j.out
set -euo pipefail
cd /scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001
/wanglab/jshen/anaconda3/bin/python -u s4_hprc_walks.py

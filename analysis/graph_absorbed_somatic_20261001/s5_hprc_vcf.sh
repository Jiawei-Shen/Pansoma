#!/bin/bash
#SBATCH -J s5_hprc_vcf
#SBATCH -p general
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 2:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s5_hprc_vcf.%j.log
# per-locus tabix queries of the d9 and full-graph deconstruct VCFs (s5_hprc_vcf.py); needs s5b + s5_index_full_vcf.sh first
/wanglab/jshen/anaconda3/bin/python /scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/s5_hprc_vcf.py

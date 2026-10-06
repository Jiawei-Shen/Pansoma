#!/bin/bash
#SBATCH -J s5_tbi
#SBATCH -p general
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH -t 2:00:00
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/graph_absorbed_somatic_20261001/s5_index_full_vcf.%j.log
# tabix index of the full-graph deconstruct VCF (the original has none); index written next to a symlink in $D
cd /scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/full_vcf_index
/opt/apps/bcftools/1.18/bcftools index --tbi --threads 2 -o hprc-v1.1-mc-grch38.raw.vcf.gz.tbi hprc-v1.1-mc-grch38.raw.vcf.gz
ls -la

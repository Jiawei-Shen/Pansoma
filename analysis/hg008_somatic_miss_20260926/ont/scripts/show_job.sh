#!/bin/bash
#SBATCH -J show_ont
#SBATCH -p general
#SBATCH -c 12
#SBATCH --mem=9G
#SBATCH -t 4:00:00
#SBATCH -o show_%j.out
cd /scratch/jshen/Github/Pansoma/tmp/somatic_miss_analysis_ont
cat show_loci.txt | xargs -P 12 -L 1 bash -c '/wanglab/jshen/anaconda3/bin/python show_locus.py $0 $1 6 > show/$0_$1.txt 2>&1'
for l in $(awk '{print $1"_"$2}' show_loci.txt); do cat show/$l.txt; echo; done > show_examples.txt

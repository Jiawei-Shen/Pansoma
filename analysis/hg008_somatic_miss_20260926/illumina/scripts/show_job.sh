#!/bin/bash
#SBATCH -J show_ilmn
#SBATCH -p general
#SBATCH -c 12
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -o show_%j.out
cd /scratch/jshen/Github/Pansoma/tmp/somatic_miss_analysis_illumina
cat show_loci.txt refgap_loci.txt | xargs -P 12 -L 1 bash -c '/wanglab/jshen/anaconda3/bin/python show_locus.py $0 $1 40 > show/$0_$1.txt 2>&1'
mkdir -p ont_show_compare  # the same loci with the ONT scripts (run from their folder, written here)
cat refgap_loci.txt | xargs -P 3 -L 1 bash -c 'cd ../somatic_miss_analysis_ont && /wanglab/jshen/anaconda3/bin/python show_locus.py $0 $1 40 > ../somatic_miss_analysis_illumina/ont_show_compare/$0_$1.txt 2>&1'
for l in $(awk '{print $1"_"$2}' show_loci.txt); do cat show/$l.txt; echo; done > show_examples.txt

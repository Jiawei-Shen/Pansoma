#!/bin/bash
# E8: index BAMs that have no .bai without writing next to the raw data: a symlink in the E8 bams/ directory and
# samtools index on the link (htslib looks for <link>.bai). usage: sbatch index_bams.sh NAME=/path/to.bam [...]
#SBATCH -p general
#SBATCH -c 8
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
set -euo pipefail
module load samtools/1.18
B=/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/bams
mkdir -p $B
echo "host $(hostname) job ${SLURM_JOB_ID:-}; $(samtools --version | head -1)"
for pair in "$@"; do
  name=${pair%%=*}; bam=${pair#*=}
  [ -s "$bam" ] || { echo "missing $bam"; exit 1; }
  ln -sfn "$bam" $B/$name.bam
  if [ -s $B/$name.bam.bai ]; then echo "$name: index exists"; continue; fi
  /usr/bin/time -f "$name index %es %MKB" samtools index -@ 8 $B/$name.bam
  samtools idxstats $B/$name.bam | awk -v n=$name '$1=="chr1"{print n, "chr1 mapped", $3}'
done

#!/bin/bash
# E8: bwa mem 0.7.17 of paired FASTQs to GRCh38 (GCA_000001405.15 no-alt analysis set, the FASTA of the labels),
# keeping only alignments on chr1 (the test chromosome; disk is 97% full), coordinate-sorted and indexed.
# usage: sbatch bwa_chr1.sh OUT_PREFIX READ_GROUP R1 R2      (OUT_PREFIX under the E8 bams/ directory)
#SBATCH -p general
#SBATCH -c 32
#SBATCH --mem=48G
#SBATCH -t 3-00:00:00
set -euo pipefail
module load bwa/0.7.17 samtools/1.18
OUT=$1; RG=$2; R1=$3; R2=$4
REF=/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
TMP=$OUT.sorttmp; mkdir -p $TMP
echo "host $(hostname) job ${SLURM_JOB_ID:-}; $(date); bwa $(bwa 2>&1 | grep Version); $(samtools --version | head -1)"
echo "bwa mem -t 28 -K 100000000 -R '$RG' $REF $R1 $R2 | samtools view -u -e 'rname==\"chr1\"' | samtools sort"
/usr/bin/time -f "bwa %es %MKB" bwa mem -t 28 -K 100000000 -R "$RG" $REF $R1 $R2 2> $OUT.bwa.log \
  | samtools view -u -e 'rname=="chr1"' - \
  | samtools sort -@ 4 -m 2G -T $TMP/part -o $OUT.chr1.bam -
samtools index -@ 4 $OUT.chr1.bam
samtools flagstat -@ 4 $OUT.chr1.bam > $OUT.chr1.flagstat.txt
tail -3 $OUT.bwa.log; cat $OUT.chr1.flagstat.txt; date

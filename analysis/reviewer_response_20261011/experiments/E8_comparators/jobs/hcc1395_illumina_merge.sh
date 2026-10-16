#!/bin/bash
# E8: HCC1395 Illumina chr1 BAM = the reads of the Pansoma tensor set (SEQC2 Fudan WGS_FD_T_1 + WGS_FD_T_2):
# WGS_FD_T_1 from bwa_chr1.sh (bwa 0.7.17, GCA no-alt) + chr1 of the existing SRR7890824 (= WGS_FD_T_2) BAM
# (bwa 0.7.17 to /scratch/sdong/refGenome/GRCh38/GRCh38.fasta: same 195 contigs and lengths, same chr1 sequence),
# which has no read group: one is added.
#SBATCH -p general
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 12:00:00
set -euo pipefail
module load samtools/1.18
B=/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/bams
echo "host $(hostname) job ${SLURM_JOB_ID:-}; $(date)"
samtools view -u -@ 4 /scratch/jshen/data/SEQC2/SRR7890824_Tumor/SRR7890824_Tumor.sorted.bam chr1 \
  | samtools addreplacerg -@ 4 -r '@RG\tID:WGS_FD_T_2\tSM:HCC1395T\tPL:ILLUMINA\tLB:WGS_FD_T_2' -o $B/HCC1395T_Illumina_FD_T_2.chr1.bam -
samtools merge -@ 8 -f -o $B/HCC1395T_Illumina.chr1.bam $B/HCC1395T_Illumina_FD_T_1.chr1.bam $B/HCC1395T_Illumina_FD_T_2.chr1.bam
samtools index -@ 8 $B/HCC1395T_Illumina.chr1.bam
samtools flagstat -@ 8 $B/HCC1395T_Illumina.chr1.bam > $B/HCC1395T_Illumina.chr1.flagstat.txt
cat $B/HCC1395T_Illumina.chr1.flagstat.txt
rm /scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/bams/HCC1395T_Illumina_FD_T_2.chr1.bam
date

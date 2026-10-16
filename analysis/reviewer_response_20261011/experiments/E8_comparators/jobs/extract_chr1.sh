#!/bin/bash
# E8: chr1 BAM of a coordinate-sorted BAM that has no index, without reading the whole file: the records are streamed
# until the first record off chr1 (chr1 is the first @SQ of these GRCh38 BAMs), written to <E8 bams>/<NAME>.bam and
# indexed. <NAME>.bam may be a symlink to the raw BAM left by index_bams.sh: mv replaces the link, never its target.
# usage: sbatch extract_chr1.sh NAME=/path/to.bam [...]
#SBATCH -p general
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 1-00:00:00
set -euo pipefail
module load samtools/1.18
B=/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/bams
mkdir -p $B
echo "host $(hostname) job ${SLURM_JOB_ID:-}; $(date); $(samtools --version | head -1)"
for pair in "$@"; do
  name=${pair%%=*}; bam=${pair#*=}
  [ -s "$bam" ] || { echo "missing $bam"; exit 1; }
  if [ -f $B/$name.bam ] && [ ! -L $B/$name.bam ] && [ -s $B/$name.bam.bai ]; then echo "$name: chr1 BAM exists"; continue; fi
  first=$(samtools view -H "$bam" | awk '$1=="@SQ"{print $2; exit}')
  [ "$first" = "SN:chr1" ] || { echo "$name: first @SQ is $first, not chr1"; exit 1; }
  samtools view -H "$bam" | grep -q "SO:coordinate" || { echo "$name: not coordinate-sorted"; exit 1; }
  set +o pipefail  # samtools view upstream stops on SIGPIPE when awk exits after chr1
  /usr/bin/time -f "$name extract %es" bash -c "samtools view -h -@ 3 '$bam' | awk '/^@/ {print; next} \$3==\"chr1\" {print; next} {exit}' | samtools view -b -@ 3 -o $B/$name.chr1.tmp.bam -"
  set -o pipefail
  samtools quickcheck $B/$name.chr1.tmp.bam
  mv -f $B/$name.chr1.tmp.bam $B/$name.bam
  samtools index -@ 8 $B/$name.bam
  echo "$name chr1 records $(samtools view -c -@ 8 $B/$name.bam); last $(samtools view $B/$name.bam | tail -1 | cut -f3,4)"
done
date

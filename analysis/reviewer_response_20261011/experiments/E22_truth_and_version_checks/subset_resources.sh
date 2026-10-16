#!/bin/bash
# E22 check 1, step 1: read-only subsets of the resources at the somatic truth positions of every set in sets.tsv,
# plus each normal's genome-wide germline records and gnomAD chr20 (the germline baseline). Outputs go to $W.
#   bash subset_resources.sh            here (4 parallel bcftools)
#   bash subset_resources.sh --slurm    as one Slurm job (-p general -c 4); log in tmp/reviewer_cpu_20261013/E22_truth_and_version_checks/
set -euo pipefail
HERE=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E22_truth_and_version_checks
T=/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E22_truth_and_version_checks
W=/scratch/jshen/data/reviewer_experiments_20261013/E22_truth_and_version_checks
BCF=/wanglab/jshen/anaconda3/envs/sambcfenv/bin/bcftools
PON=/scratch/jshen/data/Pansoma/panel_of_normal_VCFs
GNOMAD=$PON/af-only-gnomad.hg38.vcf.gz
COLORS=$PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz
if [ "${1:-}" = --slurm ]; then
  mkdir -p $T/logs
  sbatch -p general -c 4 --mem=8G -t 6:00:00 -J E22_subset -o $T/logs/E22_subset-%j.out \
    --wrap "bash $HERE/subset_resources.sh && bash $HERE/run_analysis.sh"
  exit 0
fi
mkdir -p $W
tasks=$W/tasks.txt
: > $tasks
declare -A seen
tail -n +2 $HERE/sets.tsv | while IFS=$'\t' read -r name used truth bed full normal sonly lgerm; do
  echo "$BCF view -R $truth $GNOMAD | $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' > $W/$name.gnomad.tsv.tmp && mv $W/$name.gnomad.tsv.tmp $W/$name.gnomad.tsv" >> $tasks
  echo "$BCF view -R $truth $COLORS | $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' > $W/$name.colorsdb.tsv.tmp && mv $W/$name.colorsdb.tsv.tmp $W/$name.colorsdb.tsv" >> $tasks
  echo "$BCF view -R $truth -s $normal $full | $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\t[%VAF]\t[%AD]\n' > $W/$name.full.tsv.tmp && mv $W/$name.full.tsv.tmp $W/$name.full.tsv" >> $tasks
  if [ "$lgerm" != - ]; then
    echo "$BCF view -R $truth $lgerm | $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > $W/$name.labelgerm.tsv.tmp && mv $W/$name.labelgerm.tsv.tmp $W/$name.labelgerm.tsv" >> $tasks
  fi
done
# genome-wide germline records of each normal (FILTER GERMLINE of the full benchmark VCF), one per line
tail -n +2 $HERE/sets.tsv | cut -f5,6 | sort -u | while IFS=$'\t' read -r full normal; do
  line=$(basename $full _DeepSomatic_HCC1395-model.vcf.gz)
  echo "$BCF query -i 'FILTER=\"GERMLINE\"' -s $normal -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT]\n' $full | gzip -1 > $W/$line.normal_germline.tsv.gz.tmp && mv $W/$line.normal_germline.tsv.gz.tmp $W/$line.normal_germline.tsv.gz" >> $tasks
done
echo "$BCF view -r chr20 -v snps $GNOMAD | $BCF query -f '%POS\t%REF\t%ALT\t%INFO/AF\n' | gzip -1 > $W/gnomad_chr20_snps.tsv.gz.tmp && mv $W/gnomad_chr20_snps.tsv.gz.tmp $W/gnomad_chr20_snps.tsv.gz" >> $tasks
echo "tasks: $(wc -l < $tasks)"
date
# skip tasks whose output exists (re-runnable)
grep -v -F -f <(ls $W | grep -v '\.tmp$' | grep -E '\.(tsv|tsv\.gz)$' | sed 's|^|> '$W'/|;s|$|.tmp|') $tasks > $tasks.todo || true
echo "to do: $(wc -l < $tasks.todo)"
xargs -P 4 -I{} -d '\n' bash -c '{}' < $tasks.todo
date
ls -la $W

#!/bin/bash
# E22 check 1(a): FILTER values of DeepSomatic's benchmark VCFs (mixed = somatic PASS + germline GERMLINE, and
# somatic_only) per line; writes HERE/filter_counts.tsv. Read-only, 1 core, a few minutes.
set -euo pipefail
HERE=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E22_truth_and_version_checks
BCF=/wanglab/jshen/anaconda3/envs/sambcfenv/bin/bcftools
B=/scratch/jshen/data/CASTLE/benchmarks/DeepSomatic_HCC1395-model_benchmark/vcfs
out=$HERE/filter_counts.tsv
echo -e "line\tfile\tFILTER\trecords" > $out
for l in H1437 HCC1954 H2009 HCC1395 HCC1937; do
  for f in $B/${l}_DeepSomatic_HCC1395-model.vcf.gz $B/somatic_only/${l}_DeepSomatic_HCC1395-model_somatic-only.vcf.gz; do
    $BCF query -f '%FILTER\n' $f | sort | uniq -c | awk -v L=$l -v F=$(basename $f) '{print L"\t"F"\t"$2"\t"$1}' >> $out
  done
done
cat $out

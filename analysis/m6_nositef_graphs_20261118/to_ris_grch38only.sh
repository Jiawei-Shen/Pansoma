#!/bin/bash
# (wangcluster) Copy to storage2 the chr1 rtg records of the nine GRCh38-only six-tumor SNV call sets
# (m4m6_graphs_20261110/multicancer6_<P>_SNV_linear_nositef/vcf_chr1_<TUMOR>/eval_bed: report.json and the
# raw_calls / pon_calls records), a few MB each. RIS has only their counts per stratum (per_callset_strata.tsv); the
# paired block-bootstrap intervals of Response Tables 1.5 = 2.5 = 3.15 and 3.1 against the HPRC models without the
# site features (regions.py) resample the GRCh38-only and HPRC records together, so they need the records.
# Written on RIS without access to /scratch: a directory that does not exist is reported and skipped.
# usage: [RIS=user@host] bash analysis/m6_nositef_graphs_20261118/to_ris_grch38only.sh
set -euo pipefail
RIS=${RIS:-jiawei.s@compute1-client-1.ris.wustl.edu}
A=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs
R=/scratch/jshen/data/pansoma_net_v2_runs
src=()
for p in Illumina PacBio ONT; do for t in HG008T COLO829T HCC1395T; do
  d=m4m6_graphs_20261110/multicancer6_${p}_SNV_linear_nositef/vcf_chr1_$t
  if [ -s "$R/$d/eval_bed/report.json" ]; then src+=("$R/./$d"); else echo "NOT FOUND, skipped: $R/$d"; fi
done; done
rsync -aRv --include '*/' --include 'report.json' --include '*.stats.json' \
  --include 'raw_calls/*.vcf.gz' --include 'pon_calls/*.vcf.gz' --exclude '*' "${src[@]}" "$RIS:$A/"
echo "copied ${#src[@]} of 9 call sets to $RIS:$A/m4m6_graphs_20261110"

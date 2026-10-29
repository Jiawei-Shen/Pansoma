#!/bin/bash
# tables.md of the ONT control / readnorm pair: per purity, PASS calls of chr1 scored by vcf.sh (rtg vcfeval, benchmark
# BED, before (raw) and after (pon) population filtering), from each run's vcf_chr1_<P>/eval_bed/report.txt.
# usage: bash collect.sh > tables.md
set -euo pipefail
R=/scratch/jshen/data/pansoma_net_v2_runs
echo "# E14b: COLO829T ONT purity joint SNV model, control vs per-read BQ centring (chr1, BED, PASS)"
for set in pon raw; do
  if [ $set = pon ]; then h="after population filtering"; else h="before population filtering"; fi
  echo; echo "## $h ($set)"; echo
  echo "| model | purity | P | R | F1 | best F1 | ceiling |"; echo "|---|---|---|---|---|---|---|"
  for m in control readnorm; do
    for p in P100 P75 P50 P25; do
      awk -v s=$set -v m=$m -v p=$p '$1 == s {printf "| %s | %s | %s | %s | %s | %s | %s |\n", m, p, $6, $7, $8, $10, $NF}' \
        $R/COLO829T_ONT_SNV_joint_${m}_tq/vcf_chr1_$p/eval_bed/report.txt
    done
  done
done
echo; echo "## Training"; echo
for m in control readnorm; do echo "- $m: $(grep -h '^done; best epoch' $R/COLO829T_ONT_SNV_joint_${m}_tq/train.log | sed 's/ -> .*//')"; done

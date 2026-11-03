#!/bin/bash
# ClairS-TO's COLO829T calls (P100 and the P75 / P50 / P25 in-silico mixtures, three platforms), SNV and INDEL, scored on
# RIS with E8's protocol into $PANSOMA_RUNS/comparators/clairs_to/<SET>/<KIND>/. The P100 rows repeat E8's scoring
# (a check of the RIS paths); the mixtures are new. No DeepSomatic mixture VCF is on storage2.
# Run from the repo root: bash analysis/indel_ppr8_9919d4a_20261026/comparators/submit.sh
set -euo pipefail
OUT=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs/comparators/clairs_to
SRC=/storage2/fs1/epigenome/Active/jiawei.s/outputs
for P in Illumina fiberseq ONT; do
  case $P in Illumina) d=WGS;; fiberseq) d=Pacbio;; ONT) d=ONT;; esac
  for F in P100 P75 P50 P25; do
    if [ $F = P100 ]; then set=COLO829T_$P dir=$SRC/clairs-to_COLO829T_$d
    else set=COLO829T_${P}_$F dir=$SRC/clairs-to_INDEL_COLO829_${d}_chr1_$F; fi
    for K in SNV INDEL; do
      ris/submit.sh -g 0 -n 4 -m 16 pnv2_cmp_clairs_to_${set}_$K ../analysis/indel_ppr8_9919d4a_20261026/comparators/score.sh \
          clairs_to "$set" $K "$OUT/$set/$K" "$dir/snv.vcf.gz" "$dir/indel.vcf.gz"
    done
  done
done

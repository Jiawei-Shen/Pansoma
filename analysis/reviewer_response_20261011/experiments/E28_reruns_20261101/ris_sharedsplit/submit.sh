#!/bin/bash
# E28 (3) on RIS (2026-11-03): R3.1 shared validation split, 18 SNV models, seed 0; the command the wangcluster session
# gave, with the 4-GPU SNV resources of the 7fc7396 / normbq runs (SLOTS=16 MEM=192 VAL_WORKERS=3; the GRCh38-only runs
# peaked at 147-149 GB, above snv_pipeline.sh's default 128). snv_pipeline.sh submits vcf.sh after each run.
# Run from the repo root. Job ids: submit_log.txt.
set -euo pipefail
for t in HG008 COLO829; do for p in Illumina PacBio ONT; do
  m=$p; [ $t = COLO829 ] && [ $p = PacBio ] && m=fiberseq
  for g in linear d9 d46; do
    case $g in linear) S=${t}T_${m}_linear;; d9) S=${t}T_${m};; d46) S=${t}T_${m}_d46;; esac
    GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 ris/snv_pipeline.sh ${t}_${m}_SNV_${g}_valchr20-22 $S \
      --norm train --seed 0 --val-chroms chr20 chr21 chr22
  done
done; done

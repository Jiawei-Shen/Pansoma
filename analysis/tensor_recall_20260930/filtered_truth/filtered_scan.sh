#!/bin/bash
# One grep per sample and kind over <sample>/tensors/<kind>/filtered_candidates.ndjson (patterns from patterns.py).
#SBATCH -J filt_scan
#SBATCH -p general
#SBATCH -c 12
#SBATCH --mem=8G
#SBATCH -t 6:00:00
#SBATCH -o filtered_scan_%j.out
set -uo pipefail
cd /scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930/filtered_truth
for s in HG008T_PacBio HG008T_ONT HG008T_Illumina COLO829T_Illumina COLO829T_fiberseq COLO829T_ONT; do
  for k in SNV INDEL; do
    ( LC_ALL=C grep -F -f $s.$k.patterns /scratch/jshen/data/pansoma_v2_tensors/$s/tensors/$k/filtered_candidates.ndjson > $s.$k.hits.ndjson
      echo "$s $k exit $? $(wc -l < $s.$k.hits.ndjson) lines $(date +%T)" ) &
  done
done
wait

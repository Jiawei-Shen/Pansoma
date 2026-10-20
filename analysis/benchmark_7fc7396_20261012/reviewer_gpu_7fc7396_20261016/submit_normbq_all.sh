#!/bin/bash
# The reported models with --norm dataset as of f344bca (BQ only per dataset, statistics of all chromosomes; the
# default), submitted 2026-10-19 (user: "seeds with the new --norm too, and train, test and update whatever else needs
# it"). Every run: ris/snv_pipeline.sh (train, predict chr1, vcf.sh) on the epigenome H100s, SNV on 4 GPUs (DDP,
# --batch-size 1024 per GPU), INDEL on 1. In queue order:
#   A1. d46 SNV seeds 1, 2 (seed 0 = <S>_<P>_SNV_d46_normbq of submit_normbq.sh)                     12 runs
#   A2. d9 SNV seeds 0, 1, 2 (4 GPUs, the matched recipe)                                             18
#   C.  no ch6 on d46 (--drop-planes path_count): SNV seeds 0, 1, 2; INDEL                            24
#       (2026-10-20 13:45, the user: no seeds for no ch6; the 12 SNV seed-1 / seed-2 runs and their vcf jobs were
#       cancelled while pending, so C is SNV seed 0 + INDEL, 12 runs, and the batch 76 runs)
#   B.  INDEL d9, d46 (but the two Illumina of submit_normbq.sh), linear                              16
#   A3. GRCh38-only (linear) SNV seeds 0, 1, 2 (the slowest)                                          18
# Names: <S>_<P>_<KIND>_<arm>_normbq[_s1|_s2], arm d9 / d46 / linear; no ch6 <S>_<P>_<KIND>_nopc_d46_normbq[_sK].
# Run from the repository root on the login node: bash <this file>
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
SETS=(HG008:Illumina HG008:PacBio HG008:ONT COLO829:Illumina COLO829:fiberseq COLO829:ONT)
SNV4=(env GPUS=4 SLOTS=16 MEM=192 VAL_WORKERS=3 KIND=SNV ris/snv_pipeline.sh)
INDEL1=(env GPUS=1 SLOTS=5 MEM=192 VAL_WORKERS=6 KIND=INDEL ris/snv_pipeline.sh)
NORM=(--norm dataset)
set_of() { case $2 in d9) echo "${1}T_$3";; d46) echo "${1}T_$3_d46";; linear) echo "${1}T_$3_linear";; esac; }
seed_args() { [ "$1" = 0 ] && echo "" || echo "--seed $1"; }
suffix() { [ "$1" = 0 ] && echo "" || echo "_s$1"; }

# A1. d46 SNV seeds 1, 2
for seed in 1 2; do for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
  "${SNV4[@]}" ${S}_${P}_SNV_d46_normbq$(suffix $seed) $(set_of $S d46 $P) "${NORM[@]}" $(seed_args $seed)
done; done
# A2. d9 SNV seeds 0, 1, 2
for seed in 0 1 2; do for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
  "${SNV4[@]}" ${S}_${P}_SNV_d9_normbq$(suffix $seed) $(set_of $S d9 $P) "${NORM[@]}" $(seed_args $seed)
done; done
# C. no ch6 on d46
for seed in 0 1 2; do for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
  "${SNV4[@]}" ${S}_${P}_SNV_nopc_d46_normbq$(suffix $seed) $(set_of $S d46 $P) "${NORM[@]}" --drop-planes path_count \
      $(seed_args $seed)
done; done
for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
  "${INDEL1[@]}" ${S}_${P}_INDEL_nopc_d46_normbq $(set_of $S d46 $P) "${NORM[@]}" --drop-planes path_count
done
# B. INDEL d9, d46, linear
for arm in d9 d46 linear; do for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
  [ "$arm:$P" = d46:Illumina ] && continue  # <S>_Illumina_INDEL_d46_normbq: submit_normbq.sh
  "${INDEL1[@]}" ${S}_${P}_INDEL_${arm}_normbq $(set_of $S $arm $P) "${NORM[@]}"
done; done
# A3. linear SNV seeds 0, 1, 2
for seed in 0 1 2; do for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
  "${SNV4[@]}" ${S}_${P}_SNV_linear_normbq$(suffix $seed) $(set_of $S linear $P) "${NORM[@]}" $(seed_args $seed)
done; done

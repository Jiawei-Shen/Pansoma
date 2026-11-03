#!/bin/bash
# The RIS INDEL experiments rerun with the INDEL pipeline of 9919d4a ("ppr8": labels that follow the benchmark, five
# population site features, best epoch by the post-PoN validation F1, PASS at the post-PoN validation recall 0.8 of
# the callable positives; ris/recipes.sh), submitted 2026-10-26. The user: keep each experiment's design and swap the
# recipe; also no ch6, the cross-tumor predictions and the post-PoN best-F1 scoring. In queue order:
#   S. site features of the 27 INDEL sets (CPU, general): HG008T / COLO829T x platform x d9 / d46 / GRCh38-only, and
#      the COLO829T P75 / P50 / P25 dilutions (d9), into $PANSOMA_RUNS/site_features              27 jobs
#   A. seed 0, one model per set: 6 sets x d46 / d9 / linear (1 H100 each), HG008T Illumina d46 first     18 runs
#   B. no ch6 on d46 (--drop-planes path_count), seed 0                                                  6
#   C. COLO829T purity, one model per platform on P100 + P75 + P50 + P25 (d9), tested per purity, GPUs / slots / memory
#      of the 7fc7396 joint INDEL runs                                                                    3
#   D. seeds 1 and 2 of A                                                                                36
#   E. cross-tumor: each d46 model (seeds 0-2) predicts chr1 of the other tumor of its platform (Illumina <->
#      Illumina, PacBio <-> fiberseq, ONT <-> ONT) into <run>/test_chr1_<SET>/                            18 predictions
# Memory: the general queue had ~29,000 pending jobs on 2026-10-26 and 32 GB requests did not start; the jobs submitted
# then were changed in place (bmod with the same -a docker(...), -G and LSF_DOCKER_VOLUMES as ris/submit.sh) to 16 GB
# (site features of ONT sets 32 GB), the values below; the vcf.sh jobs that snv_pipeline.sh / joint_pipeline.sh submit
# (32 GB there) were changed to 16 GB the same way. vcf.sh needs rtg 4 GB plus graph_vcf.
# Every run / prediction is scored twice by vcf.sh: at the checkpoint's call rule (post-PoN recall 0.8) into
# vcf_chr1[_SET]/, and at the post-PoN best-F1 threshold (--threshold postpon) into vcf_chr1[_SET]_postpon/.
# Names: <S>_<P>_INDEL_<arm>_ppr8[_s1|_s2] (arm d46 / d9 / linear), <S>_<P>_INDEL_nopc_d46_ppr8,
# COLO829_<P>_INDEL_joint_ppr8. --norm dataset (the default): every model here is new, nothing pairs with a
# pre-2026-10-13 checkpoint.
# Run from the repository root on the login node: bash analysis/indel_ppr8_9919d4a_20261026/submit.sh [S|A|B|C|D|E ...]
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
STEPS=${*:-S A B C D E}
SETS=(HG008:Illumina HG008:PacBio HG008:ONT COLO829:Illumina COLO829:fiberseq COLO829:ONT)
ARMS=(d46 d9 linear)
INDEL1=(env GPUS=1 SLOTS=5 MEM=192 VAL_WORKERS=6 KIND=INDEL ris/snv_pipeline.sh)
set_of() { case $2 in d9) echo "${1}T_$3";; d46) echo "${1}T_$3_d46";; linear) echo "${1}T_$3_linear";; esac; }
seed_args() { [ "$1" = 0 ] && echo "" || echo "--seed $1"; }
suffix() { [ "$1" = 0 ] && echo "" || echo "_s$1"; }
sitef() { echo "done(pnv2_sitef_$1)"; }
has() { [[ " $STEPS " == *" $1 "* ]]; }
mem_of() { case $1 in *ONT*) echo 32;; *) echo 16;; esac; }  # GB of a site-feature build
# the post-PoN best-F1 scoring of a run (or of one set of a joint run / cross prediction: --set SET) after JOB is done
vcf_postpon() {  # vcf_postpon JOB RUN [SET]
  if [ -n "${3:-}" ]; then
    ris/submit.sh -g 0 -n 4 -m 16 -w "done($1)" "pnv2_vcfpp_${2}_$3" vcf.sh --set "$3" "$2" INDEL postpon _postpon
  else
    ris/submit.sh -g 0 -n 4 -m 16 -w "done($1)" "pnv2_vcfpp_$2" vcf.sh "$2" INDEL postpon _postpon
  fi
}
one() {  # one NAME SET [train options]: train + predict chr1 (+ vcf.sh) and the post-PoN best-F1 scoring
  local name=$1 set=$2; shift 2
  WAIT=$(sitef "$set") "${INDEL1[@]}" "$name" "$set" "$@"
  vcf_postpon "pnv2_$name" "$name"
}

# S. site features (HG008T Illumina d46 first; its test job 425052 was submitted by hand)
if has S; then
  for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
    for set in $(set_of $S d46 $P) $(set_of $S d9 $P) $(set_of $S linear $P); do
      [ "$set" = HG008T_Illumina_d46 ] && continue
      ris/submit.sh -g 0 -n 4 -m $(mem_of "$set") "pnv2_sitef_$set" site_features.sh "$set"
    done
  done
  for P in Illumina fiberseq ONT; do for d in P75 P50 P25; do
    ris/submit.sh -g 0 -n 4 -m $(mem_of $P) "pnv2_sitef_COLO829T_${P}_$d" site_features.sh "COLO829T_${P}_$d"
  done; done
fi
# A. seed 0
if has A; then
  for arm in "${ARMS[@]}"; do for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
    one ${S}_${P}_INDEL_${arm}_ppr8 $(set_of $S $arm $P)
  done; done
fi
# B. no ch6 on d46
if has B; then
  for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
    one ${S}_${P}_INDEL_nopc_d46_ppr8 $(set_of $S d46 $P) --drop-planes path_count
  done
fi
# C. COLO829T purity joint models (d9)
if has C; then
  for P in Illumina fiberseq ONT; do
    case $P in Illumina) G=(GPUS=1 SLOTS=5 MEM=128);; fiberseq) G=(GPUS=2 SLOTS=14 MEM=384);; ONT) G=(GPUS=2 SLOTS=14 MEM=512);; esac
    sets=COLO829T_$P,COLO829T_${P}_P75,COLO829T_${P}_P50,COLO829T_${P}_P25
    w=$(for s in ${sets//,/ }; do printf '%s && ' "$(sitef $s)"; done); w=${w% && }
    name=COLO829_${P}_INDEL_joint_ppr8
    env "${G[@]}" WAIT="$w" VAL_WORKERS=3 KIND=INDEL ris/joint_pipeline.sh $name "$sets"
    for s in ${sets//,/ }; do vcf_postpon "pnv2_$name" $name $s; done
  done
fi
# D. seeds 1 and 2
if has D; then
  for seed in 1 2; do for arm in "${ARMS[@]}"; do for SP in "${SETS[@]}"; do S=${SP%%:*} P=${SP#*:}
    one ${S}_${P}_INDEL_${arm}_ppr8$(suffix $seed) $(set_of $S $arm $P) $(seed_args $seed)
  done; done; done
fi
# E. cross-tumor predictions of the d46 models
if has E; then
  for seed in 0 1 2; do for pair in Illumina:Illumina PacBio:fiberseq ONT:ONT; do H=${pair%%:*} C=${pair#*:}
    for dir in "HG008_${H}:COLO829T_${C}_d46" "COLO829_${C}:HG008T_${H}_d46"; do
      name=${dir%%:*}_INDEL_d46_ppr8$(suffix $seed) set=${dir#*:}
      job=pnv2_xpred_${name}_$set
      ris/submit.sh -q epigenome -s gmodel=NVIDIAH10080GBHBM3 -n 4 -m 64 -w "done(pnv2_$name) && $(sitef $set)" \
          "$job" predict_set.sh "$name" "$set" INDEL
      ris/submit.sh -g 0 -n 4 -m 16 -w "done($job)" "pnv2_vcf_${name}_$set" vcf.sh --set "$set" "$name" INDEL
      vcf_postpon "$job" "$name" "$set"
    done
  done; done
fi

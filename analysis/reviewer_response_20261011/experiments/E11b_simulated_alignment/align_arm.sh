#!/bin/bash
# E11b step 3: align the simulated reads to one reference arm with vg giraffe, convert the GAM to GAF (graph-space
# view) and surject it to GRCh38 (BAM). Usage: align_arm.sh ARM PLATFORM[,PLATFORM...]
#   ARM: linear (GRCh38-only graph, vg 1.66) | d9 (HPRC v1.1 d9, vg 1.65) | d46 (HPRC v2.1 d46, vg 1.77)
#   PLATFORM: illumina (paired, short-read minimizers) | hifi (-b hifi) | ont (-b r10), long-read minimizers
# Indexes and vg versions are those of the production GAMs (memory notes graph-ablation-gams / hprc-v2-graph):
# the linear and d9 minimizer files are format v10 (vg <= 1.66), the d46 ones v11 (vg 1.77 rebuild).
set -euo pipefail
ARM=$1
PLATFORMS=$2
THREADS=${SLURM_CPUS_PER_TASK:-16}
SIM=/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment/sim
OUT=/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment/aln/${ARM}
mkdir -p "${OUT}"

case "${ARM}" in
    linear) VG=/scratch/wzhang/bin/vg_v1.66.0; P=/scratch/jshen/data/Linear_Reference_VG_Indexes/hg38_linear
            GBZ=${P}.giraffe.gbz; DIST=${P}.dist; SURJ_REF=(); SURJ_LEN=0 ;;
    d9)     VG=/scratch/wzhang/bin/vg_v1.65.0; P=/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9
            GBZ=${P}.gbz; DIST=${P}.dist; SURJ_REF=(-n GRCh38); SURJ_LEN=0 ;;
    d46)    VG=/scratch/jshen/bin/vg_v1.77.0; P=/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/hprc-v2.1-mc-grch38.d46
            GBZ=${P}.gbz; DIST=${P}.dist; SURJ_REF=(-n GRCh38); SURJ_LEN=1 ;;
    chm13)  # control: CHM13-only linear graph built like the GRCh38-only arm (chm13_index.sbatch); truth = CHM13 lift
            VG=/scratch/wzhang/bin/vg_v1.66.0; P=/scratch/jshen/data/reviewer_experiments_20261013/E11b_simulated_alignment/chm13_linear/chm13_linear
            GBZ=${P}.giraffe.gbz; DIST=${P}.dist; SURJ_REF=(); SURJ_LEN=0 ;;
    *) echo "unknown arm ${ARM}" >&2; exit 2 ;;
esac
echo "[INFO] host $(hostname) job ${SLURM_JOB_ID:-} arm ${ARM}: $(${VG} version | head -1)"

for PLAT in ${PLATFORMS//,/ }; do
    case "${PLAT}" in
        illumina) MODE=shortread; READS=(-f "${SIM}/illumina/sim_R1.fq.gz" -f "${SIM}/illumina/sim_R2.fq.gz"); PRESET=(); PAIR=(-i); RL=short ;;
        hifi)     MODE=longread; READS=(-f "${SIM}/hifi/sim.fq.gz"); PRESET=(-b hifi); PAIR=(); RL=long ;;
        ont)      MODE=longread; READS=(-f "${SIM}/ont/sim.fq.gz"); PRESET=(-b r10); PAIR=(); RL=long ;;
        *) echo "unknown platform ${PLAT}" >&2; exit 2 ;;
    esac
    MIN=${P}.${MODE}.withzip.min; ZIP=${P}.${MODE}.zipcodes
    GAM=${OUT}/${PLAT}.gam
    if [[ ! -s ${GAM} ]]; then
        CMD=("${VG}" giraffe -Z "${GBZ}" -d "${DIST}" -m "${MIN}" -z "${ZIP}" "${READS[@]}" "${PRESET[@]}" -t "${THREADS}")
        echo "[INFO] $(date +%T) ${CMD[*]} > ${GAM}"
        /usr/bin/time -v "${CMD[@]}" > "${GAM}.tmp" 2> "${OUT}/${PLAT}.giraffe.log"
        mv "${GAM}.tmp" "${GAM}"
    fi
    if [[ ! -s ${OUT}/${PLAT}.gaf.gz ]]; then
        echo "[INFO] $(date +%T) ${VG} convert -G ${GAM} ${GBZ} | gzip"
        "${VG}" convert -t "${THREADS}" -G "${GAM}" "${GBZ}" | gzip -1 > "${OUT}/${PLAT}.gaf.gz.tmp"
        mv "${OUT}/${PLAT}.gaf.gz.tmp" "${OUT}/${PLAT}.gaf.gz"
    fi
    if [[ ! -s ${OUT}/${PLAT}.bam ]]; then
        SL=(); [[ ${SURJ_LEN} == 1 ]] && SL=(-D "${RL}")
        CMD=("${VG}" surject -x "${GBZ}" "${SURJ_REF[@]}" "${SL[@]}" "${PAIR[@]}" -t "${THREADS}" -b "${GAM}")
        echo "[INFO] $(date +%T) ${CMD[*]} > ${OUT}/${PLAT}.bam"
        /usr/bin/time -v "${CMD[@]}" > "${OUT}/${PLAT}.bam.tmp" 2> "${OUT}/${PLAT}.surject.log"
        mv "${OUT}/${PLAT}.bam.tmp" "${OUT}/${PLAT}.bam"
    fi
    grep -E "Maximum resident|Elapsed" "${OUT}/${PLAT}.giraffe.log" "${OUT}/${PLAT}.surject.log" || true
done
echo "[INFO] $(date +%T) done"

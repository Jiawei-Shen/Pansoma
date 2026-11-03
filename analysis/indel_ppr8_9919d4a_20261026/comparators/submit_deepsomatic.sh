#!/bin/bash
# DeepSomatic tumor-only where no call set exists: the COLO829T P75 / P50 / P25 in-silico mixtures (the BAMs ClairS-TO
# ran on), three platforms, one H100 each on the epigenome queue (run_deepsomatic.sh), then the SNV / INDEL scoring
# (score.sh, CPU) into $PANSOMA_RUNS/comparators/deepsomatic/<SET>/<KIND>/. The user's P100 runs
# (outputs/deepsomatic_COLO829T_*, E8's inputs) are scored as they are, a check against E8.
# Run from the repo root: bash analysis/indel_ppr8_9919d4a_20261026/comparators/submit_deepsomatic.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
RUNS=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs
OUT=$RUNS/comparators/deepsomatic
MIX=/storage2/fs1/epigenome/Active/jiawei.s/Dataset/COLO829T_Tumor_Normal_inSilico_Mix
SRC=/storage2/fs1/epigenome/Active/jiawei.s/outputs
score() {  # SET VCF [WAIT]
  for K in SNV INDEL; do
    ris/submit.sh -g 0 -n 4 -m 16 ${3:+-w "$3"} pnv2_cmp_deepsomatic_${1}_$K \
        ../analysis/indel_ppr8_9919d4a_20261026/comparators/score.sh deepsomatic "$1" $K "$OUT/$1/$K" "$2"
  done
}
for P in Illumina fiberseq ONT; do
  case $P in
    Illumina) d=WGS bam=COLO829_WGS/MIX_%s_SMHTCOLO829.mixed.sorted.bam model=WGS_TUMOR_ONLY;;
    fiberseq) d=Pacbio bam=COLO829_Pacbio/MIX_%s_SMHTCOLO829_fiberseq.mixed.sorted.bam model=PACBIO_TUMOR_ONLY;;
    ONT) d=ONT bam=COLO829_ONT/MIX_%s_SMHTCOLO829_ONT.mixed.sorted.bam model=ONT_TUMOR_ONLY;;
  esac
  score COLO829T_$P "$SRC/deepsomatic_COLO829T_$d/COLO829T_$d.chr1.vcf.gz"
  for F in P75 P50 P25; do
    set=COLO829T_${P}_$F f=$(echo $F | tr P p)
    LSF_DOCKER_VOLUMES="/storage2/fs1/epigenome/Active:/storage2/fs1/epigenome/Active" LSF_DOCKER_SHM_SIZE=16g \
      bsub -G compute-epigenome-condo -q epigenome -J pnv2_ds_$set -g /jiawei.s/pansoma \
        -gpu "num=1:gmodel=NVIDIAH10080GBHBM3:j_exclusive=yes" -Q 99 -n 16 \
        -R "gpuhost span[hosts=1] rusage[mem=192GB]" -M 192GB -a 'docker(google/deepsomatic:1.9.0-gpu)' \
        -oo "$RUNS/logs/pnv2_ds_$set.log" bash "$HERE/run_deepsomatic.sh" "$set" "$MIX/$(printf "$bam" $f)" $model "$OUT/$set/calls"
    score $set "$OUT/$set/calls/$set.chr1.vcf.gz" "done(pnv2_ds_$set)"
  done
done

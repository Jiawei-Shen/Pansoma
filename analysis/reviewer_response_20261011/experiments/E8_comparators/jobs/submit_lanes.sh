#!/bin/bash
# E8: submit the ClairS-TO / DeepSomatic chr1 jobs in 6 lanes (each lane a chain: at most 6 E8 jobs run at once).
# Lane jobs follow each other with afterany; prerequisites (image pull 387136, BAM index, bwa) with afterok;
# --kill-on-invalid-dep=yes cancels a job whose prerequisite failed, so its lane moves on.
# usage: bash submit_lanes.sh   (once; writes lanes.tsv)
set -euo pipefail
D=/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E8_comparators
TSV=/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E8_comparators/datasets.tsv
IMG=387136; BWA_COLO=387161; MERGE_HCC=387163; IDX_HCC=387164
cd $D
field() { awk -F'\t' -v d=$1 -v c=$2 '$1==d{print $c}' $TSV; }
sub() {  # NAME SCRIPT DEPENDENCY ARGS...
  local name=$1 script=$2 dep=$3; shift 3
  local j; j=$(sbatch --parsable -J $name -o $D/logs/%x-%j.out --kill-on-invalid-dep=yes --dependency=$dep $script "$@")
  echo -e "$j\t$name\t$dep" >> $D/lanes.tsv; echo $j
}
ct() {  # DATASET DEPENDENCY -> job id
  local extra; extra=$(field $1 6); [ "$extra" = "-" ] && extra=""
  sub e8_ct_$1 $D/run_clairs_to.sh $2 $1 $(field $1 4) $(field $1 5) $(field $1 3) $extra
}
ds() { sub e8_ds_$1 $D/run_deepsomatic.sh $2 $1 $(field $1 4) $(field $1 7) $(field $1 3); }
echo -e "job\tname\tdependency" > $D/lanes.tsv

# lane D: ClairS-TO HG008T x3, COLO829T fiberseq / ONT, then the purity mixtures
j=$(ct HG008T_Illumina afterok:$IMG)
for d in HG008T_PacBio HG008T_ONT COLO829T_fiberseq COLO829T_ONT COLO829T_ONT_P75 COLO829T_ONT_P50 COLO829T_ONT_P25 \
         COLO829T_fiberseq_P75 COLO829T_fiberseq_P50 COLO829T_fiberseq_P25 COLO829T_Illumina_P75 COLO829T_Illumina_P50 COLO829T_Illumina_P25; do
  j=$(ct $d afterany:$j,afterok:$IMG)
done
# lane C: index HCC1395 (running) -> ClairS-TO HCC1395 PacBio / ONT -> index the CASTLE BAMs -> ClairS-TO CASTLE
j=$(ct HCC1395T_PacBio afterok:$IMG:$IDX_HCC)
j=$(ct HCC1395T_ONT afterany:$j,afterok:$IMG:$IDX_HCC)
IDX_CASTLE=$(sub e8_index_CASTLE $D/index_bams.sh afterany:$j \
  HCC1954T_PacBio=/scratch/jshen/data/CASTLE/HCC1954/PacBio/1954_Tumor_HiFi.GRCh38.sorted.bam \
  HCC1954T_ONT=/scratch/jshen/data/CASTLE/HCC1954/ONT/1954_Tumor_ONT.GRCh38.sorted.bam \
  HCC1954T_Illumina=/scratch/jshen/data/CASTLE/HCC1954/Illumina/1954_Illumina_Tumor.GRCh38.addreadgroup.100x.sorted.bam \
  H1437T_ONT=/scratch/jshen/data/CASTLE/H1437/ONT/1437_Tumor_ONT.GRCh38.sorted.bam \
  H1437T_Illumina=/scratch/jshen/data/CASTLE/H1437/Illumina/1437_Illumina_Tumor.GRCh38.addreadgroup.100x.sorted.bam \
  H2009T_PacBio=/scratch/jshen/data/CASTLE/H2009/PacBio/2009_Tumor_HiFi.GRCh38.sorted.bam \
  H2009T_ONT=/scratch/jshen/data/CASTLE/H2009/ONT/2009_Tumor_ONT.GRCh38.sorted.bam \
  H2009T_Illumina=/scratch/jshen/data/CASTLE/H2009/Illumina/2009_Illumina_Tumor.GRCh38.addreadgroup.100x.sorted.bam)
j=$IDX_CASTLE
for d in HCC1954T_PacBio HCC1954T_ONT HCC1954T_Illumina H1437T_PacBio H1437T_ONT H1437T_Illumina H2009T_PacBio H2009T_ONT H2009T_Illumina; do
  j=$(ct $d afterany:$j,afterok:$IMG:$IDX_CASTLE)
done
# lane E: DeepSomatic HG008T PacBio, HG008T Illumina, COLO829T ONT, HCC1395 PacBio, CASTLE
j=$(ds HG008T_PacBio afterok:$IMG)
j=$(ds HG008T_Illumina afterany:$j,afterok:$IMG)
j=$(ds COLO829T_ONT afterany:$j,afterok:$IMG)
j=$(ds HCC1395T_PacBio afterany:$j,afterok:$IMG:$IDX_HCC)
for d in HCC1954T_PacBio H1437T_PacBio H2009T_PacBio HCC1954T_Illumina H1437T_Illumina; do
  j=$(ds $d afterany:$j,afterok:$IMG:$IDX_CASTLE)
done
# lane F: DeepSomatic HG008T ONT, COLO829T fiberseq, HCC1395 ONT, CASTLE
j=$(ds HG008T_ONT afterok:$IMG)
j=$(ds COLO829T_fiberseq afterany:$j,afterok:$IMG)
j=$(ds HCC1395T_ONT afterany:$j,afterok:$IMG:$IDX_HCC)
for d in HCC1954T_ONT H1437T_ONT H2009T_ONT H2009T_Illumina; do
  j=$(ds $d afterany:$j,afterok:$IMG:$IDX_CASTLE)
done
# lane A: bwa COLO829T Illumina (running) -> ClairS-TO, DeepSomatic COLO829T Illumina -> DeepSomatic Illumina mixtures
j=$(ct COLO829T_Illumina afterok:$IMG:$BWA_COLO)
j=$(ds COLO829T_Illumina afterany:$j,afterok:$IMG:$BWA_COLO)
for d in COLO829T_Illumina_P75 COLO829T_Illumina_P50 COLO829T_Illumina_P25; do j=$(ds $d afterany:$j,afterok:$IMG); done
# lane B: bwa + merge HCC1395 Illumina (running) -> ClairS-TO, DeepSomatic -> DeepSomatic long-read mixtures
j=$(ct HCC1395T_Illumina afterok:$IMG:$MERGE_HCC)
j=$(ds HCC1395T_Illumina afterany:$j,afterok:$IMG:$MERGE_HCC)
for d in COLO829T_ONT_P75 COLO829T_ONT_P50 COLO829T_ONT_P25 COLO829T_fiberseq_P75 COLO829T_fiberseq_P50 COLO829T_fiberseq_P25; do
  j=$(ds $d afterany:$j,afterok:$IMG)
done
column -t $D/lanes.tsv

#!/bin/bash
# Seed-0 checkpoints on their own chr1 (check against the RIS scores) and all 18 on HCC1395T chr1 / chr20.
set -euo pipefail
W=/scratch/jshen/Github/Pansoma/tmp/ris_ppr8_hcc1395_20261102
VCF=/scratch/jshen/Github/Pansoma/machine_learning/pansoma_net_v2/jobs/vcf.sh
R=/scratch/jshen/data/pansoma_net_v2_runs/ris_ppr8_20261026
sub() {  # RUN SET CHROM
  local j v
  j=$(sbatch --parsable -J p_$1_$2_$3 -o $W/logs/pred_$1_$2_$3_%j.out $W/predict.sh $1 $2 $3)
  v=$(sbatch --parsable -p general -x maotai,guinness,meow,miao -c 4 --mem=2G -t 6:00:00 -J v_$1_$2_$3 \
      -o $W/logs/vcf_$1_$2_$3_%j.out --dependency=afterok:$j --export=ALL,TEST=test_$3_$2,SET_NAME=$2 \
      $VCF ris_ppr8_20261026/$1 INDEL $R/$1/vcf_$3_$2)
  printf "%s\t%s\t%s\t%s\t%s\n" $j $v $1 $2 $3 | tee -a $W/jobs.tsv
}
for s in HG008 COLO829; do for p in Illumina PacBio ONT; do
  m=$p; [ $s = COLO829 ] && [ $p = PacBio ] && m=fiberseq
  own=${s}T_${m}_d46
  sub ${s}_${m}_INDEL_d46_ppr8 $own chr1
  for sd in "" _s1 _s2; do for c in chr1 chr20; do sub ${s}_${m}_INDEL_d46_ppr8$sd HCC1395T_${p}_d46 $c; done; done
done; done

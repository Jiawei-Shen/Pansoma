#!/bin/bash
# usage: submit_train.sh DATA F RUN [DEPENDENCY] [SEED] [PLAT] [MEM]  -- train.sh then the four tests (vcf_code.sh of tmp/snv_sitefeat_20261103);
# vcf outputs <run>/vcf_<chrom>_<set without _snvaf05>; ids appended to jobs.tsv
set -euo pipefail
DATA=$1 F=$2 RUN=$3 DEP=${4:-} SEED=${5:-0} PLAT=${6:-PacBio} MEM=${7:-40G}
CO=COLO829T_${PLAT}_d46; [ $PLAT = PacBio ] && CO=COLO829T_fiberseq_d46; HG=HG008T_${PLAT}_d46
VMEM=3G; [ $PLAT = ONT ] && VMEM=8G
W=/scratch/jshen/Github/Pansoma/tmp/snv_af05_20261103; V=/scratch/jshen/Github/Pansoma/tmp/snv_sitefeat_20261103/vcf_code.sh
R=/scratch/jshen/data/pansoma_net_v2_runs/$RUN
HCC=HCC1395T_${PLAT}_d46; [ $DATA = af05 ] && HCC=HCC1395T_${PLAT}_d46_snvaf05
TJ=$(sbatch --parsable -J t_$RUN --mem=$MEM ${DEP:+--dependency=afterok:$DEP --kill-on-invalid-dep=yes} -o $W/logs/train_${RUN}_%j.out $W/train.sh $DATA $F $RUN $SEED $PLAT)
line="$(date '+%F %T')	$RUN	train $TJ (plat $PLAT, data $DATA, F $F, seed $SEED, mem $MEM${DEP:+, after $DEP})"
for x in $HCC:chr1:HCC1395T_${PLAT}_d46 $HG:chr1:$HG $CO:chr1:$CO $HCC:chr20:HCC1395T_${PLAT}_d46; do
  IFS=: read s c o <<< "$x"
  j=$(SET_NAME=$s TEST=test_$c sbatch --parsable --export=ALL -J v_${RUN}_${o}_$c --mem=$VMEM --dependency=afterok:$TJ --kill-on-invalid-dep=yes \
      -o $W/logs/vcf_${RUN}_${o}_${c}_%j.out $V $RUN SNV $R/vcf_${c}_$o)
  line="$line vcf_${o}_$c $j"
done
echo "$line" | tee -a $W/jobs.tsv

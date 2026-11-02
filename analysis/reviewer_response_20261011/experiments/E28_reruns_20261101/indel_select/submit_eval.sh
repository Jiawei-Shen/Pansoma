#!/bin/bash
# E28 (2): chr1 / chr20 test predictions + ppr8 VCF and rtg for the reruns whose two selections keep different epochs.
# usage: submit_eval.sh RUN WHICH...   (WHICH: best or excl; RUN e.g. release_s1_Illumina_INDEL)
set -euo pipefail
W=/scratch/jshen/Github/Pansoma/tmp/rerun_20261101/indel_select
VCF=/scratch/jshen/Github/Pansoma/tmp/indel_gap_20261020/jobs/vcf_v3.sh
CODE=/scratch/jshen/Github/Pansoma/tmp/indel_gap_20261020/code_final/machine_learning
RUN=$1; shift; p=${RUN#release_s?_}; p=${p%_INDEL}; co=COLO829T_$p; [ $p = PacBio ] && co=COLO829T_fiberseq
for which in "$@"; do
  ck=best.pth; [ $which = excl ] && ck=best_excl_offpath.pth
  for sc in HCC1395T_${p}_d46:chr1 HCC1395T_${p}_d46:chr20 HG008T_${p}_d46:chr1 ${co}_d46:chr1; do
    s=${sc%:*}; c=${sc#*:}; t=test_${which}_${c}_$s
    j=$(sbatch --parsable -J pred_${which}_${RUN}_${s}_$c $W/predict_ckpt.sh selcheck_20261101/$RUN $ck $s $c $t)
    v=$(sbatch --parsable -p general -x maotai,guinness,meow,miao -c 4 --mem=8G -t 6:00:00 -J vcf_${which}_${RUN}_${s}_$c \
        -o $W/logs/vcf_${which}_${RUN}_${s}_${c}_%j.out --dependency=afterok:$j \
        --export=ALL,CODE=$CODE,THRESHOLD=postpon-recall,TARGET_RECALL=0.8,TEST=$t,SET_NAME=$s \
        $VCF selcheck_20261101/$RUN INDEL /scratch/jshen/data/pansoma_net_v2_runs/selcheck_20261101/$RUN/vcf_ppr8_${which}_${c}_$s)
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" $j $v $RUN $which $s $c | tee -a $W/eval_jobs.tsv
  done
done

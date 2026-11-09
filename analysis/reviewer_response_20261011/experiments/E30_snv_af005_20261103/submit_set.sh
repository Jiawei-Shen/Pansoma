#!/bin/bash
# usage: submit_set.sh SET PARTITION NODE MEM [AFTER]   (AFTER: a job id the prepare waits for, afterany) -- prepare -> run tasks (no finalize) -> fallback (afternotok) -> finalize
# -> site features (afterok finalize) for <SET>_snvaf05; job ids appended to jobs.tsv
set -euo pipefail
S=$1 PART=$2 NODE=$3 MEM=$4 AFTER=${5:-}
W=/scratch/jshen/Github/Pansoma/tmp/snv_af05_20261103; SF=/scratch/jshen/Github/Pansoma/tmp/split_finalize_20261003
R=/scratch/jshen/data/pansoma_v2_tensors/${S}_snvaf05
P=$(/wanglab/jshen/anaconda3/bin/python -c "import json; print(json.load(open('/scratch/jshen/data/pansoma_v2_tensors/$S/run/config.json'))['processes'])")
J1=$(sbatch --parsable -J af05_${S}_prepare ${AFTER:+--dependency=afterany:$AFTER} -p $PART -c 1 --mem=16G -t 6:00:00 -w $NODE -o $W/logs/prepare_${S}_%j.out $W/$S/prepare.sh)
J2=$(sbatch --parsable -J af05_${S}_run --dependency=afterok:$J1 --kill-on-invalid-dep=yes -p $PART -w $NODE -c $P --mem=$MEM -t 3-00:00:00 -o $W/logs/run_${S}_%j.out $SF/run_tasks.sh $R/run)
J3=$(sbatch --parsable -J af05_${S}_fallback --dependency=afternotok:$J2 --kill-on-invalid-dep=yes -p $PART -w $NODE -c $P --mem=$MEM -t 3-00:00:00 -o $W/logs/run_${S}_%j.out $SF/fallback_tasks.sh $R/run)
J4=$(sbatch --parsable -J af05_${S}_finalize --dependency=afterok:$J2?afterok:$J3 -p $PART -c 8 --mem=20G -t 2-00:00:00 -w $NODE -o $W/logs/finalize_${S}_%j.out $SF/finalize.sh $R/run)
J5=$(sbatch --parsable -J af05_${S}_sitef --dependency=afterok:$J4 -p $PART -w $NODE -o $W/logs/sitef_${S}_%j.out $W/build_sitef.sh ${S}_snvaf05)
printf "%s\t%s\tprepare %s run %s (-p %s -w %s -c %s --mem=%s) fallback %s finalize %s sitef %s\n" "$(date '+%F %T')" "$S" $J1 $J2 $PART $NODE $P $MEM $J3 $J4 $J5 | tee -a $W/jobs.tsv

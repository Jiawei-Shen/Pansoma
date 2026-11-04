#!/bin/bash
# usage: watch.sh PLAT TRAINJOB VCFJOBS(comma list): wait until the training fails or every vcf job ends, then compare.py
PLAT=$1; TR=$2; V=$3
cd /scratch/jshen/Github/Pansoma/tmp/snv_sitefeat_20261103
while true; do
  st=$(sacct -n -X -P -j $TR,$V -o JobID,State 2>/dev/null)
  tr=$(echo "$st" | awk -F'|' -v j=$TR '$1==j{print $2}')
  case "$tr" in FAILED*|CANCELLED*|TIMEOUT*|OUT_OF_ME*|NODE_FAIL*|PREEMPTED*) echo "TRAIN $TR $tr"; tail -5 logs/train_${PLAT}_$TR.out; exit 1;; esac
  if [ -n "$st" ] && ! echo "$st" | grep -v "^$TR|" | grep -Eq "\|(PENDING|RUNNING|REQUEUED|CONFIGURING|COMPLETING)"; then break; fi
  sleep 120
done
echo "$st"; grep "^done; best epoch" /scratch/jshen/data/pansoma_net_v2_runs/multicancer6_${PLAT}_SNV_sitef_ppr9/train.log
/wanglab/jshen/anaconda3/bin/python compare.py $PLAT

#!/bin/bash
# Submits the six INDEL models (train_indel_lab.sh) on tequila and, after each, the call sets of chr1 of HCC1395T, HG008T
# and COLO829T and chr20 of HCC1395T (tmp/indel_gap_20261020/jobs/vcf_v3.sh with code_final): PASS ppr8
# (THRESHOLD=postpon-recall TARGET_RECALL=0.8 -> vcf_ppr8_<chrom>_<tumor>); ONT also pp, its PASS rule
# (THRESHOLD=postpon -> vcf_pp_<chrom>_<tumor>). A model whose run has best.pth and both predictions is skipped.
# GPU and time limit per model below (24 GB slice: Illumina 3-11 h, PacBio 20-23 h, ONT 36-37 h for the paired models).
# usage: [DRY=1] [ONLY="PacBio:d9 ..."] bash submit_indel_lab.sh   (login node; ids appended to jobs_indel_lab.tsv)
set -euo pipefail
A=/scratch/jshen/Github/Pansoma/analysis/m6_nositef_graphs_20261118; L=/scratch/jshen/Github/Pansoma/tmp/m6_nositef_indel_20261124/logs
VI=/scratch/jshen/Github/Pansoma/tmp/indel_gap_20261020/jobs/vcf_v3.sh
source /scratch/jshen/Github/Pansoma/tmp/m4m6_graphs_20261110/sets.sh
sb() { if [ "${DRY:-0}" = 1 ]; then echo "DRY$RANDOM"; else sbatch --parsable "$@"; fi; }
# model:gres:time limit
PLAN=${PLAN:-"ONT:d46:h100:2-00:00:00 ONT:d9:24gb:2-00:00:00 PacBio:d46:24gb:1-06:00:00 PacBio:d9:h100:1-06:00:00
  Illumina:d46:h100:14:00:00 Illumina:d9:24gb:14:00:00"}
for x in $PLAN; do
  IFS=: read -r P G GRES TL <<< "$x"
  [ -z "${ONLY:-}" ] || [[ " $ONLY " == *" $P:$G "* ]] || continue
  RUN=m6_nositef_graphs_20261118/multicancer6_${P}_INDEL_${G}_nositef; R=/scratch/jshen/data/pansoma_net_v2_runs/$RUN
  if [ -s $R/best.pth ] && ls $R/test_chr20/*.predictions.ndjson.gz > /dev/null 2>&1; then echo "$RUN exists, skipped"; continue; fi
  TJ=$(sb -J t_m6_${P}_INDEL_${G}_nositef --gres=gpu:$GRES:1 -t $TL -o $L/train_m6_${P}_INDEL_${G}_nositef_%j.out $A/train_indel_lab.sh $P $G)
  line="$RUN train $TJ (gpu:$GRES:1, -t $TL)"
  for x2 in HCC1395T:chr1 HG008T:chr1 COLO829T:chr1 HCC1395T:chr20; do
    t=${x2%:*} c=${x2#*:}
    if [ $G = d9 ]; then s=$(set_name INDEL $P d9 $t); else p=$P; [ $t = COLO829T ] && [ $P = PacBio ] && p=fiberseq; s=${t}_${p}_d46; fi
    rules="ppr8:postpon-recall"; [ $P = ONT ] && rules+=" pp:postpon"
    for r in $rules; do
      n=${r%:*} th=${r#*:}; tr=""; [ $n = ppr8 ] && tr=0.8
      j=$(CODE=/scratch/jshen/Github/Pansoma/tmp/indel_gap_20261020/code_final/machine_learning THRESHOLD=$th TARGET_RECALL=$tr \
          SET_NAME=$s TEST=test_$c sb --export=ALL -J v${n}_m6_${P}_INDEL_${G}_nositef_${t}_$c --mem=8G --dependency=afterok:$TJ --kill-on-invalid-dep=yes \
          -o $L/vcf_${n}_m6_${P}_INDEL_${G}_nositef_${t}_${c}_%j.out $VI $RUN INDEL $R/vcf_${n}_${c}_$t)
      line+=" vcf_${n}_${c}_$t $j"
    done
  done
  echo "$(date '+%F %T')	$line" | tee -a $([ "${DRY:-0}" = 1 ] && echo /dev/null || echo $A/jobs_indel_lab.tsv)
done

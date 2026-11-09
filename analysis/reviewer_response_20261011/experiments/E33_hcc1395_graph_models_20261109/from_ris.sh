#!/bin/bash
# Step 3 of 3 (wangcluster), when the RIS pnv2_pred20_* jobs are done: copy each finished run back (without last.pth)
# into pansoma_net_v2_runs/hcc1395_graphs_af05_20261109/, cancel the held tequila copies, and submit the chr1 / chr20
# VCF jobs with the same vcf_code.sh as the tequila runs (graph VCF -> GRCh38 VCF -> PoN -> vcfeval in the SEQC2 BED).
# A run whose training or chr20 prediction is not finished is skipped; rerun the script later for it.
# usage: [RIS=user@host] bash analysis/reviewer_response_20261011/experiments/E33_hcc1395_graph_models_20261109/from_ris.sh
set -euo pipefail
RIS=${RIS:-jiawei.s@compute1-client-1.ris.wustl.edu}
A=/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs
L=/scratch/jshen/data/pansoma_net_v2_runs/hcc1395_graphs_af05_20261109
W=/scratch/jshen/Github/Pansoma/tmp/hcc1395_linear_d9_20261104/train_af05
V=/scratch/jshen/Github/Pansoma/tmp/snv_sitefeat_20261103/vcf_code.sh
SOCK=$HOME/.ssh/ris_mux_$$
ssh -fNM -S "$SOCK" "$RIS"
trap 'ssh -S "$SOCK" -O exit "$RIS" 2>/dev/null || true' EXIT
E="ssh -S $SOCK"
# held tequila trainings (474803 474807 474811 474815 474819) and their VCF jobs: replaced by the RIS runs
scancel 474803 474804 474805 474807 474808 474809 474811 474812 474813 474815 474816 474817 474819 474820 474821 2>/dev/null || true
for x in HCC1395T_ONT_SNV_linear_normfix:HCC1395T_ONT_linear HCC1395T_Illumina_SNV_d9_lnf25w1:HCC1395T_Illumina_d9_snvaf05 \
         HCC1395T_Illumina_SNV_d46_lnf25w1:HCC1395T_Illumina_d46_snvaf05 HCC1395T_ONT_SNV_d46_normfix:HCC1395T_ONT_d46 \
         HCC1395T_ONT_SNV_d9_normfix:HCC1395T_ONT_d9; do
  N=${x%%:*} S=${x#*:}
  if grep -q "hcc1395_graphs_af05_20261109/$N	from RIS" $W/jobs.tsv; then echo "$N: already copied and submitted"; continue; fi
  if ! $E "$RIS" "grep -q '^done; best epoch' $A/$N/train.log && ls $A/$N/test_chr1/$S.*.SNV.predictions.ndjson.gz $A/$N/test_chr20/$S.*.SNV.predictions.ndjson.gz" > /dev/null 2>&1; then
    echo "$N: not finished on RIS (training or a prediction missing), skipped"; continue
  fi
  rsync -a --exclude last.pth -e "$E" "$RIS:$A/$N" "$L/"
  line="$(date '+%F %T')	hcc1395_graphs_af05_20261109/$N	from RIS (E33 from_ris.sh)"
  for c in chr1 chr20; do
    j=$(SET_NAME=$S TEST=test_$c sbatch --parsable --export=ALL -J v_${N}_$c --mem=8G -o $W/logs/vcf_${N}_${c}_%j.out \
        $V hcc1395_graphs_af05_20261109/$N SNV $L/$N/vcf_$c)
    line="$line vcf_$c $j"
  done
  echo "$line" | tee -a $W/jobs.tsv
done

#!/bin/bash
# Submit the 2026-11-01 reruns (user: run all three). Order on the two tequila H100s: R1.1 masked retraining pairs, then
# the four INDEL release reruns that trained on an H100, then the R3.1 shared-split models. The five INDEL release
# reruns that trained on a 24gb MIG slice run there (at most 3 at a time). Each SNV training gets a CPU VCF + rtg job.
set -euo pipefail
W=/scratch/jshen/Github/Pansoma/tmp/rerun_20261101
VCF=/scratch/jshen/Github/Pansoma/machine_learning/pansoma_net_v2/jobs/vcf.sh
: > $W/jobs.tsv
snv() {  # NAME SET MODE ROOT
  local j v
  j=$(sbatch --parsable -J $1 $W/snv/snv_train.sh $1 $2 $3 0)
  v=$(sbatch --parsable -J vcf_$1 -o $W/snv/logs/vcf_$1_%j.out --dependency=afterok:$j --export=ALL,SET_NAME=$2 $VCF $4/$1 SNV)
  printf "%s\t%s\t%s\t%s\t%s\n" $j $v $1 $2 $3 >> $W/jobs.tsv
}
# R1.1: v2.1 SNV, seed 0, with and without the off-node mask, same 1-GPU setup
for sp in HG008T:HG008:Illumina HG008T:HG008:PacBio HG008T:HG008:ONT COLO829T:COLO829:Illumina COLO829T:COLO829:fiberseq COLO829T:COLO829:ONT; do
  IFS=: read t s p <<< "$sp"
  for m in mask nomask; do snv ${s}_${p}_SNV_d46_offnode_$m ${t}_${p}_d46 $m offnode_retrain_20261101; done
done
# R1.8: INDEL release reruns (args of the originals: micro-batches 1 = H100, 4 = 24gb MIG)
for r in s0_ONT s1_ONT s2_PacBio s2_ONT; do
  seed=${r%%_*}; seed=${seed#s}; p=${r#*_}
  j=$(sbatch --parsable -J sel_release_${r} -p gpu -w tequila --gres=gpu:h100:1 -c 8 --mem=40G -t 3-00:00:00 \
      $W/indel_select/train_select.sh release_${r}_INDEL release $p $seed 1)
  printf "%s\t-\trelease_%s_INDEL\th100\tselect\n" $j $r >> $W/jobs.tsv
done
j=$(sbatch --parsable -J sel_release_mig --array=1-5%3 -p gpu -w tequila --gres=gpu:24gb:1 -c 8 --mem=40G -t 3-00:00:00 \
    -o $W/indel_select/logs/%x_%A_%a.out $W/indel_select/array_mig.sh)
printf "%s\t-\tarray 1-5 (mig_tasks.tsv)\t24gb\tselect\n" $j >> $W/jobs.tsv
# R3.1: shared split (validation chr20-22), seed 0, three references x six datasets
for sp in HG008T:HG008:Illumina HG008T:HG008:PacBio HG008T:HG008:ONT COLO829T:COLO829:Illumina COLO829T:COLO829:fiberseq COLO829T:COLO829:ONT; do
  IFS=: read t s p <<< "$sp"
  snv ${s}_${p}_SNV_linear_valchr20-22 ${t}_${p}_linear shared sharedsplit_20261101
  snv ${s}_${p}_SNV_d9_valchr20-22 ${t}_${p} shared sharedsplit_20261101
  snv ${s}_${p}_SNV_d46_valchr20-22 ${t}_${p}_d46 shared sharedsplit_20261101
done
cat $W/jobs.tsv

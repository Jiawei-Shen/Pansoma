#!/bin/bash
# ONT dilutions without the read-source BQ cue (R1.5 / R2.4 / R3.10 [待补: ONT]), clean pair: control and readnorm,
# same code, hardware, recipe and seed (run_joint.sh). 1. SNV chunk caches chr2-22 of the three ONT dilution sets
# (trained on RIS so far; P100 included, existing chunks are skipped); 2. the two models, one H100 each; 3. vcf.sh per
# set and model (chr1, BED and no BED, PoN, rtg vcfeval) -> <run>/vcf_chr1_<P100|P75|P50|P25>.
set -euo pipefail
H=/scratch/jshen/Github/Pansoma/tmp/ont_readnorm_20261026; R=/scratch/jshen/data/pansoma_net_v2_runs
D=/scratch/jshen/data/pansoma_v2_tensors; P=/wanglab/jshen/anaconda3/bin/python
J=/scratch/jshen/Github/Pansoma/machine_learning/pansoma_net_v2/jobs
deps=""
for s in COLO829T_ONT COLO829T_ONT_P75 COLO829T_ONT_P50 COLO829T_ONT_P25; do
  id=$(KINDS=SNV sbatch --parsable -p general -w ${NODE:-meow} -c 4 --mem=8G -t 6:00:00 -J chunks_$s -o $H/logs/%x-%j.out \
       --export=ALL,KINDS=SNV --wrap "$P $J/build_chunks_from_manifest.py $D/$s/tensors $R/index_cache/chunks")
  deps="$deps:$id"; echo "$(date '+%F %T') chunk cache SNV $s $id" >> $H/jobs.txt
done
for MODE in control readnorm; do
  NAME=COLO829T_ONT_SNV_joint_${MODE}_tq
  id=$(sbatch --parsable -p gpu --gres=gpu:h100:1 -c 8 --mem=72G -t 2-00:00:00 --dependency=afterok$deps \
       -J $NAME -o $H/logs/%x-%j.out $H/run_joint.sh $NAME $MODE)
  echo "$(date '+%F %T') $NAME ($MODE) $id (after$deps)" >> $H/jobs.txt
  for s in COLO829T_ONT COLO829T_ONT_P75 COLO829T_ONT_P50 COLO829T_ONT_P25; do
    sfx=${s##*_}; [ $s = COLO829T_ONT ] && sfx=P100
    v=$(SET_NAME=$s sbatch --parsable --mem=16G --dependency=afterok:$id -J vcf_${NAME}_$sfx \
        -o $H/logs/%x-%j.out --export=ALL,SET_NAME=$s $J/vcf.sh $NAME SNV $R/$NAME/vcf_chr1_$sfx)
    echo "$(date '+%F %T') vcf $NAME $sfx $v (after $id)" >> $H/jobs.txt
  done
done
tail -11 $H/jobs.txt

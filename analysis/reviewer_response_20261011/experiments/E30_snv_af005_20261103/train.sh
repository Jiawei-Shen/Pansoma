#!/bin/bash
# The six-tumor PacBio SNV site-feature model (recipe of tmp/snv_sitefeat_20261103/train_sitef.sh: 50.4 M, 20 epochs,
# pow:0.75, lr 1e-4, batch 1024 on one H100, --site-features, --best-by postpon, --call-recall 0.9) with
# --pon-negative-fraction F (each epoch every somatic and untagged tensor, a fresh fraction F of the PoN-tagged non /
# germline ones, weight 1/F). Code: code_ml/ (git archive of the commit with that option).
# usage: sbatch -J NAME train.sh DATA F OUT
#   DATA af05: HCC1395T / HCC1954T / H1437T / H2009T PacBio_d46_snvaf05 (SNV AF floor 0.05) + HG008T_PacBio_d46 and
#              COLO829T_fiberseq_d46 as they are (user 2026-11-03), site features /scratch/jshen/data/snv_af05_20261103
#   DATA af08: the six sets of multicancer6_PacBio_SNV_sitef_ppr9 (floor 0.08), its site features (the control)
#SBATCH -p gpu
#SBATCH -w tequila
#SBATCH --gres=gpu:h100:1
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 3-00:00:00
set -euo pipefail
DATA=$1 F=$2 OUT=/scratch/jshen/data/pansoma_net_v2_runs/$3
D=/scratch/jshen/data/pansoma_v2_tensors; R=/scratch/jshen/data/pansoma_net_v2_runs; P=/wanglab/jshen/anaconda3/bin/python
if [ $DATA = af05 ]; then
  SF=/scratch/jshen/data/snv_af05_20261103/site_features; HCC=HCC1395T_PacBio_d46_snvaf05
  SETS="HCC1395T_PacBio_d46_snvaf05 HCC1954T_PacBio_d46_snvaf05 H1437T_PacBio_d46_snvaf05 H2009T_PacBio_d46_snvaf05 HG008T_PacBio_d46 COLO829T_fiberseq_d46"
elif [ $DATA = af08 ]; then
  SF=/scratch/jshen/data/snv_sitefeat_20261103/site_features; HCC=HCC1395T_PacBio_d46
  SETS="HCC1395T_PacBio_d46 HCC1954T_PacBio_d46 H1437T_PacBio_d46 H2009T_PacBio_d46 HG008T_PacBio_d46 COLO829T_fiberseq_d46"
else echo "DATA af05 | af08"; exit 1; fi
T=(); for s in $SETS; do T+=($D/$s/tensors); done
cd /scratch/jshen/Github/Pansoma/tmp/snv_af05_20261103/code_ml/machine_learning
mkdir -p $OUT
echo "host $(hostname) job ${SLURM_JOB_ID:-} gpus ${CUDA_VISIBLE_DEVICES:-} code $(cat ../git_head.txt) data $DATA F $F $(date)"; nvidia-smi --query-gpu=index,name,memory.total --format=csv
if ! grep -q "^done; best epoch" $OUT/train.log 2>/dev/null; then
  $P -m pansoma_net_v2.train --tensors "${T[@]}" --kinds SNV --output $OUT --cache-dir $R/index_cache \
    --test-chroms --set-test-chroms $D/$HCC/tensors chr1 chr20 --set-test-chroms $D/HG008T_PacBio_d46/tensors chr1 \
    --set-test-chroms $D/COLO829T_fiberseq_d46/tensors chr1 --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 \
    --class-weights pow:0.75 --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth \
    --site-features $SF --best-by postpon --call-recall 0.9 --pon-negative-fraction $F 2>&1 | grep --line-buffered -v 'UserWarning\|warnings.warn\|x = torch.from_numpy\|isnan'
fi
grep -q "^done; best epoch" $OUT/train.log
$P -m pansoma_net_v2.predict --checkpoint $OUT/best.pth --tensors $D/$HCC/tensors $D/HG008T_PacBio_d46/tensors $D/COLO829T_fiberseq_d46/tensors \
  --kinds SNV --chroms chr1 --output $OUT/test_chr1 --cache-dir $R/index_cache --num-workers 8
$P -m pansoma_net_v2.predict --checkpoint $OUT/best.pth --tensors $D/$HCC/tensors --kinds SNV --chroms chr20 \
  --output $OUT/test_chr20 --cache-dir $R/index_cache --num-workers 8
echo "done $(date)"

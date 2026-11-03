#!/bin/bash
# multicancer6 SNV model with population site features (tmp/snv_sitefeat_20261103): the recipe and data of
# multicancer6_<PLAT>_SNV (jobs/multicancer_20261016) plus --site-features (SNV files of build_sitef.sh),
# --best-by postpon (best epoch by the post-PoN validation somatic F1) and --call-recall 0.9 (PASS at the post-PoN
# validation recall 0.9 of the callable positives, the SNV recall target computed after the PoN); then the chr1 tests
# (HCC1395T, HG008T, COLO829T) and HCC1395T chr20. Code: code/ (ca50ab0 + the SNV site-feature builder).
# usage: sbatch -J NAME train_sitef.sh PLAT     (PLAT Illumina | PacBio | ONT)
#SBATCH -p gpu
#SBATCH -w tequila
#SBATCH --gres=gpu:h100:1
#SBATCH -c 8
#SBATCH --mem=36G
#SBATCH -t 3-00:00:00
set -euo pipefail
PLAT=$1
D=/scratch/jshen/data/pansoma_v2_tensors; R=/scratch/jshen/data/pansoma_net_v2_runs; P=/wanglab/jshen/anaconda3/bin/python
SF=/scratch/jshen/data/snv_sitefeat_20261103/site_features
CO=COLO829T_${PLAT}_d46; [ $PLAT = PacBio ] && CO=COLO829T_fiberseq_d46
OUT=$R/multicancer6_${PLAT}_SNV_sitef_ppr9
T=(); for s in HCC1395T_${PLAT}_d46 HCC1954T_${PLAT}_d46 H1437T_${PLAT}_d46 H2009T_${PLAT}_d46 HG008T_${PLAT}_d46 $CO; do T+=($D/$s/tensors); done
cd /scratch/jshen/Github/Pansoma/tmp/snv_sitefeat_20261103/code/machine_learning
mkdir -p $OUT
echo "host $(hostname) job ${SLURM_JOB_ID:-} gpus ${CUDA_VISIBLE_DEVICES:-} $(date)"; nvidia-smi --query-gpu=index,name,memory.total --format=csv
if ! grep -q "^done; best epoch" $OUT/train.log 2>/dev/null; then
  $P -m pansoma_net_v2.train --tensors "${T[@]}" --kinds SNV --output $OUT --cache-dir $R/index_cache \
    --test-chroms --set-test-chroms $D/HCC1395T_${PLAT}_d46/tensors chr1 chr20 --set-test-chroms $D/HG008T_${PLAT}_d46/tensors chr1 \
    --set-test-chroms $D/$CO/tensors chr1 --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 \
    --class-weights pow:0.75 --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth \
    --site-features $SF --best-by postpon --call-recall 0.9 2>&1 | grep --line-buffered -v 'UserWarning\|warnings.warn\|x = torch.from_numpy\|isnan'
fi
grep -q "^done; best epoch" $OUT/train.log
$P -m pansoma_net_v2.predict --checkpoint $OUT/best.pth --tensors $D/HCC1395T_${PLAT}_d46/tensors $D/HG008T_${PLAT}_d46/tensors $D/$CO/tensors \
  --kinds SNV --chroms chr1 --output $OUT/test_chr1 --cache-dir $R/index_cache --num-workers 8
$P -m pansoma_net_v2.predict --checkpoint $OUT/best.pth --tensors $D/HCC1395T_${PLAT}_d46/tensors --kinds SNV --chroms chr20 \
  --output $OUT/test_chr20 --cache-dir $R/index_cache --num-workers 8
echo "done $(date)"

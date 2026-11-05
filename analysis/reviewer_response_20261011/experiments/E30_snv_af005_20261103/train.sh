#!/bin/bash
# The six-tumor PacBio SNV site-feature model (recipe of tmp/snv_sitefeat_20261103/train_sitef.sh: 50.4 M, 20 epochs,
# pow:0.75, lr 1e-4, batch 1024 on one H100, --site-features, --best-by postpon, --call-recall 0.9) with
# --pon-negative-fraction F (each epoch every somatic and untagged tensor, a fresh fraction F of the PoN-tagged non /
# germline ones, weight 1/F). Code: code_ml/ (git archive of the commit with that option).
# usage: sbatch -J NAME train.sh DATA F OUT [SEED] [PLAT] [TUMORS] [MODE]   (SEED: train --seed, default 0; PLAT PacBio (default) |
#        Illumina | ONT; COLO829T = COLO829T_fiberseq_d46 for PacBio, COLO829T_<PLAT>_d46 otherwise; TUMORS 6 (default) |
#        4: the four CASTLE sets only (multicancer4: HG008T and COLO829T held out entirely, as multicancer4_<PLAT>_SNV);
#        F 1.0 = every negative, the E29 site-feature recipe; MODE pon (default): --pon-negative-fraction F, code_ml/ (96e2d6c);
#        MODE all: --negative-fraction F (a fresh F of all non / germline tensors each epoch), code_ml2/ (629418f);
#        MODE lowaf: --non-fraction F --keep-non-af KEEP (non tensors with AF >= KEEP every epoch, the lower-AF non a
#        fresh F with loss weight 1/F, germline and somatic every epoch; KEEP default 0.100001 = AF <= 0.1 sampled))
#   DATA af05: HCC1395T / HCC1954T / H1437T / H2009T <PLAT>_d46_snvaf05 (SNV AF floor 0.05) + HG008T_<PLAT>_d46 and
#              COLO829T as they are (user 2026-11-03), site features /scratch/jshen/data/snv_af05_20261103
#   DATA af08: the six <PLAT>_d46 sets of multicancer6_<PLAT>_SNV (floor 0.08), site features of tmp/snv_sitefeat_20261103
#SBATCH -p gpu
#SBATCH -w tequila
#SBATCH --gres=gpu:h100:1
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 3-00:00:00
set -euo pipefail
DATA=$1 F=$2 OUT=/scratch/jshen/data/pansoma_net_v2_runs/$3 SEED=${4:-0} PLAT=${5:-PacBio} TUMORS=${6:-6} MODE=${7:-pon}
case $MODE in pon) CODE=code_ml; SAMPLING="--pon-negative-fraction $F";; all) CODE=code_ml2; SAMPLING="--negative-fraction $F";;
  lowaf) CODE=code_ml2; SAMPLING="--non-fraction $F --keep-non-af ${KEEP:-0.100001}";; *) echo "MODE pon | all | lowaf"; exit 1;; esac
CO=COLO829T_${PLAT}_d46; [ $PLAT = PacBio ] && CO=COLO829T_fiberseq_d46; HG=HG008T_${PLAT}_d46
D=/scratch/jshen/data/pansoma_v2_tensors; R=/scratch/jshen/data/pansoma_net_v2_runs; P=/wanglab/jshen/anaconda3/bin/python
if [ $DATA = af05 ]; then
  SF=/scratch/jshen/data/snv_af05_20261103/site_features; HCC=HCC1395T_${PLAT}_d46_snvaf05
  SETS="$HCC HCC1954T_${PLAT}_d46_snvaf05 H1437T_${PLAT}_d46_snvaf05 H2009T_${PLAT}_d46_snvaf05 $HG $CO"
elif [ $DATA = af08 ]; then
  SF=/scratch/jshen/data/snv_sitefeat_20261103/site_features; HCC=HCC1395T_${PLAT}_d46
  SETS="$HCC HCC1954T_${PLAT}_d46 H1437T_${PLAT}_d46 H2009T_${PLAT}_d46 $HG $CO"
else echo "DATA af05 | af08"; exit 1; fi
[ $TUMORS = 4 ] && SETS=$(echo $SETS | cut -d' ' -f1-4)
T=(); for s in $SETS; do T+=($D/$s/tensors); done
HELD="--set-test-chroms $D/$HG/tensors chr1 --set-test-chroms $D/$CO/tensors chr1"; [ $TUMORS = 4 ] && HELD=""
cd /scratch/jshen/Github/Pansoma/tmp/snv_af05_20261103/$CODE/machine_learning
mkdir -p $OUT
echo "host $(hostname) job ${SLURM_JOB_ID:-} gpus ${CUDA_VISIBLE_DEVICES:-} code $(cat ../git_head.txt) plat $PLAT tumors $TUMORS sets $SETS F $F mode $MODE seed $SEED $(date)"; nvidia-smi --query-gpu=index,name,memory.total --format=csv
if ! grep -q "^done; best epoch" $OUT/train.log 2>/dev/null; then
  $P -m pansoma_net_v2.train --tensors "${T[@]}" --kinds SNV --output $OUT --cache-dir $R/index_cache \
    --test-chroms --set-test-chroms $D/$HCC/tensors chr1 chr20 $HELD \
    --depths 3 3 27 3 --dims 96 192 384 768 --epochs 20 --batch-size 1024 --lr 1e-4 \
    --class-weights pow:0.75 --non-fraction 1.0 --block v2 --ignore-reasons residual_partial_somatic_truth allele_partial_somatic_truth \
    --site-features $SF --best-by postpon --call-recall 0.9 $SAMPLING --seed $SEED 2>&1 | grep --line-buffered -v 'UserWarning\|warnings.warn\|x = torch.from_numpy\|isnan'
fi
grep -q "^done; best epoch" $OUT/train.log
$P -m pansoma_net_v2.predict --checkpoint $OUT/best.pth --tensors $D/$HCC/tensors $D/$HG/tensors $D/$CO/tensors \
  --kinds SNV --chroms chr1 --output $OUT/test_chr1 --cache-dir $R/index_cache --num-workers 8
$P -m pansoma_net_v2.predict --checkpoint $OUT/best.pth --tensors $D/$HCC/tensors --kinds SNV --chroms chr20 \
  --output $OUT/test_chr20 --cache-dir $R/index_cache --num-workers 8
echo "done $(date)"

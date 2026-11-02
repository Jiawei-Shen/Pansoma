#!/bin/bash
#SBATCH -p gpu
#SBATCH -w tequila
#SBATCH --gres=gpu:12gb:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -J q6_mask_infer
#SBATCH -o /scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/infer/logs/%x_%A_%a.out
# R1.1 cross-node context, inference only: chr1 re-prediction of a seed-0 checkpoint with (mask) / without (nomask)
# the off-node mask, isolated code copy q6/code (repo pansoma_net_v2 at 90d7c0f + offnode.py), then jobs/vcf.sh.
set -euo pipefail
Q=/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6
read SRC SET MODE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $Q/infer/tasks.tsv)
R=/scratch/jshen/data/pansoma_net_v2_runs; D=/scratch/jshen/data/pansoma_v2_tensors; P=/wanglab/jshen/anaconda3/bin/python
OUTRUN=offnode_mask_20261030/${SRC}__${MODE}; O=$R/$OUTRUN
mkdir -p $O
ln -sfn $R/ris_7fc7396_20261012/$SRC/best.pth $O/best.pth
echo "task ${SLURM_ARRAY_TASK_ID}: $SRC $SET $MODE host $(hostname) job ${SLURM_JOB_ID:-}"
cd $Q/code
export OMP_NUM_THREADS=4; nvidia-smi --query-gpu=index,name,memory.total --format=csv || true
EXTRA=""; [ $MODE = mask ] && EXTRA="--mask-off-node $Q/side"
if [ ! -s $O/test_chr1/$SET.tensors.SNV.predictions.ndjson.gz ]; then
  /usr/bin/time -f "predict %es %MKB" $P -m pansoma_net_v2.predict --checkpoint $O/best.pth --tensors $D/$SET/tensors --kinds SNV \
    --chroms chr1 --output $O/test_chr1 --cache-dir $R/index_cache --num-workers 6 --batch-size 256 $EXTRA
fi
[ -s $O/vcf_chr1/eval_nobed/report.json ] || SET_NAME=$SET bash /scratch/jshen/Github/Pansoma/machine_learning/pansoma_net_v2/jobs/vcf.sh $OUTRUN SNV $O/vcf_chr1
echo "task done $(date)"

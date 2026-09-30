#!/bin/bash
# NOT SUBMITTED (no Slurm from this session). Exact v1 truth-level curve on HG008 Illumina chr1 (all 361,394 v1 testing
# tensors, incl. 16,018 off-reference) for the published model (e067) and, as a cross-check, e053.
# CPU measured here: 6.6-9 tensors/s on 4 threads (~12 h per checkpoint). GPU: expect well under 30 min per checkpoint;
# the job reads ~36 GB of int8 tensors per checkpoint. Memory: model 2.4 GB checkpoint + one mmap'd shard page cache.
#SBATCH --job-name=v1_chr1_full
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
S=/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head
# copy $S to a persistent place first if this session's scratch may be cleaned; outputs go to $S/v1_full
PY=/wanglab/jshen/anaconda3/bin/python
cd /scratch/jshen/Github/Pansoma/machine_learning/pansoma_net
$PY -u "$S/v1_full_chr1_gpu_test.py" --out "$S/v1_full" --tag e067 --batch 256
$PY -u "$S/v1_full_chr1_gpu_test.py" --out "$S/v1_full" --tag e053 --batch 256 \
  --ckpt /scratch/jshen/Github/GoogleNet/HG008T_GIAB_AF-HPRC_CE_Large_Model_V2_weight200/model_e053_f1_0.1733.pth
# score (CPU, ~1 min): exact v1 rows next to the v2 rows, no PoN and with the two PoN rules
$PY "$S/score_v1_full.py" "$S/v1_full/v1_chr1_test_p_true.e067.npy" "$S/v1_full/v1_chr1_test_p_true.e053.npy"

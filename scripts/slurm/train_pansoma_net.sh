#!/usr/bin/env bash
#SBATCH --job-name=pansoma_train
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:h100:3
#SBATCH --mem=300G
#SBATCH --time=72:00:00
#SBATCH --output=/scratch/jshen/Log_info/pansoma_train.%j.out
#SBATCH --error=/scratch/jshen/Log_info/pansoma_train.%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

: "${PANSOMA_NET_DIR:=${PROJECT_ROOT}/machine_learning/pansoma_net}"
: "${TRAIN_DATA_PATHS_FILE:=train_data_dir.txt}"
: "${VAL_DATA_PATHS_FILE:=val_data_dir.txt}"
: "${OUT_DIR:=./Pansoma_Model}"
: "${BATCH_SIZE:=256}"
: "${LR:=0.0001}"
: "${EPOCHS:=80}"
: "${NPROC_PER_NODE:=3}"
: "${CUDA_VISIBLE_DEVICES:=0,1,2}"

mkdir -p /scratch/jshen/Log_info

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export CUDA_VISIBLE_DEVICES

echo "Host: $(hostname)"
echo "Start: $(date)"
nvidia-smi || true

cd "${PANSOMA_NET_DIR}"

torchrun --nproc_per_node="${NPROC_PER_NODE}" scripts/train_6channels_npy_pansoma_DDP.py \
  --ddp \
  --train_data_paths_file "${TRAIN_DATA_PATHS_FILE}" \
  --val_data_paths_file "${VAL_DATA_PATHS_FILE}" \
  --batch_size "${BATCH_SIZE}" \
  --lr "${LR}" \
  -o "${OUT_DIR}" \
  --training_data_ratio 0.3 \
  --epochs "${EPOCHS}"

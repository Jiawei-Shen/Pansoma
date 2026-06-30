#!/usr/bin/env bash
#SBATCH --job-name=pansoma_infer
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/jshen/Log_info/pansoma_infer.%j.out
#SBATCH --error=/scratch/jshen/Log_info/pansoma_infer.%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

: "${PANSOMA_NET_DIR:=${PROJECT_ROOT}/machine_learning/pansoma_net}"
: "${INPUT_DIR:?Set INPUT_DIR to tensor shard directory}"
: "${CKPT:?Set CKPT to model checkpoint}"
: "${OUT_PREFIX:?Set OUT_PREFIX to output prefix}"
: "${MAP_JSON:?Set MAP_JSON to node mapping JSON}"
: "${VARIANT_SUMMARY:?Set VARIANT_SUMMARY to variant_summary.ndjson}"
: "${BATCH_SIZE:=32}"
: "${NUM_WORKERS:=4}"

mkdir -p "$(dirname "${OUT_PREFIX}")" /scratch/jshen/Log_info

cd "${PANSOMA_NET_DIR}"

python scripts/test_5channels_npy_pansoma.py \
  --input_dir "${INPUT_DIR}" \
  --ckpt "${CKPT}" \
  --out_prefix "${OUT_PREFIX}" \
  --input_mode shard \
  --num_workers "${NUM_WORKERS}" \
  --map_json "${MAP_JSON}" \
  --variant_summary "${VARIANT_SUMMARY}" \
  --batch_size "${BATCH_SIZE}" \
  --gpu-normalize

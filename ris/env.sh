# Sourced by the ris/ job scripts inside the pytorch/pytorch:2.8.0-cuda12.6-cudnn9-devel container.
# The image has torch 2.8 (CUDA 12.6), numpy and zstandard; timm and pysam live in a user base on storage2
# (setup_env.sh), so the stock image is enough.
export REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PANSOMA_RUNS=${PANSOMA_RUNS:-/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs}
export PANSOMA_DATA=${PANSOMA_DATA:-/storage2/fs1/epigenome/Active/jiawei.s/Dataset}
export PYTHONUSERBASE=$PANSOMA_RUNS/env/pyuser
export PYTHONNOUSERSITE=
export PATH=$PYTHONUSERBASE/bin:/opt/conda/bin:$PATH
export RTG=$PANSOMA_RUNS/env/rtg-tools-3.13/rtg
# Triton / inductor caches on storage2, shared by the jobs
export TRITON_CACHE_DIR=$PANSOMA_RUNS/env/cache/triton
export TORCHINDUCTOR_CACHE_DIR=$PANSOMA_RUNS/env/cache/inductor
mkdir -p "$TRITON_CACHE_DIR" "$TORCHINDUCTOR_CACHE_DIR"

#!/bin/bash
# The GPU tests of pansoma_net_v2 (GPU encoding = CPU, a bf16 step learns, compile + channels_last = eager).
#   ris/submit.sh pnv2_gpu_check gpu_check.sh
set -euo pipefail
source "$(dirname "$0")/env.sh"
echo "host $(hostname) CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"
cd "$REPO/machine_learning"
python -c "import torch; print('torch', torch.__version__, torch.cuda.device_count(), torch.cuda.get_device_name(0),
'bf16', torch.cuda.is_bf16_supported())"
python -m unittest -v pansoma_net_v2.tests.test_gpu 2>&1 | tail -12

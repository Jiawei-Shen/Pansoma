#!/bin/bash
# Once: timm and pysam into $PYTHONUSERBASE, then the CPU unit tests. A CPU job:
#   ris/submit.sh -g 0 -n 4 -m 16 pnv2_setup setup_env.sh
set -euo pipefail
source "$(dirname "$0")/env.sh"
hostname; which python; python -V
python -m pip install --user --no-cache-dir "timm==1.0.*" zstandard pysam
python -c "import torch, timm, zstandard, pysam, numpy; print('torch', torch.__version__, 'timm', timm.__version__,
'zstandard', zstandard.__version__, 'pysam', pysam.__version__, 'numpy', numpy.__version__)"
cd "$REPO/machine_learning"
PANSOMA_TEST_TENSORS=/nonexistent python -m unittest discover -s pansoma_net_v2/tests -t . 2>&1 | tail -5

#!/bin/bash
# After set_call_threshold.sh: the changed checkpoint against the unchanged copy. Every tensor of model_state_dict must
# be equal, and every other entry but val.somatic_at_recall_postpon_callable[call_recall] and call_threshold_manual.
# Usage: check_checkpoint.sh RUN ORIGINAL_RUN. A CPU job:
#   ris/submit.sh -g 0 -n 2 -m 16 pnv2_chkthr_RUN ../analysis/colo829t_dilution_newpipe_20261112/check_checkpoint.sh \
#       COLO829_Illumina_SNV_joint_newpipe_e5 COLO829_Illumina_SNV_joint_newpipe
set -euo pipefail
source "$(dirname "$0")/../../ris/env.sh"
python - "$PANSOMA_RUNS/$1/best.pth" "$PANSOMA_RUNS/$2/best.pth" <<'PY'
import sys
import torch


def same(a, b):
    if isinstance(a, torch.Tensor) or isinstance(b, torch.Tensor):
        return isinstance(a, torch.Tensor) and isinstance(b, torch.Tensor) and a.dtype == b.dtype and torch.equal(a, b)
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    return a == b or (a != a and b != b)


new, old = (torch.load(f, map_location="cpu", weights_only=False) for f in sys.argv[1:3])
key = f"{old['call_recall']:g}"
manual = new.pop("call_threshold_manual")
record = new["val"]["somatic_at_recall_postpon_callable"][key]
assert same(record["validation"], old["val"]["somatic_at_recall_postpon_callable"][key]), "validation record"
new["val"]["somatic_at_recall_postpon_callable"][key] = record["validation"]
differ = [k for k in sorted(set(new) | set(old)) if k not in new or k not in old or not same(new[k], old[k])]
n = len(old["model_state_dict"])
print(f"{sys.argv[1]}: {n} weight tensors equal: {same(new['model_state_dict'], old['model_state_dict'])}; "
      f"other entries that differ: {differ or 'none'}; threshold {record['validation']['threshold']} -> "
      f"{record['threshold']} ({manual['note']})")
assert not differ
PY

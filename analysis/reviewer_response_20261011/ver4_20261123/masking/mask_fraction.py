"""Share of the tensor cells that --mask-off-node zeroes, chr1 of the nine test sets of Response Figure 1.2 (ver4).

Per set: an evenly strided sample of N chr1 tensors (all labels) and an evenly strided sample of up to N chr1 tensors
labelled somatic (index label 1). For each sample, with the mask of the patched tree (offnode.own_mask, the function
the masked predictions ran) and the side file the predictions read:
  tensors_with_masked  share of the tensors with at least one non-empty cell outside the candidate's own node visit
  cells_nonempty       cells with any of the 8 channels non-zero
  masked_share         share of those cells the mask sets to 0
  covered_masked_share the same over the cells with a read base or gap (channel 4 > 0)
Read-only. usage: python mask_fraction.py [N=3000]  -> mask_fraction.tsv / mask_fraction.json next to this script
Run on a compute node (tmp/ver4_20261123/masking/mask_fraction.sh)."""
import json
import sys
from pathlib import Path

import numpy as np

W = Path("/scratch/jshen/Github/Pansoma/tmp/ver4_20261123/masking")
sys.path.insert(0, str(W / "code_ml3_offnode/machine_learning"))
from pansoma_net_v2 import offnode  # noqa: E402
from pansoma_net_v2.data import TensorDataset, load_parts  # noqa: E402

D = Path("/scratch/jshen/data/pansoma_v2_tensors")
SIDE = Path("/scratch/jshen/data/reviewer_experiments_20261013/ver4_20261123/masking/side")
CACHE = "/scratch/jshen/data/pansoma_net_v2_runs/ver4_masking_20261123/index_cache"
HERE = Path(__file__).resolve().parent
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
sets = [l.split("\t")[0] for l in (W / "sets.tsv").read_text().splitlines() if l.strip()]
rows = []
for s in sets:
    (index, positions), = load_parts([str(D / s / "tensors")], ["SNV"], CACHE, chroms=["chr1"], labelled=False)
    positions = np.asarray(positions, np.int64)
    with np.load(SIDE / f"{s}.SNV.side.npz") as z:
        assert len(z["node"]) == len(index) and (z["node"] == index.arrays["node"]).all(), s
        side = {k: z[k] for k in ("start", "slots", "span", "node_len")}
    somatic = positions[index.arrays["label"][positions] == 1]
    for group, pos in (("all", positions), ("somatic", somatic)):
        pick = pos[np.unique(np.linspace(0, len(pos) - 1, min(N, len(pos))).astype(np.int64))]
        ds = TensorDataset([(index, pick)])
        t = tm = cells = masked = cov = covm = 0
        for i in range(len(ds)):
            x = ds[i][0].numpy()
            p = pick[i]
            own = offnode.own_mask(x[None], [side["start"][p]], [side["slots"][p]], [side["span"][p]], [side["node_len"][p]])[0]
            nonempty = (x != 0).any(0)
            covered = x[4] > 0
            m = int((nonempty & ~own).sum())
            t += 1; tm += m > 0; cells += int(nonempty.sum()); masked += m
            cov += int(covered.sum()); covm += int((covered & ~own).sum())
        rows.append(dict(set=s, group=group, chr1_tensors=int(len(pos)), sampled=t, tensors_with_masked=round(tm / t, 4),
                         cells_nonempty=cells, masked_cells=masked, masked_share=round(masked / max(cells, 1), 4),
                         covered_masked_share=round(covm / max(cov, 1), 4)))
        print(json.dumps(rows[-1]), flush=True)
keys = list(rows[0])
(HERE / "mask_fraction.tsv").write_text("\n".join(["\t".join(keys)] + ["\t".join(str(r[k]) for k in keys) for r in rows]) + "\n")
(HERE / "mask_fraction.json").write_text(json.dumps(rows, indent=1) + "\n")

"""Pick chr1 target nodes of a real set for a small --debug-rows build: every node holding a somatic (label 1) SNV or
INDEL tensor plus an evenly strided sample of other chr1 tensor nodes. usage: pick_nodes.py SET N_OTHER OUT"""
import json, sys
from pathlib import Path
import numpy as np
CACHE = Path("/scratch/jshen/data/pansoma_net_v2_runs/index_cache")
s, n_other, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
nodes = set()
for kind in ("SNV", "INDEL"):
    f = sorted(CACHE.glob(f"{s}.tensors.{kind}.*.npz"), key=lambda p: p.stat().st_mtime)[-1]
    z = np.load(f)
    meta = json.loads(str(z["meta"]))
    c1 = meta["chroms"].index("chr1")
    m = z["chrom"] == c1
    node, lab = z["node"][m], z["label"][m]
    som = set(node[lab == 1].tolist())
    other = np.unique(node[lab != 1])
    pick = other[np.unique(np.linspace(0, len(other) - 1, n_other).astype(int))]
    nodes |= som | set(pick.tolist())
    print(kind, f.name, "somatic nodes", len(som), "other picked", len(pick))
Path(out).write_text("\n".join(map(str, sorted(nodes))) + "\n")
print("nodes", len(nodes))

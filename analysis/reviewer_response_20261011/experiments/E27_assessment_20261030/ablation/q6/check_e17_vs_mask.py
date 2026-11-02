"""Per-tensor: E1.7 geometric rule (start < 50 or L - end < 50) vs own_node_mask having masked covered cells; strided
sample of chr1 shard 0. usage: python check_e17_vs_mask.py N SET KIND"""
import json, sqlite3, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from node_mask import own_node_mask
D = Path("/scratch/jshen/data/pansoma_v2_tensors")
n_sample, s, kind = int(sys.argv[1]), sys.argv[2], sys.argv[3]
kd = D / s / "tensors" / kind
gi = json.load(open(kd / "manifest.json"))["graph_index"]["path"]
recs = {}
with open(kd / "chr1_variant_summary.ndjson") as f:
    for line in f:
        r = json.loads(line)
        if r["shard_file"] != "chr1_shard_00000_data.npy":
            break
        recs[r["index_within_shard"]] = (r["node_id"], r["start"], r["end"], r["site_layout"]["insertion_slots"], r["site_layout"]["span"])
n = len(recs)
idx = np.unique(np.linspace(0, n - 1, n_sample).astype(int)).tolist()
con = sqlite3.connect(f"file:{gi}?mode=ro", uri=True)
nodes = sorted({recs[i][0] for i in idx})
lens = dict(con.execute(f"select node_id, length(seq) from nodes where node_id in ({','.join('?' * len(nodes))})", nodes).fetchall())
# whole-shard E1.7 share needs every node length
allnodes = sorted({v[0] for v in recs.values()})
for k in range(0, len(allnodes), 50000):
    ch = allnodes[k:k + 50000]
    lens.update(con.execute(f"select node_id, length(seq) from nodes where node_id in ({','.join('?' * len(ch))})", ch).fetchall())
geo_all = np.mean([st < 50 or lens[nd] - en < 50 for nd, st, en, _, _ in recs.values()])
x_all = np.load(kd / "chr1_shard_00000_data.npy", mmap_mode="r")
tab = {}
for i in idx:
    nd, st, en, slots, span = recs[i]
    x = np.asarray(x_all[i]); L = lens[nd]
    geo = st < 50 or L - en < 50
    m = int(((x[4] > 0) & ~own_node_mask(x, st, L, slots)).sum())
    tab[(geo, m > 0)] = tab.get((geo, m > 0), 0) + 1
print(f"{s} {kind}: shard-0 tensors {n}, E1.7 rule share over all shard-0 tensors {geo_all:.4f}; sample {len(idx)} (geo, masked>0): {tab}")

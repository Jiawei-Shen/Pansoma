"""Compare the own-node rule (q6/node_mask.own_node_mask, plus the site span kept) with the per-column node
coordinates of a --debug-rows build. Node lengths are read from the graph-index sqlite for the build's nodes only.
usage: check_mask.py BUILD_KIND_DIR GRAPH_SQLITE"""
import json, sqlite3, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6")
from node_mask import own_node_mask
d, gi = Path(sys.argv[1]), sys.argv[2]
recs = [json.loads(l) for l in open(d / "variant_summary.ndjson")]
ids = sorted({r["node_id"] for r in recs})
con = sqlite3.connect(f"file:{gi}?mode=ro", uri=True)
lens = {}
for k in range(0, len(ids), 900):
    ch = ids[k:k + 900]
    lens.update(con.execute(f"select node_id, length(seq) from nodes where node_id in ({','.join('?' * len(ch))})", ch).fetchall())
shards = {}
n_t = n_cells = n_bad = n_off = n_span_other = n_t_off = 0
bad = []
orient = {"+": 0, "-": 0}
for r in recs:
    f = d / f"shard_{r['shard_index']:05d}_data.npy"
    if f not in shards:
        shards[f] = np.load(f)
    x = shards[f][r["index_within_shard"]]
    lay = r["site_layout"]
    own = own_node_mask(x, r["start"], lens[r["node_id"]], lay["insertion_slots"])
    a = r["anchor_column"]
    assert a == x.shape[2] // 2
    site = set(range(a, a + lay["insertion_slots"] + lay["span"]))
    n_t += 1
    orient[r.get("orientation", "+")] = orient.get(r.get("orientation", "+"), 0) + 1
    t_off = 0
    for ri, row in enumerate(r["rows"]):
        for ci, c in enumerate(row["columns"]):
            if c is None:
                continue
            truth = c["node_id"] == r["node_id"] and c["mapping_index"] == row["anchor_mapping_index"]
            if ci in site and not truth:
                n_span_other += 1
            n_cells += 1
            n_off += not truth
            t_off += not truth
            if truth != own[ri, ci]:
                n_bad += 1
                if len(bad) < 5:
                    bad.append((r["candidate_id"], ri, ci, c, bool(own[ri, ci])))
    n_t_off += t_off > 0
print(f"{d}: tensors {n_t} (orientation {orient}), with off-node cells {n_t_off}; covered cells {n_cells}, off own node visit "
      f"{n_off} ({n_off / max(n_cells, 1):.3f}); rule disagreements {n_bad}; site-span cells on a further node {n_span_other}")
for b in bad:
    print("  disagreement", b)

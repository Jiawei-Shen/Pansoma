"""Categorize own-node rule disagreements in a --debug-rows build. usage: categorize.py KIND_DIR GRAPH_SQLITE"""
import json, sqlite3, sys
from collections import Counter
from pathlib import Path
import numpy as np
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/ablation/q6/code")
from pansoma_net_v2.offnode import own_mask
d, gi = Path(sys.argv[1]), sys.argv[2]
recs = [json.loads(l) for l in open(d / "variant_summary.ndjson")]
ids = sorted({r["node_id"] for r in recs})
con = sqlite3.connect(f"file:{gi}?mode=ro", uri=True)
lens = {}
for k in range(0, len(ids), 900):
    ch = ids[k:k + 900]
    lens.update(con.execute(f"select node_id, length(seq) from nodes where node_id in ({','.join('?' * len(ch))})", ch).fetchall())
shards, cat, n_cells, n_off = {}, Counter(), 0, 0
tens_bad = 0
for r in recs:
    f = d / f"shard_{r['shard_index']:05d}_data.npy"
    if f not in shards:
        shards[f] = np.load(f)
    x = shards[f][r["index_within_shard"]]
    lay = r["site_layout"]
    own = own_mask(x[None], [r["start"]], [lay["insertion_slots"]], [lay["span"]], [lens[r["node_id"]]])[0]
    a = r["anchor_column"]; site_end = a + lay["insertion_slots"] + lay["span"]
    bad_here = 0
    for ri, row in enumerate(r["rows"]):
        for ci, c in enumerate(row["columns"]):
            if c is None:
                continue
            truth = c["node_id"] == r["node_id"] and c["mapping_index"] == row["anchor_mapping_index"]
            n_cells += 1; n_off += not truth
            if truth != bool(own[ri, ci]):
                bad_here += 1
                where = "left" if ci < a else ("site" if ci < site_end else "right")
                kind = "boundary" if "boundary" in c else ("offset" if "offset" in c else "/".join(sorted(c)))
                same_node = c["node_id"] == r["node_id"]
                cat[(where, kind, f"g5={int(x[5, ri, ci])}", f"op={int(x[4, ri, ci])}", "truth_own" if truth else "truth_off",
                     "same_node" if same_node else "other_node")] += 1
    tens_bad += bad_here > 0
print(f"{d}: covered {n_cells}, off {n_off}, disagreements {sum(cat.values())} in {tens_bad} of {len(recs)} tensors")
for k, v in cat.most_common(15):
    print(" ", v, k)

"""Per-tensor off-node mask parameters of a merged kind directory, in KindIndex order (manifest chromosome order, summary
line order): node_id, start (site offset on the node, candidate-node forward orientation), insertion_slots, span, node
length (graph-index sqlite). Checked against the index cache's node array when one is given.
usage: build_side.py KIND_DIR OUT.npz [INDEX_CACHE.npz]"""
import json, re, sqlite3, sys, time
from pathlib import Path
import numpy as np

kd, out = Path(sys.argv[1]), Path(sys.argv[2])
cache = Path(sys.argv[3]) if len(sys.argv) > 3 else None
m = json.loads((kd / "manifest.json").read_text())
gi = m["graph_index"]["path"]
R_NODE = re.compile(rb'"node_id": (\d+)')
R_START = re.compile(rb'"start": (-?\d+)')
R_SLOTS = re.compile(rb'"insertion_slots": (\d+)')
R_SPAN = re.compile(rb'"span": (\d+)')
R_ANCH = re.compile(rb'"anchor_column": (\d+)')
node, start, slots, span = [], [], [], []
t0 = time.time()
checked = 0
for chrom in m["chromosomes"]:
    info = m["chromosomes"][chrom]
    n0 = len(node)
    with open(kd / info["summary"], "rb") as f:
        for line in f:
            i_node = R_NODE.search(line)
            node.append(int(i_node.group(1)))
            start.append(int(R_START.search(line, i_node.end()).group(1)))
            slots.append(int(R_SLOTS.search(line).group(1)))
            span.append(int(R_SPAN.search(line).group(1)))
            if int(R_ANCH.search(line).group(1)) != 50:
                raise SystemExit(f"anchor_column != 50 in {chrom}")
            if checked < 200 and len(node) % 997 == 1:  # spot-check the regexes against a full parse
                r = json.loads(line)
                assert (r["node_id"], r["start"], r["site_layout"]["insertion_slots"], r["site_layout"]["span"]) == \
                    (node[-1], start[-1], slots[-1], span[-1]), (chrom, r["candidate_id"])
                checked += 1
    if len(node) - n0 != info["tensors"]:
        raise SystemExit(f"{chrom}: {len(node) - n0} summary lines, manifest {info['tensors']}")
print(f"summaries: {len(node):,} tensors, {checked} spot checks, {time.time() - t0:.0f} s", flush=True)
node = np.asarray(node, np.int64)
ids = np.unique(node)
con = sqlite3.connect(f"file:{gi}?mode=ro", uri=True)
lens = {}
for k in range(0, len(ids), 900):
    ch = ids[k:k + 900].tolist()
    lens.update(con.execute(f"select node_id, length(seq) from nodes where node_id in ({','.join('?' * len(ch))})", ch).fetchall())
missing = [i for i in ids.tolist() if i not in lens]
if missing:
    raise SystemExit(f"{len(missing)} nodes without a length, e.g. {missing[:5]}")
node_len = np.array([lens[i] for i in node.tolist()], np.int32)
print(f"node lengths: {len(ids):,} nodes, {time.time() - t0:.0f} s", flush=True)
if cache is not None:
    z = np.load(cache)
    zn = z["node"]
    assert len(zn) == len(node) and (zn == node).all(), "index cache node order differs"
    print("index cache node array: identical", flush=True)
start = np.asarray(start, np.int32); slots = np.asarray(slots, np.int16); span = np.asarray(span, np.int16)
assert (start >= 0).all() and (start < node_len).all()
np.savez(out, node=node, start=start, slots=slots, span=span, node_len=node_len,
         meta=np.array(json.dumps(dict(kind_dir=str(kd.resolve()), graph_index=gi, chroms=list(m["chromosomes"]),
                                       tensors=int(len(node))))))
print(f"wrote {out}; start==0: {(start == 0).mean():.3f}, start+span>node_len: {(start + span > node_len).mean():.4f}")

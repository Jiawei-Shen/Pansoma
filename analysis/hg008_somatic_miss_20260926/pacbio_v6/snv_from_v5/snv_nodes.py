"""Discovery statistics of the GRCh38 node holding each SNV miss (and SNV controls): why (not) a target node?"""
import json, mmap, re, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/scratch/jshen/Github/Pansoma")
from indexed_gam_pipeline_v2.tensor_postprocessing.reference_path import ReferencePath

HERE = Path(__file__).resolve().parent
RUN = Path("/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125")
path = ReferencePath("/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path")
misses = [m for m in json.load(open(HERE / "nearby.json")) if m["kind"] == "SNP"]
nodes = {int(k.split(":")[0]) for m in misses for k in m["keys"]}
record = re.compile(rb'"(\d+)":\s*\{\s*"perfect":\s*(\d+),\s*"not_perfect":\s*(\d+),\s*"max_read_length":\s*(\d+)')
stats = {}
with open(RUN / "discovery/node_stats.json", "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as data:
    for m in record.finditer(data):
        n = int(m.group(1))
        if n in nodes:
            stats[n] = dict(perfect=int(m.group(2)), not_perfect=int(m.group(3)))
targets = np.fromfile(RUN / "discovery/target_nodes.txt", sep=" ", dtype=np.int64)
out = {}
for m in misses:
    (key,) = m["keys"] or [None]
    if key is None:
        continue
    n = int(key.split(":")[0])
    s = stats.get(n, dict(perfect=0, not_perfect=0))
    total = s["perfect"] + s["not_perfect"]
    at = min(int(np.searchsorted(targets, n)), len(targets) - 1)
    out[m["truth_id"]] = dict(node=n, node_length=int(path.lengths[n]), **s,
                              imperfect_fraction=s["not_perfect"] / total if total else None,
                              target=bool(targets[at] == n))
json.dump(out, open(HERE / "snv_nodes.json", "w"))
print(len(out), "SNV misses with node stats")

"""GRCh38 node of each SNV truth locus analysed (misses + controls): node length, discovery stats, target or not."""
import json, mmap, re, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reads as R  # noqa: E402

HERE = Path(__file__).resolve().parent
path = R.ReferencePath(R.REFPATH)
loci = [m for m in R.loci() if m["kind"] == "SNP"]
nodes = {int(k.split(":")[0]) for m in loci for k in m["keys"]}
record = re.compile(rb'"(\d+)":\s*\{\s*"perfect":\s*(\d+),\s*"not_perfect":\s*(\d+),\s*"max_read_length":\s*(\d+)')
stats = {}
with open(R.RUN / "discovery/node_stats.json", "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as data:
    for m in record.finditer(data):
        n = int(m.group(1))
        if n in nodes:
            stats[n] = dict(perfect=int(m.group(2)), not_perfect=int(m.group(3)))
targets = set()
for f in sorted((R.RUN / "v3_run/parts").glob("nodes_*.txt")):
    targets.update(np.fromfile(f, sep=" ", dtype=np.int64).tolist())
out = {}
for m in loci:
    if len(m["keys"]) != 1:
        continue
    n = int(m["keys"][0].split(":")[0])
    s = stats.get(n, dict(perfect=0, not_perfect=0))
    out[m["truth_id"]] = dict(node=n, node_length=int(path.lengths[n]), **s, target=n in targets)
json.dump(out, open(HERE / "snv_nodes.json", "w"))
print(len(out), "SNV loci with node stats")

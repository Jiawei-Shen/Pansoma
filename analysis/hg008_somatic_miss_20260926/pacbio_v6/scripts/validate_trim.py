"""Trimmed decoding (TRIM_BP=600) must reproduce the whole-record results of reads_part*.jsonl exactly.

Re-analyses N already finished loci (whole-record decode) with trimming and compares every field
(tensor labels excluded: they are a lookup, not decoding).
"""
import json, os, random, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
os.environ["TRIM_BP"] = "600"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import reads as R

def strip(r):
    r = dict(r)
    r["alt_candidates"] = {c: d["reads"] for c, d in r["alt_candidates"].items()}
    return r

done = [json.loads(l) for f in sorted(Path(".").glob("reads_part*.jsonl")) for l in open(f)]
rng = random.Random(3)
sample = rng.sample(done, min(int(sys.argv[1]) if len(sys.argv) > 1 else 60, len(done)))
ids = {r["truth_id"] for r in sample}
work = sorted([m for m in R.loci() if m["truth_id"] in ids], key=lambda m: (m["chrom"], m["pos0"]))
targets = R.NodeSet(np.concatenate([np.fromfile(f, sep=" ", dtype=np.int64) for f in sorted((R.RUN / "v6_run/parts").glob("nodes_*.txt")) + sorted((R.RUN / "v6_run/parts").glob("supplement_0[0-9]/nodes_*.txt"))]))
with ProcessPoolExecutor(16, initializer=R.init, initargs=(targets, {})) as pool:
    trimmed = {r["truth_id"]: r for part in pool.map(R.analyse, [[m] for m in work]) for r in part}
diff = 0
for r in sample:
    a, b = strip(r), strip(trimmed[r["truth_id"]])
    if a != b:
        diff += 1
        print("DIFF", r["chrom"], r["vcf_pos"], {k: (a[k], b[k]) for k in a if a[k] != b.get(k)})
print(f"{len(sample)} loci compared, {diff} differ")

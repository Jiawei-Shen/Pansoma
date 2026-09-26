"""Trimmed decoding (TRIM_BP=600) must reproduce whole-record decoding (TRIM_BP=0) on sampled ONT loci.

Both runs use analyse() of reads.py; every output field is compared (tensor labels excluded: a lookup).
"""
import json, os, random, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import reads as R


def init(targets, trim):
    R.init(targets, {})
    R.TRIM_BP = trim


def strip(r):
    r = dict(r)
    r["alt_candidates"] = {c: d["reads"] for c, d in r["alt_candidates"].items()}
    return r


rng = random.Random(3)
loci = R.loci()
work = sorted(rng.sample(loci, int(sys.argv[1]) if len(sys.argv) > 1 else 24), key=lambda m: (m["chrom"], m["pos0"]))
targets = R.NodeSet(np.concatenate([np.fromfile(f, sep=" ", dtype=np.int64) for f in sorted((R.RUN / "v3_run/parts").glob("nodes_*.txt"))]))
out = {}
for trim in (600, 0):
    with ProcessPoolExecutor(int(os.environ.get("SLURM_CPUS_PER_TASK", 16)), initializer=init, initargs=(targets, trim)) as pool:
        out[trim] = {r["truth_id"]: r for part in pool.map(R.analyse, [[m] for m in work]) for r in part}
    print("done trim", trim, flush=True)
diff = 0
for tid in out[0]:
    a, b = strip(out[0][tid]), strip(out[600][tid])
    if a != b:
        diff += 1
        print("DIFF", a["chrom"], a["vcf_pos"], {k: (a[k], b[k]) for k in a if a[k] != b.get(k)})
print(f"{len(out[0])} loci compared, {diff} differ")

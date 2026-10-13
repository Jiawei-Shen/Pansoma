"""Check a copied graph index before graph_vcf / linear_vcf use it (a copy can be truncated; scp writes the final names).
Usage (in the container: ris/submit.sh -g 0 -n 2 -m 16 NAME py.sh tools/check_graph_index.py SQLITE REF_PATH_DIR SET_TENSORS):
- SQLITE: `nodes` row count equals the set's SNV manifest graph_index.nodes, the table and graph_metadata read;
- REF_PATH_DIR (tools.graph_prep ref-path-scan): meta.json reads, every .npy loads (header + size consistent),
  path_nodes / path_reverse lengths equal the sum of the contigs' nodes, walks.ndjson has meta.json's walk count.
Exit 0 when everything matches."""
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np

sqlite_path, ref_dir, tensors = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
problems = []
want_nodes = json.load(open(tensors / "SNV" / "manifest.json"))["graph_index"]["nodes"]
c = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
n = c.execute("select count(*) from nodes").fetchone()[0]
cur = c.execute("select * from graph_metadata")
rows, cols = cur.fetchall(), [d[0] for d in cur.description]
print(f"{sqlite_path}: {n:,} nodes (manifest {want_nodes:,}); graph_metadata {len(rows)} rows, columns {cols[:8]}")
if n != want_nodes:
    problems.append(f"sqlite nodes {n} != manifest {want_nodes}")
m = json.load(open(ref_dir / "meta.json"))
contig_nodes = sum(x["nodes"] for x in m["contigs"])
for f in sorted(ref_dir.glob("*.npy")):
    a = np.load(f, mmap_mode="r")
    print(f"{f.name}: shape {a.shape} {a.dtype}")
    if f.stem in ("path_nodes", "path_reverse") and a.shape[0] != contig_nodes:
        problems.append(f"{f.name} length {a.shape[0]} != contig nodes {contig_nodes}")
walks = sum(1 for _ in open(ref_dir / "walks.ndjson"))
print(f"walks.ndjson: {walks} lines (meta {m.get('walks')}); contigs {len(m['contigs'])}, contig nodes {contig_nodes:,}")
if m.get("walks") is not None and walks != m["walks"]:
    problems.append(f"walks {walks} != meta {m['walks']}")
print("OK" if not problems else "PROBLEMS: " + "; ".join(problems))
sys.exit(1 if problems else 0)

"""E9b check: the CHM13 coordinates of prep_graph.py (vg paths -A -S CHM13 on the d9 GBZ) against the d9 GRCh38_CHM13
path_index sqlite (its CHM13 rows are complete; node ids TEXT), on the nodes of the projected calls of the given TSVs
plus random nodes.

    python check_chm13_sqlite.py <run>.projected.tsv [...] [--random 300] [--output check.json]
"""
import argparse
import json
import random
import sqlite3
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath  # noqa: E402

SQLITE = "/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13.path_index.sqlite"
CHM13 = "/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection/graph/d9/chm13_path"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv", nargs="*")
    ap.add_argument("--random", type=int, default=300)
    ap.add_argument("--max-call-nodes", type=int, default=600)
    ap.add_argument("--output")
    args = ap.parse_args()
    rp = ReferencePath(CHM13)
    nodes = set()
    for t in args.tsv:
        with open(t) as f:
            cols = f.readline().rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(cols, line.rstrip("\n").split("\t")))
                if r["filter"] == "PASS":
                    nodes.add(int(r["node"]))
    rng = random.Random(20261013)
    call_nodes = sorted(nodes)
    if len(call_nodes) > args.max_call_nodes:
        call_nodes = sorted(rng.sample(call_nodes, args.max_call_nodes))
    present = np.flatnonzero(np.asarray(rp.lengths) > 0)
    rand_nodes = sorted(int(x) for x in rng.sample(list(present[::97]), args.random))
    db = sqlite3.connect(f"file:{SQLITE}?mode=ro", uri=True)
    result = {}
    for label, group in (("call_nodes", call_nodes), ("random_nodes", rand_nodes)):
        n = dict(nodes=len(group), both_on_chm13=0, neither=0, only_sqlite=0, only_prep=0, start_equal=0,
                 start_differs=0, multi_visit=0, examples=[])
        for node in group:
            rows = db.execute("SELECT chrom, start0, end0, orientation FROM path_coords WHERE assembly = 'CHM13' "
                              "AND node_id = ?", (str(node),)).fetchall()
            on_prep = int(rp.visits[node]) >= 1
            if rows and on_prep:
                n["both_on_chm13"] += 1
                if int(rp.visits[node]) > 1 or len(rows) > 1:
                    n["multi_visit"] += 1
                chrom = rp.contigs[int(rp.chrom[node])]
                starts = {(c, int(s)) for c, s, _, _ in rows}
                if (chrom, int(rp.start0[node])) in starts:
                    n["start_equal"] += 1
                else:
                    n["start_differs"] += 1
                    if len(n["examples"]) < 10:
                        n["examples"].append(dict(node=node, prep=[chrom, int(rp.start0[node]), int(rp.lengths[node])],
                                                  sqlite=[list(r) for r in rows[:3]]))
            elif rows:
                n["only_sqlite"] += 1
            elif on_prep:
                n["only_prep"] += 1
            else:
                n["neither"] += 1
        result[label] = n
    text = json.dumps(result, indent=2)
    print(text)
    if args.output:
        Path(args.output).write_text(text + "\n")


if __name__ == "__main__":
    main()

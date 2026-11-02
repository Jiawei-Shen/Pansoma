"""Own-node column mask of v4 tensors, reconstructed from the tensor alone (+ the summary's start / site_layout and the
node length), and its check against the per-cell graph coordinates of a --debug-rows build.

Rule (indexed_gam_pipeline_v4/candidates.py site_window / decode_alignment / left_align_indels):
  a cell consumes one graph base  <=>  channel 5 (graph base) in 1..5 (A C G T N); I, G (no-insertion slot) and the
  boundary part of a complex edit have graph base '-' (6). Rows are in forward-node orientation of the candidate node.
  left flank  (k < anchor):            own  <=>  #consuming cells in [k, anchor) <= start
  insertion slots [anchor, anchor+s):  own (the site's own boundary)
  right       (k >= anchor + s):       own  <=>  start + #(consuming or empty) cells in [anchor+s, k) < node length
  (empty cells on the right are counted: coverage starting inside the span is padded with empty cells at its true
  offset; empty cells after a read end are empty anyway).
usage: python node_mask.py <typed build dir with --debug-rows (SNV or INDEL)> <graph.sqlite>
"""
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np


def own_node_mask(x, start, node_len, slots):
    """x: (8, R, W) int8 -> (R, W) bool, True where the cell lies on the candidate's own node visit."""
    _, R, W = x.shape
    a = W // 2
    covered = x[4] > 0
    consuming = (x[5] >= 1) & (x[5] <= 5)
    own = np.zeros((R, W), bool)
    left = consuming[:, :a][:, ::-1].cumsum(1)[:, ::-1]          # #consuming in [k, a)
    own[:, :a] = left <= start
    own[:, a:a + slots] = True
    r0 = a + slots
    cons_r = consuming[:, r0:] | ~covered[:, r0:]
    before = np.concatenate([np.zeros((R, 1), int), cons_r.cumsum(1)[:, :-1]], 1)  # #in [r0, k)
    own[:, r0:] = start + before < node_len
    return own


def main():
    d, gi = Path(sys.argv[1]), sys.argv[2]
    con = sqlite3.connect(f"file:{gi}?mode=ro", uri=True)
    lens = dict(con.execute("select node_id, length(seq) from nodes").fetchall())
    shards = {}
    n_t = n_cells = n_bad = n_out = n_span_other = 0
    bad_examples = []
    for line in open(d / "variant_summary.ndjson"):
        r = json.loads(line)
        f = d / f"shard_{r['shard_index']:05d}_data.npy"
        if f not in shards:
            shards[f] = np.load(f)
        x = shards[f][r["index_within_shard"]]
        lay = r["site_layout"]
        own = own_node_mask(x, r["start"], lens[r["node_id"]], lay["insertion_slots"])
        a = r["anchor_column"]
        site = range(a, a + lay["insertion_slots"] + lay["span"])
        n_t += 1
        for ri, row in enumerate(r["rows"]):
            for ci, c in enumerate(row["columns"]):
                if c is None:
                    assert x[4, ri, ci] == 0
                    continue
                truth = c["node_id"] == r["node_id"] and c["mapping_index"] == row["anchor_mapping_index"]
                if ci in site and not truth:
                    n_span_other += 1  # multi-node deletion span on a further node
                n_cells += 1
                n_out += not truth
                if truth != own[ri, ci]:
                    n_bad += 1
                    if len(bad_examples) < 5:
                        bad_examples.append((r["candidate_id"], ri, ci, c, bool(own[ri, ci])))
    print(f"{d}: tensors {n_t}, covered cells {n_cells}, cells off the own node visit {n_out} "
          f"({n_out / max(n_cells, 1):.3f}), rule disagreements {n_bad}, site-span cells on a further node {n_span_other}")
    for e in bad_examples:
        print("  disagreement", e)


if __name__ == "__main__":
    main()

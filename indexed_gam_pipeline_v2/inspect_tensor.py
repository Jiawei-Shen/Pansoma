#!/usr/bin/env python3
"""Dump a few rows of every tensor in an output folder as aligned read/graph text."""
import argparse
import json
from pathlib import Path

import numpy as np

BASES = {0: ".", 1: "A", 2: "C", 3: "G", 4: "T", 5: "N", 6: "-"}
STRANDS = {0: "?", 1: "+", 2: "-"}


def export(folder, output, per_class=3):
    folder, output = Path(folder), Path(output)
    lines = ["# '.' = padding; '-' = alignment gap; '^' = site columns; allele = channel 2, the site allele the row carries"]
    for line in (folder / "variant_summary.ndjson").read_text().splitlines():
        meta = json.loads(line)
        tensor = np.load(folder / f"shard_{meta['shard_index']:05d}_data.npy", mmap_mode="r")[meta["index_within_shard"]]
        alleles = ", ".join(f"{label}={cid}" for label, cid in meta.get("allele_labels", {}).items())
        lines.append(f"\n{meta['candidate_id']} site_coverage={meta.get('site_coverage')} "
                     f"blocks={meta.get('site_counts')} alleles: {alleles}")
        lo, hi = meta["candidate_columns"]
        lines.append("anchor " + " " * lo + "^" * (hi - lo))
        counts = {}
        labels = [g["allele"] for g in meta.get("row_groups", []) for _ in range(g["start_row"], g["end_row"])]
        for ri in range(meta["selected_alignments"]):
            block = labels[ri] if ri < len(labels) else "row"
            counts[block] = counts.get(block, 0) + 1
            if counts[block] > per_class:
                continue
            covered = tensor[7, ri] != 0
            strand = STRANDS[int(tensor[7, ri][covered][0])] if covered.any() else "?"
            lines.append(f"row={ri} block={block} strand={strand}")
            lines.append("read   " + "".join(BASES[int(b)] for b in tensor[0, ri]))
            lines.append("graph  " + "".join(BASES[int(b)] for b in tensor[5, ri]))
            lines.append("allele " + "".join(BASES[int(b)] for b in tensor[2, ri]))
            lines.append("paths  " + " ".join(str(int(c)) for c in tensor[6, ri]))
    output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder")
    parser.add_argument("output")
    parser.add_argument("--per-class", type=int, default=3, help="rows shown per row block (A1.., REF, OTHER)")
    args = parser.parse_args()
    export(args.folder, args.output, args.per_class)

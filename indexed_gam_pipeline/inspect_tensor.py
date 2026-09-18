#!/usr/bin/env python3
"""Export inspectable paired read/reference rows from a candidate-v2/v3 output folder."""
import argparse
import json
from pathlib import Path

import numpy as np

BASES = {0: ".", 1: "A", 2: "C", 3: "G", 4: "T", 5: "N", 6: "-"}


def export(folder, output, per_class=3):
    folder, output = Path(folder), Path(output)
    manifest = json.loads((folder / "manifest.json").read_text())
    if manifest["tensor_format_version"] not in ("indexed-gam-candidate-v2", "indexed-gam-candidate-v3"):
        raise ValueError("Only candidate-v2/v3 tensors are supported")
    lines = ["# '.' = padding; '-' = alignment gap; '^' = candidate region"]
    for line in (folder / "variant_summary.ndjson").read_text().splitlines():
        meta = json.loads(line)
        tensor = np.load(folder / f"shard_{meta['shard_index']:05d}_data.npy", mmap_mode="r")[meta["index_within_shard"]]
        lines.append(f"\n{meta['candidate_id']} coverage={meta['coverage']} "
                     f"ALT/REF/other={meta['alt_count']}/{meta['ref_count']}/{meta['other_count']}")
        lo, hi = meta["candidate_columns"]
        lines.append("anchor " + " "*lo + "^"*(hi-lo))
        counts = {}
        for ri in range(meta["selected_alignments"]):
            support = meta.get("rows", [{}]*meta["selected_alignments"])[ri].get("support", "row")
            counts[support] = counts.get(support, 0)+1
            if counts[support] > per_class:
                continue
            lines.append(f"row={ri} support={support}")
            lines.append("read   " + "".join(BASES[int(b)] for b in tensor[0,ri]))
            lines.append("graph  " + "".join(BASES[int(b)] for b in tensor[5,ri]))
            if tensor.shape[0] == 7:
                lines.append("walks  " + " ".join(str(int(c)) for c in tensor[6,ri]))
    output.write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder")
    parser.add_argument("output")
    parser.add_argument("--per-class", type=int, default=3)
    args = parser.parse_args()
    export(args.folder, args.output, args.per_class)

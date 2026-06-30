#!/usr/bin/env python3
"""
Build a new JSON from:
  - merged node JSON (list of node dicts OR dict with "nodes"),
  - GFA file,
  - node IDs loaded from a binary .idx file (provided loader).

Rules:
  * If node_id exists in merged JSON, copy its FULL record as-is.
  * Otherwise, read its sequence from GFA S-lines and output ONLY:
        {"node_id": "<id>", "sequence": "<seq>"}
  * Output is a JSON LIST of node records, in the order read from .idx.

Usage:
  python build_nodes_from_idx.py \
      --merged merged_nodes.json[.gz] \
      --gfa graph.gfa[.gz] \
      --idx nodes.idx \
      --out new_nodes.json
"""

import argparse
import gzip
import json
import struct
import sys
from typing import Any, Dict, List, Set

# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers

def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def load_merged_map(path: str, id_key: str = "node_id") -> Dict[str, Dict[str, Any]]:
    """
    Accepts:
      - a JSON list of node dicts, or
      - a dict with key "nodes" -> list of node dicts
    Returns: { node_id(str) : node_dict }
    """
    with open_maybe_gzip(path, "rt") as f:
        data = json.load(f)

    nodes = data.get("nodes") if isinstance(data, dict) else data
    if not isinstance(nodes, list):
        raise ValueError("Merged JSON must be a list of nodes or a dict with 'nodes'.")

    out: Dict[str, Dict[str, Any]] = {}
    for n in nodes:
        if isinstance(n, dict) and id_key in n:
            out[str(n[id_key])] = n
    return out

def scan_gfa_for_sequences(gfa_path: str, wanted_ids: Set[str]) -> Dict[str, str]:
    """
    Stream the GFA once and collect sequences for wanted node IDs.
    GFA S-line format: S <node_id> <sequence> ...
    """
    seqs: Dict[str, str] = {}
    if not wanted_ids:
        return seqs

    with open_maybe_gzip(gfa_path, "rt") as f:
        for line in f:
            if not line or line[0] != "S":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            nid, seq = parts[1], parts[2]
            if nid in wanted_ids:
                seqs[nid] = "" if seq == "*" else seq
                if len(seqs) == len(wanted_ids):  # all found
                    break
    return seqs

# ─────────────────────────────────────────────────────────────────────────────
# Provided .idx loader (unchanged)

def load_index(idx_path):
    """Loads node IDs from a binary .idx file. (Unchanged)"""
    node_index = {}
    try:
        with open(idx_path, "rb") as f:
            blocks_num_bytes = f.read(4)
            if not blocks_num_bytes or len(blocks_num_bytes) < 4:
                print(f"Error: Could not read blocks_num from {idx_path}.", file=sys.stderr)
                return {}
            blocks_num, = struct.unpack("<I", blocks_num_bytes)
            for i in range(blocks_num):
                # Read only the relevant parts of the header to get the block_id
                header_data = f.read(22)  # Read the full header
                if len(header_data) < 22:
                    break
                block_id, _, _, _, metadata_len = struct.unpack("<I Q I I H", header_data)
                if metadata_len > 0:
                    f.read(metadata_len)  # Skip metadata
                node_index[block_id] = True  # We only care about the existence of the node ID
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}
    return node_index

# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    ap = argparse.ArgumentParser(description="Copy full records from merged JSON if present; otherwise add node_id+sequence from GFA; node IDs come from .idx.")
    ap.add_argument("--merged", required=True, help="Merged node JSON (list or dict with 'nodes'). Supports .gz")
    ap.add_argument("--gfa", required=True, help="GFA file (.gfa or .gfa.gz)")
    ap.add_argument("--idx", required=True, help="Binary .idx file containing node IDs")
    ap.add_argument("--out", required=True, help="Output JSON (list of nodes)")
    ap.add_argument("--id-key", default="node_id", help="Key for node id in JSON (default: node_id)")
    ap.add_argument("--seq-key", default="sequence", help="Key for sequence in JSON (default: sequence)")
    args = ap.parse_args()

    # Load IDs from .idx (preserve iteration order from the dict insertion)
    idx_map = load_index(args.idx)
    node_ids: List[str] = [str(k) for k in idx_map.keys()]

    # Load merged nodes
    merged_map = load_merged_map(args.merged, id_key=args.id_key)

    # Determine which IDs are missing and need GFA lookup
    missing: Set[str] = {nid for nid in node_ids if nid not in merged_map}

    # Scan GFA for missing sequences
    gfa_seqs = scan_gfa_for_sequences(args.gfa, missing)

    # Build output in .idx order
    out_nodes: List[Dict[str, Any]] = []
    for nid in node_ids:
        if nid in merged_map:
            out_nodes.append(merged_map[nid])  # full record as-is
        else:
            out_nodes.append({
                args.id_key: nid,
                args.seq_key: gfa_seqs.get(nid, "")
            })

    # Write JSON list
    with open(args.out, "w") as fo:
        json.dump(out_nodes, fo, ensure_ascii=False)

    # Summary
    copied_full = len(node_ids) - len(missing)
    added_from_gfa = sum(1 for nid in node_ids if nid in gfa_seqs)
    missing_seq = sum(1 for nid in node_ids if (nid not in merged_map and nid not in gfa_seqs))
    print(f"[summary] total:{len(node_ids)} copied_full:{copied_full} added_from_gfa:{added_from_gfa} missing_seq:{missing_seq}", file=sys.stderr)

if __name__ == "__main__":
    main()

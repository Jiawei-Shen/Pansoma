#!/usr/bin/env python3
"""
Build a new JSON from:
  - merged node JSON (list of node dicts OR dict with "nodes"),
  - GFA file,
  - node IDs loaded from one or more binary .idx files (provided loader).

Rules:
  * If node_id exists in merged JSON, copy its FULL record as-is (optionally add chrom if missing).
  * Otherwise, read its sequence from GFA S-lines and output ONLY:
        {"node_id": "<id>", "sequence": "<seq>"}  (optionally add chrom if requested)
  * Output is a JSON LIST of node records, in the order of node IDs read from idx files.
    If multiple --idx are provided, we use the UNION of node IDs, preserving first-seen order
    across the idx files (idx1 order first, then new IDs from idx2, ...).

Usage:
  python build_nodes_from_idx.py \
      --merged merged_nodes.json[.gz] \
      --gfa graph.gfa[.gz] \
      --idx nodes1.idx nodes2.idx \
      --out new_nodes.json \
      --add-chrom chr1
"""

import argparse
import gzip
import json
import struct
import sys
import os
from typing import Any, Dict, List, Set, Optional


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
# Provided .idx loader (latest 30B entries)
#
# Latest .idx entry (30B): <I Q I I H I I>
#   nid, offset, block_size, n_records, flags, R, C
#
# We only need nid, but we parse strictly and preserve file order.
# ─────────────────────────────────────────────────────────────────────────────

IDX_ENTRY_PACK_LATEST = struct.Struct("<I Q I I H I I")
IDX_ENTRY_SIZE_LATEST = IDX_ENTRY_PACK_LATEST.size  # 30


def load_index_ids(idx_path: str) -> List[str]:
    """Loads node IDs from a binary .idx file (latest format: 30B entries), preserving file order."""
    ids: List[str] = []
    try:
        with open(idx_path, "rb") as f:
            hdr = f.read(4)
            if not hdr or len(hdr) < 4:
                print(f"Error: Could not read blocks_num from {idx_path}.", file=sys.stderr)
                return []
            blocks_num, = struct.unpack("<I", hdr)

            # Strict size check for latest format
            try:
                f.seek(0, os.SEEK_END)
                size = f.tell()
                expected = 4 + blocks_num * IDX_ENTRY_SIZE_LATEST
                if size != expected:
                    print(
                        f"Error: IDX size mismatch for latest format: file={size}, expected={expected} "
                        f"(count={blocks_num}, entry={IDX_ENTRY_SIZE_LATEST}).",
                        file=sys.stderr
                    )
                    return []
            except Exception:
                # If size check fails for some reason, proceed without it (still parse carefully)
                pass

            f.seek(4, os.SEEK_SET)
            for i in range(blocks_num):
                rec = f.read(IDX_ENTRY_SIZE_LATEST)
                if len(rec) != IDX_ENTRY_SIZE_LATEST:
                    print(f"Error: Truncated IDX at entry {i+1} in {idx_path}.", file=sys.stderr)
                    return []
                block_id, _off, _blk_sz, _nrec, _flags, _R, _C = IDX_ENTRY_PACK_LATEST.unpack(rec)
                ids.append(str(block_id))
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return []
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return []
    return ids


def union_node_ids_from_multiple_idx(idx_paths: List[str]) -> List[str]:
    """
    Union node IDs across multiple idx files, preserving first-seen order:
      - iterate idx_paths in the given order
      - within each idx, preserve file order
      - keep only first occurrence of a node id across all idx files
    """
    seen: Set[str] = set()
    out: List[str] = []
    for p in idx_paths:
        ids = load_index_ids(p)
        for nid in ids:
            if nid not in seen:
                seen.add(nid)
                out.append(nid)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    ap = argparse.ArgumentParser(
        description="Copy full records from merged JSON if present; otherwise add node_id+sequence from GFA; node IDs come from one or more .idx files."
    )
    ap.add_argument("--merged", required=True, help="Merged node JSON (list or dict with 'nodes'). Supports .gz")
    ap.add_argument("--gfa", required=True, help="GFA file (.gfa or .gfa.gz)")

    # UPDATED: accept multiple idx files
    ap.add_argument("--idx", required=True, nargs="+", help="One or more binary .idx files (union of node IDs is used)")

    ap.add_argument("--out", required=True, help="Output JSON (list of nodes)")
    ap.add_argument("--id-key", default="node_id", help="Key for node id in JSON (default: node_id)")
    ap.add_argument("--seq-key", default="sequence", help="Key for sequence in JSON (default: sequence)")

    # NEW: manually add chrom if missing
    ap.add_argument(
        "--add-chrom",
        default=None,
        help='If set (e.g. "chr1"), add this chrom value to each output node record ONLY when chrom key is missing.'
    )
    ap.add_argument(
        "--chrom-key",
        default="chrom",
        help='Chromosome key name to use when --add-chrom is provided (default: "chrom").'
    )

    args = ap.parse_args()

    # Load IDs from idx files (union, preserve order)
    node_ids: List[str] = union_node_ids_from_multiple_idx(args.idx)
    if not node_ids:
        print("ERROR: No node IDs loaded from --idx inputs.", file=sys.stderr)
        sys.exit(2)

    # Load merged nodes
    merged_map = load_merged_map(args.merged, id_key=args.id_key)

    # Determine which IDs are missing and need GFA lookup
    missing: Set[str] = {nid for nid in node_ids if nid not in merged_map}

    # Scan GFA for missing sequences
    gfa_seqs = scan_gfa_for_sequences(args.gfa, missing)

    # Helper: add chrom only if requested and missing
    def maybe_add_chrom(rec: Dict[str, Any]) -> Dict[str, Any]:
        if args.add_chrom is None:
            return rec
        if args.chrom_key not in rec:
            rec[args.chrom_key] = args.add_chrom
        return rec

    # Build output in union order
    out_nodes: List[Dict[str, Any]] = []
    for nid in node_ids:
        if nid in merged_map:
            # full record as-is (but optionally add chrom if missing)
            rec = merged_map[nid]
            maybe_add_chrom(rec)
            out_nodes.append(rec)
        else:
            rec = {
                args.id_key: nid,
                args.seq_key: gfa_seqs.get(nid, "")
            }
            maybe_add_chrom(rec)
            out_nodes.append(rec)

    # Write JSON list
    with open(args.out, "w", encoding="utf-8") as fo:
        json.dump(out_nodes, fo, ensure_ascii=False)

    # Summary
    copied_full = len(node_ids) - len(missing)
    added_from_gfa = sum(1 for nid in node_ids if nid in gfa_seqs)
    missing_seq = sum(1 for nid in node_ids if (nid not in merged_map and nid not in gfa_seqs))
    print(
        f"[summary] total:{len(node_ids)} copied_full:{copied_full} added_from_gfa:{added_from_gfa} missing_seq:{missing_seq}",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Build a JSON list of nodes from:
  - GFA file
  - one or more binary .idx files
  - optional input JSON containing richer node records

Behavior:
  * Read node IDs from all --idx files
  * Take the UNION of node IDs, preserving first-seen order
  * Scan the GFA and collect sequence for every wanted node ID
  * If --input_json is provided and a node exists there:
      - copy its full record
      - remove excluded fields (default: genomead_af)
      - if sequence is missing, fill from GFA
      - if --add-chrom is set and chrom is missing, add it
  * Otherwise output:
        {"node_id": "<id>", "sequence": "<seq>"}
    and, if requested:
        {"node_id": "<id>", "sequence": "<seq>", "chrom": "<chr>"}

Notes:
  * If a node ID from idx is not found in the GFA, its sequence is set to "".
  * input_json may be:
      - a JSON list of node dicts
      - a dict with key "nodes" -> list of node dicts
  * genomead_af is excluded by default.

Usage:
  python build_nodes_from_idx_gfa_plus_json.py \
      --gfa graph.gfa.gz \
      --idx a.idx b.idx \
      --input_json merged_nodes.json.gz \
      --out out_nodes.json \
      --add-chrom chr1
"""

import argparse
import gzip
import json
import os
import struct
import sys
from typing import Any, Dict, List, Set


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers

def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


# ─────────────────────────────────────────────────────────────────────────────
# Input JSON loader

def load_input_json_map(path: str, id_key: str = "node_id") -> Dict[str, Dict[str, Any]]:
    """
    Accept:
      - a JSON list of node dicts
      - a dict with key "nodes" -> list of node dicts

    Return:
      { str(node_id): node_record_dict }
    """
    log(f"[json] loading input_json: {path}")
    with open_maybe_gzip(path, "rt") as f:
        data = json.load(f)

    if isinstance(data, dict):
        nodes = data.get("nodes")
    else:
        nodes = data

    if not isinstance(nodes, list):
        raise ValueError("input_json must be a list of node dicts or a dict with key 'nodes'.")

    out: Dict[str, Dict[str, Any]] = {}
    skipped = 0
    for rec in nodes:
        if not isinstance(rec, dict):
            skipped += 1
            continue
        if id_key not in rec:
            skipped += 1
            continue
        out[str(rec[id_key])] = rec

    log(f"[json] loaded records={len(out)} skipped={skipped}")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# IDX loader
#
# Latest .idx entry (30B): <I Q I I H I I>
#   nid, offset, block_size, n_records, flags, R, C
# We only need nid, but parse strictly and preserve file order.
# ─────────────────────────────────────────────────────────────────────────────

IDX_ENTRY_PACK_LATEST = struct.Struct("<I Q I I H I I")
IDX_ENTRY_SIZE_LATEST = IDX_ENTRY_PACK_LATEST.size  # 30


def load_index_ids(idx_path: str) -> List[str]:
    """Load node IDs from a binary .idx file, preserving file order."""
    ids: List[str] = []
    try:
        log(f"[idx] reading: {idx_path}")
        with open(idx_path, "rb") as f:
            hdr = f.read(4)
            if len(hdr) != 4:
                log(f"Error: Could not read blocks_num from {idx_path}.")
                return []

            (blocks_num,) = struct.unpack("<I", hdr)

            try:
                f.seek(0, os.SEEK_END)
                size = f.tell()
                expected = 4 + blocks_num * IDX_ENTRY_SIZE_LATEST
                if size != expected:
                    log(
                        f"Error: IDX size mismatch for latest format: file={size}, expected={expected} "
                        f"(count={blocks_num}, entry={IDX_ENTRY_SIZE_LATEST}) in {idx_path}."
                    )
                    return []
            except Exception:
                pass

            f.seek(4, os.SEEK_SET)
            for i in range(blocks_num):
                rec = f.read(IDX_ENTRY_SIZE_LATEST)
                if len(rec) != IDX_ENTRY_SIZE_LATEST:
                    log(f"Error: Truncated IDX at entry {i+1} in {idx_path}.")
                    return []
                block_id, _off, _blk_sz, _nrec, _flags, _R, _C = IDX_ENTRY_PACK_LATEST.unpack(rec)
                ids.append(str(block_id))

        log(f"[idx] done: {idx_path}  blocks={blocks_num} ids={len(ids)}")

    except FileNotFoundError:
        log(f"Error: IDX file not found: {idx_path}")
        return []
    except Exception as e:
        log(f"Unexpected error while reading IDX file {idx_path}: {e}")
        return []

    return ids


def union_node_ids_from_multiple_idx(idx_paths: List[str]) -> List[str]:
    """
    Union node IDs across multiple idx files, preserving first-seen order.
    """
    seen: Set[str] = set()
    out: List[str] = []

    log(f"[idx] start loading {len(idx_paths)} idx file(s)")
    for idx_path in idx_paths:
        ids = load_index_ids(idx_path)
        before = len(out)
        for nid in ids:
            if nid not in seen:
                seen.add(nid)
                out.append(nid)
        added = len(out) - before
        log(f"[idx] union updated: +{added} new ids, total={len(out)}")

    return out


# ─────────────────────────────────────────────────────────────────────────────
# GFA scan

def scan_gfa_for_sequences(gfa_path: str, wanted_ids: Set[str]) -> Dict[str, str]:
    """
    Stream the GFA once and collect sequences for wanted node IDs.
    GFA S-line format: S <node_id> <sequence> ...
    """
    seqs: Dict[str, str] = {}
    if not wanted_ids:
        return seqs

    log(f"[gfa] scanning: {gfa_path}  wanted_ids={len(wanted_ids)}")
    with open_maybe_gzip(gfa_path, "rt") as f:
        for line in f:
            if not line or line[0] != "S":
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            nid, seq = parts[1], parts[2]
            if nid in wanted_ids and nid not in seqs:
                seqs[nid] = "" if seq == "*" else seq
                if len(seqs) == len(wanted_ids):
                    break

    log(f"[gfa] done: found_seq={len(seqs)} / {len(wanted_ids)}")
    return seqs


# ─────────────────────────────────────────────────────────────────────────────
# Record builder helpers

def normalize_exclude_fields(raw_fields: List[str]) -> Set[str]:
    return {x.strip() for x in raw_fields if x and x.strip()}


def build_record(
    nid: str,
    seqs: Dict[str, str],
    input_map: Dict[str, Dict[str, Any]],
    id_key: str,
    seq_key: str,
    chrom_key: str,
    add_chrom: str,
    exclude_fields: Set[str],
) -> Dict[str, Any]:
    """
    Build one output record.
    """
    gfa_seq = seqs.get(nid, "")

    if nid in input_map:
        # copy full record
        src = input_map[nid]
        rec = dict(src)

        # remove excluded fields exactly by key name
        for k in list(rec.keys()):
            if k in exclude_fields:
                rec.pop(k, None)

        # ensure node_id key exists and matches requested output key
        rec[id_key] = nid

        # fill sequence if missing / null / empty
        if seq_key not in rec or rec[seq_key] is None or rec[seq_key] == "":
            rec[seq_key] = gfa_seq

        # optionally add chrom if missing
        if add_chrom is not None and chrom_key not in rec:
            rec[chrom_key] = add_chrom

        return rec

    # fallback: GFA-only record
    rec = {
        id_key: nid,
        seq_key: gfa_seq,
    }
    if add_chrom is not None:
        rec[chrom_key] = add_chrom
    return rec


# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    ap = argparse.ArgumentParser(
        description="Build JSON node records from GFA using node IDs loaded from one or more .idx files, optionally enriching from input_json."
    )
    ap.add_argument("--gfa", required=True, help="Input GFA file (.gfa or .gfa.gz)")
    ap.add_argument(
        "--idx",
        required=True,
        nargs="+",
        help="One or more binary .idx files; node IDs are unioned in first-seen order",
    )
    ap.add_argument(
        "--input_json",
        default=None,
        help="Optional JSON/JSON.gz with richer node records (list or dict with 'nodes')",
    )
    ap.add_argument("--out", required=True, help="Output JSON file")
    ap.add_argument(
        "--add-chrom",
        default=None,
        help='If set (e.g. "chr1"), add {"chrom": "..."} to output records only when missing',
    )
    ap.add_argument(
        "--id-key",
        default="node_id",
        help='Output node ID key name (default: "node_id")',
    )
    ap.add_argument(
        "--seq-key",
        default="sequence",
        help='Output sequence key name (default: "sequence")',
    )
    ap.add_argument(
        "--chrom-key",
        default="chrom",
        help='Output chrom key name when --add-chrom is used (default: "chrom")',
    )
    ap.add_argument(
        "--exclude-fields",
        nargs="*",
        default=["genomead_af"],
        help='Fields to remove from copied input_json records (default: genomead_af)',
    )
    ap.add_argument(
        "--indent",
        type=int,
        default=None,
        help="Optional JSON indent for pretty output",
    )

    args = ap.parse_args()

    exclude_fields = normalize_exclude_fields(args.exclude_fields)

    log("[main] step 1/4: loading node IDs from idx")
    node_ids = union_node_ids_from_multiple_idx(args.idx)
    if not node_ids:
        log("ERROR: No node IDs loaded from --idx inputs.")
        sys.exit(2)

    input_map: Dict[str, Dict[str, Any]] = {}
    if args.input_json:
        log("[main] step 2/4: loading input_json")
        input_map = load_input_json_map(args.input_json, id_key=args.id_key)
    else:
        log("[main] step 2/4: no input_json provided; will build records from GFA only")

    log(f"[main] step 3/4: scanning GFA for {len(node_ids)} unique node IDs")
    wanted_ids = set(node_ids)
    seqs = scan_gfa_for_sequences(args.gfa, wanted_ids)

    log(f"[main] step 4/4: writing output JSON -> {args.out}")
    out_nodes: List[Dict[str, Any]] = []
    copied_from_input = 0
    built_from_gfa_only = 0

    for nid in node_ids:
        if nid in input_map:
            copied_from_input += 1
        else:
            built_from_gfa_only += 1

        rec = build_record(
            nid=nid,
            seqs=seqs,
            input_map=input_map,
            id_key=args.id_key,
            seq_key=args.seq_key,
            chrom_key=args.chrom_key,
            add_chrom=args.add_chrom,
            exclude_fields=exclude_fields,
        )
        out_nodes.append(rec)

    with open(args.out, "w", encoding="utf-8") as fo:
        json.dump(out_nodes, fo, ensure_ascii=False, indent=args.indent)

    found_seq = sum(1 for nid in node_ids if nid in seqs)
    missing_seq = len(node_ids) - found_seq
    log(
        f"[summary] total:{len(node_ids)} copied_from_input:{copied_from_input} "
        f"built_from_gfa_only:{built_from_gfa_only} found_seq:{found_seq} missing_seq:{missing_seq}"
    )


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Compare node IDs from strict-format .idx file(s)
against node IDs extracted from the graph_nodes_json column in a TSV.

IDX format:
  header:
    <I   blocks_num
  each entry (30 bytes):
    <I Q I I H I I
    nid, offset, block_size, n_records, flags, R, C

TSV:
  must contain a column named graph_nodes_json
  each cell is a JSON list like:
    [{"id":"90223","sequence":"AAAA"}, {"id":"90224","sequence":"GA"}]

Outputs:
  *.idx_ids.txt
  *.tsv_graph_node_ids.txt
  *.overlap.txt
  *.only_in_idx.txt
  *.only_in_tsv.txt
  *.summary.txt
"""

import argparse
import csv
import gzip
import json
import os
import struct
import sys
from typing import List, Set, Tuple


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


# ─────────────────────────────────────────────────────────────────────────────
# IDX loader (strict)

IDX_ENTRY_PACK_LATEST = struct.Struct("<I Q I I H I I")
IDX_ENTRY_SIZE_LATEST = IDX_ENTRY_PACK_LATEST.size  # 30


def load_index_ids(idx_path: str) -> List[str]:
    ids: List[str] = []
    try:
        log(f"[idx] reading: {idx_path}")
        with open(idx_path, "rb") as f:
            hdr = f.read(4)
            if len(hdr) != 4:
                log(f"ERROR: Could not read blocks_num from {idx_path}")
                return []

            (blocks_num,) = struct.unpack("<I", hdr)

            f.seek(0, os.SEEK_END)
            size = f.tell()
            expected = 4 + blocks_num * IDX_ENTRY_SIZE_LATEST
            if size != expected:
                log(
                    f"ERROR: IDX size mismatch: file={size}, expected={expected} "
                    f"(count={blocks_num}, entry={IDX_ENTRY_SIZE_LATEST}) in {idx_path}"
                )
                return []

            f.seek(4, os.SEEK_SET)
            for i in range(blocks_num):
                rec = f.read(IDX_ENTRY_SIZE_LATEST)
                if len(rec) != IDX_ENTRY_SIZE_LATEST:
                    log(f"ERROR: Truncated IDX at entry {i+1} in {idx_path}")
                    return []
                nid, _off, _blk_sz, _nrec, _flags, _R, _C = IDX_ENTRY_PACK_LATEST.unpack(rec)
                ids.append(str(nid))

        log(f"[idx] done: {idx_path}  blocks={blocks_num} ids={len(ids)}")

    except FileNotFoundError:
        log(f"ERROR: IDX file not found: {idx_path}")
        return []
    except Exception as e:
        log(f"ERROR: Unexpected error while reading IDX file {idx_path}: {e}")
        return []

    return ids


def union_node_ids_from_multiple_idx(idx_paths: List[str]) -> List[str]:
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
# TSV graph_nodes_json loader

def load_graph_node_ids_from_tsv(tsv_path: str) -> Tuple[Set[str], List[str], int, int]:
    """
    Returns:
      - set of unique node ids
      - ordered unique node ids (first-seen)
      - number of rows processed
      - number of rows with parse problems
    """
    seen: Set[str] = set()
    ordered: List[str] = []
    rows_processed = 0
    rows_bad = 0

    with open_maybe_gzip(tsv_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("TSV appears empty or has no header.")

        if "graph_nodes_json" not in reader.fieldnames:
            raise ValueError("TSV does not contain a 'graph_nodes_json' column.")

        for row in reader:
            rows_processed += 1
            raw = row.get("graph_nodes_json", "")

            if raw is None:
                continue
            raw = raw.strip()
            if raw == "":
                continue

            try:
                arr = json.loads(raw)
                if not isinstance(arr, list):
                    rows_bad += 1
                    continue

                for item in arr:
                    if not isinstance(item, dict):
                        continue
                    nid = item.get("id")
                    if nid is None:
                        continue
                    nid = str(nid)
                    if nid not in seen:
                        seen.add(nid)
                        ordered.append(nid)

            except Exception:
                rows_bad += 1
                continue

    return seen, ordered, rows_processed, rows_bad


# ─────────────────────────────────────────────────────────────────────────────
# Output helpers

def write_id_list(path: str, ids: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for nid in ids:
            f.write(f"{nid}\n")


def write_summary(path: str, summary: List[Tuple[str, int]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for k, v in summary:
            f.write(f"{k}\t{v}\n")


# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    ap = argparse.ArgumentParser(
        description="Compare strict-format .idx node IDs against node IDs extracted from graph_nodes_json in a TSV."
    )
    ap.add_argument(
        "--idx",
        required=True,
        nargs="+",
        help="One or more .idx files"
    )
    ap.add_argument(
        "--tsv",
        required=True,
        help="Input TSV containing graph_nodes_json column"
    )
    ap.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix"
    )
    args = ap.parse_args()

    log("[main] step 1/3: loading node IDs from idx")
    idx_ids_ordered = union_node_ids_from_multiple_idx(args.idx)
    if not idx_ids_ordered:
        log("ERROR: No node IDs loaded from --idx inputs")
        sys.exit(2)
    idx_ids_set = set(idx_ids_ordered)

    log("[main] step 2/3: loading node IDs from TSV graph_nodes_json")
    tsv_ids_set, tsv_ids_ordered, rows_processed, rows_bad = load_graph_node_ids_from_tsv(args.tsv)
    if not tsv_ids_ordered:
        log("ERROR: No node IDs extracted from TSV graph_nodes_json")
        sys.exit(3)

    log("[main] step 3/3: comparing")
    overlap = [nid for nid in idx_ids_ordered if nid in tsv_ids_set]
    only_in_idx = [nid for nid in idx_ids_ordered if nid not in tsv_ids_set]
    only_in_tsv = [nid for nid in tsv_ids_ordered if nid not in idx_ids_set]

    out_prefix = args.out_prefix
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    write_id_list(out_prefix + ".idx_ids.txt", idx_ids_ordered)
    write_id_list(out_prefix + ".tsv_graph_node_ids.txt", tsv_ids_ordered)
    write_id_list(out_prefix + ".overlap.txt", overlap)
    write_id_list(out_prefix + ".only_in_idx.txt", only_in_idx)
    write_id_list(out_prefix + ".only_in_tsv.txt", only_in_tsv)

    summary = [
        ("idx_unique_ids", len(idx_ids_ordered)),
        ("tsv_graph_node_unique_ids", len(tsv_ids_ordered)),
        ("overlap", len(overlap)),
        ("only_in_idx", len(only_in_idx)),
        ("only_in_tsv", len(only_in_tsv)),
        ("tsv_rows_processed", rows_processed),
        ("tsv_rows_bad_json", rows_bad),
    ]
    write_summary(out_prefix + ".summary.txt", summary)

    log(
        f"[summary] idx_unique_ids={len(idx_ids_ordered)} "
        f"tsv_graph_node_unique_ids={len(tsv_ids_ordered)} "
        f"overlap={len(overlap)} "
        f"only_in_idx={len(only_in_idx)} "
        f"only_in_tsv={len(only_in_tsv)} "
        f"tsv_rows_processed={rows_processed} "
        f"tsv_rows_bad_json={rows_bad}"
    )


if __name__ == "__main__":
    main()
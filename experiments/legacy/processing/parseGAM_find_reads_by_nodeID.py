#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
import time
import json
import csv
import concurrent.futures
import vg_pb2
from typing import Dict, Set, Iterable, Tuple
from google.protobuf.json_format import MessageToDict  # optional import; not required below

# ─────────────────────────────────────────────────────────────────────────────
# Basic varint + GAM group parsing (works for gz or plain .gam)

def read_varint(stream) -> int:
    result = 0
    shift = 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError("Unexpected EOF while reading varint")
        byte = b[0]
        result |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            break
        shift += 7
    return result

def is_gzipped(filename: str) -> bool:
    with open(filename, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def parse_gam_file_groups(filename: str, expected_tag: str = "GAM") -> Iterable[Tuple[int, list]]:
    open_func = gzip.open if is_gzipped(filename) else open
    with open_func(filename, 'rb') as f:
        group_number = 0
        while True:
            try:
                group_count = read_varint(f)
            except EOFError:
                break
            group_number += 1
            if group_count == 0:
                continue
            try:
                tag_size = read_varint(f)
                tag = f.read(tag_size).decode("utf-8")
            except (EOFError, UnicodeDecodeError):
                break
            if tag != expected_tag:
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                        f.seek(msg_size, os.SEEK_CUR)
                    except EOFError:
                        break
                continue

            messages = []
            for _ in range(group_count - 1):
                try:
                    msg_size = read_varint(f)
                    msg_bytes = f.read(msg_size)
                    if len(msg_bytes) == msg_size:
                        messages.append(msg_bytes)
                except EOFError:
                    break
            yield (group_number, messages)

# ─────────────────────────────────────────────────────────────────────────────
# Per-group worker

def _process_group_for_nodes(args):
    group_number, messages, target_nodes, min_mapq = args
    hits: Dict[int, Set[str]] = {nid: set() for nid in target_nodes}  # restrict keys to requested nodes
    total_reads = 0

    aln = vg_pb2.Alignment()
    for msg_bytes in messages:
        aln.Clear()
        aln.ParseFromString(msg_bytes)

        # filter by mapping quality if requested
        if aln.mapping_quality < min_mapq:
            continue

        # collect which requested nodes this read touches (avoid duplicates within the read)
        touched = set()
        if aln.path.mapping:
            for m in aln.path.mapping:
                nid = m.position.node_id
                if nid and nid in hits:
                    touched.add(nid)

        if not touched:
            continue

        total_reads += 1
        read_name = aln.name if aln.name else f"READ_{total_reads}"
        for nid in touched:
            hits[nid].add(read_name)

    # Trim empty sets to save bandwidth
    hits = {nid: names for nid, names in hits.items() if names}
    return group_number, hits, total_reads

# ─────────────────────────────────────────────────────────────────────────────
# Utilities

def parse_nodes_list(nodes_csv: str) -> Set[int]:
    out = set()
    if not nodes_csv:
        return out
    for tok in nodes_csv.split(','):
        tok = tok.strip()
        if not tok:
            continue
        try:
            out.add(int(tok))
        except ValueError:
            pass
    return out

def parse_nodes_file(path: str) -> Set[int]:
    out = set()
    if not path:
        return out
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            # allow comments after a space or '#'
            s = s.split('#', 1)[0].strip()
            if not s:
                continue
            tok = s.split()[0]
            try:
                out.add(int(tok))
            except ValueError:
                continue
    return out

# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    ap = argparse.ArgumentParser(
        description="Extract read names aligned to specified node IDs from a GAM (gz or plain). "
                    "Outputs node_id, count, reads as JSON/TSV/CSV."
    )
    ap.add_argument("gam", help="Input .gam or .gam.gz")
    ap.add_argument("--nodes", help="Comma-separated node IDs", default="")
    ap.add_argument("--nodes_file", help="File with node IDs (one per line)", default="")
    ap.add_argument("--min_mapq", type=int, default=0, help="Minimum mapping quality to include (default: 0)")
    ap.add_argument("--threads", type=int, default=4, help="Worker processes (default: 4)")
    ap.add_argument("--max_pending", type=int, default=16, help="Max groups in flight (default: 16)")
    ap.add_argument("--format", choices=["json", "tsv", "csv"], default="json", help="Output format (default: json)")
    ap.add_argument("--output", help="Output file path; default picks extension from --format")

    args = ap.parse_args()

    if not os.path.isfile(args.gam):
        sys.exit(f"Input not found: {args.gam}")

    targets = parse_nodes_list(args.nodes) | parse_nodes_file(args.nodes_file)
    if not targets:
        sys.exit("No node IDs provided. Use --nodes and/or --nodes_file.")

    out_path = args.output
    if not out_path:
        base = os.path.splitext(os.path.basename(args.gam))[0]
        suf = {"json": ".reads_for_nodes.json", "tsv": ".reads_for_nodes.tsv", "csv": ".reads_for_nodes.csv"}[args.format]
        out_path = base + suf

    t0 = time.time()
    print(f"[INFO] Targets: {len(targets)} nodes")
    print(f"[INFO] Parsing: {args.gam}  (threads={args.threads}, min_mapq={args.min_mapq})")

    # master accumulator: keep only requested nodes to bound memory
    all_hits: Dict[int, Set[str]] = {nid: set() for nid in targets}
    total_reads_seen = 0
    groups_seen = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
        pending = []
        for group_number, messages in parse_gam_file_groups(args.gam):
            fut = ex.submit(_process_group_for_nodes, (group_number, messages, targets, args.min_mapq))
            pending.append(fut)

            if len(pending) >= args.max_pending:
                done, not_done = concurrent.futures.wait(pending, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    groups_seen += 1
                    _, partial_hits, group_reads = fut.result()
                    total_reads_seen += group_reads
                    for nid, names in partial_hits.items():
                        all_hits[nid].update(names)
                pending = list(not_done)

        for fut in concurrent.futures.as_completed(pending):
            groups_seen += 1
            _, partial_hits, group_reads = fut.result()
            total_reads_seen += group_reads
            for nid, names in partial_hits.items():
                all_hits[nid].update(names)

    elapsed = time.time() - t0
    print(f"[INFO] Done. Groups: {groups_seen}, reads scanned: {total_reads_seen:,} in {elapsed:.1f}s")
    print(f"[INFO] Writing: {out_path} ({args.format.upper()})")

    # Emit output
    if args.format == "json":
        # { "node_id": {"count": N, "reads": [..]}, ... }
        out_obj = {}
        for nid in sorted(all_hits.keys()):
            reads = sorted(all_hits[nid])
            out_obj[str(nid)] = {"count": len(reads), "reads": reads}
        with open(out_path, "w") as f:
            json.dump(out_obj, f, indent=2)
    else:
        # TSV/CSV with columns: node_id, count, reads
        # Put all reads for a node into one cell (joined by comma).
        # csv module will quote the field when necessary.
        dialect = "excel-tab" if args.format == "tsv" else "excel"
        with open(out_path, "w", newline="") as f:
            writer = csv.writer(f, dialect=dialect)
            writer.writerow(["node_id", "count", "reads"])
            for nid in sorted(all_hits.keys()):
                reads = sorted(all_hits[nid])
                writer.writerow([nid, len(reads), ",".join(reads)])

    # quick summary
    nonempty = sum(1 for s in all_hits.values() if s)
    total_found = sum(len(s) for s in all_hits.values())
    print(f"[INFO] Nodes requested: {len(all_hits)}; nodes with ≥1 read: {nonempty}; total unique read hits: {total_found:,}")
    print("[INFO] Done.")

if __name__ == "__main__":
    main()

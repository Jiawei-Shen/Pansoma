#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import sys
import concurrent.futures
import vg_pb2
import json
import pickle
from typing import Dict, Tuple, Iterable, Optional
from google.protobuf.json_format import MessageToDict


def read_varint(stream):
    result = 0
    shift = 0
    while True:
        b = stream.read(1)
        if len(b) == 0:
            raise EOFError("Unexpected EOF while reading varint")
        byte = b[0]
        result |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            break
        shift += 7
    return result


def is_gzipped(filename):
    with open(filename, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def parse_gam_file_groups(filename, expected_tag="GAM"):
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


def is_on_chromosome(alignment, chrom_name):
    return any(rp.name == chrom_name for rp in alignment.refpos)


def process_group_serialized(args):
    group_number, messages, chrom_name = args
    total = 0
    perfect = 0
    not_perfect = 0
    node_counts = {}  # node_id -> {"perfect": X, "not_perfect": Y}

    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)

        # Skip low mapping quality reads
        if alignment.mapping_quality <= 5:
            continue

        if chrom_name and not is_on_chromosome(alignment, chrom_name):
            continue

        total += 1

        for mapping in alignment.path.mapping:
            if not mapping.position.node_id:
                continue

            node_id = mapping.position.node_id
            is_perfect = True

            for edit in mapping.edit:
                # Perfect if no inserted sequence and from_length > 0 for all edits
                if edit.sequence or edit.from_length != edit.to_length:
                    is_perfect = False
                    break

            if node_id not in node_counts:
                node_counts[node_id] = {"perfect": 0, "not_perfect": 0}

            if is_perfect:
                node_counts[node_id]["perfect"] += 1
                perfect += 1
            else:
                node_counts[node_id]["not_perfect"] += 1
                not_perfect += 1

    return group_number, node_counts, total, perfect, not_perfect


def process_groups_pipeline(filename, threads, max_pending=10, chrom_name=""):
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        pending_futures = []
        for group in parse_gam_file_groups(filename):
            args = (group[0], group[1], chrom_name)
            future = executor.submit(process_group_serialized, args)
            pending_futures.append(future)

            if len(pending_futures) >= max_pending:
                done, not_done = concurrent.futures.wait(
                    pending_futures,
                    return_when=concurrent.futures.FIRST_COMPLETED
                )
                for completed in done:
                    yield completed.result()
                pending_futures = list(not_done)

        for future in concurrent.futures.as_completed(pending_futures):
            yield future.result()


# ─────────────────────────────────────────────────────────────────────────────
# NEW: GFA length parsing

def parse_gfa_lengths(gfa_path: str) -> Dict[int, int]:
    """
    Parse a GFA (v1-style) file and return a mapping: node_id(int) -> length(int).

    Strategy:
      - Only process 'S' segment lines.
      - Prefer LN:i:<len> tag if present.
      - Else, if sequence is not '*', use len(sequence).
      - If name is not an integer, the line is skipped (since this pipeline uses numeric node IDs).
    """
    if not gfa_path:
        return {}
    open_func = gzip.open if gfa_path.endswith(".gz") else open
    node_len: Dict[int, int] = {}
    total_s = 0
    with open_func(gfa_path, "rt") as fh:
        for line in fh:
            if not line or line[0] != 'S':
                continue
            # GFA1 S-line: S <name> <sequence> [opt fields...]
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            total_s += 1
            name = parts[1]
            seq = parts[2]
            # Optional tags start at parts[3:]
            length_val: Optional[int] = None
            # Prefer LN tag if available
            for field in parts[3:]:
                # LN:i:<len>
                if field.startswith("LN:i:"):
                    try:
                        length_val = int(field.split(":", 2)[2])
                    except ValueError:
                        pass
                    break
            if length_val is None and seq != "*":
                length_val = len(seq)
            # Convert name to int node_id if possible
            try:
                node_id = int(name)
            except ValueError:
                # Not a numeric node id; skip for this pipeline
                continue
            if length_val is None:
                # Unknown length, default to 0 (could change to None if preferred)
                length_val = 0
            node_len[node_id] = length_val
    if not node_len and total_s > 0:
        print("[warn] Parsed GFA but found no numeric node IDs in S-lines.")
    return node_len


def merge_lengths(counts: Dict[int, Dict[str, int]], lengths: Dict[int, int]) -> None:
    """
    In-place: add 'length' to each node entry in counts.
    Missing nodes get length=0.
    """
    missing = 0
    for nid, rec in counts.items():
        ln = lengths.get(nid, 0)
        rec["length"] = int(ln)
        if ln == 0 and nid not in lengths:
            missing += 1
    if lengths:
        print(f"Added lengths for {len(counts) - missing} nodes (missing={missing}).")


# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Parse a GAM file, count perfect/non-perfect mappings per node, optionally enrich with GFA node lengths.")
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=4, help="Number of worker processes (default: 4)")
    parser.add_argument("--max_pending", type=int, default=16, help="Max number of groups in flight (default: 16)")
    parser.add_argument("--chr", default="", help="Chromosome name to filter on (default: \"\")")
    parser.add_argument("--milestone", type=int, default=100_000_000, help="Number of reads between progress updates")
    parser.add_argument("--output_format", choices=["json", "pickle"], default="json", help="Output format (default: json)")
    parser.add_argument("--output", help="Output file name (overrides default)")
    parser.add_argument("--gfa", default="", help="Optional GFA file to add node lengths")
    args = parser.parse_args()
    start_time = time.perf_counter()

    total_count = 0
    perfect_count = 0
    not_perfect_count = 0
    milestone = args.milestone
    node_to_counts: Dict[int, Dict[str, int]] = {}

    for group_number, node_reads, t, p, np in process_groups_pipeline(
            args.filename, args.threads, args.max_pending, chrom_name=args.chr):

        total_count += t
        perfect_count += p
        not_perfect_count += np

        for node_id, counts in node_reads.items():
            if node_id not in node_to_counts:
                node_to_counts[node_id] = {"perfect": 0, "not_perfect": 0}
            node_to_counts[node_id]["perfect"] += counts["perfect"]
            node_to_counts[node_id]["not_perfect"] += counts["not_perfect"]

        if total_count >= milestone:
            elapsed = time.perf_counter() - start_time
            print("\nMilestone reached:")
            print(f"  Total reads processed: {total_count}")
            print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
            print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
            print(f"  Elapsed time: {elapsed:.2f} seconds.")
            milestone += args.milestone

    # Optional: enrich with GFA lengths
    if args.gfa:
        print(f"Parsing GFA for lengths: {args.gfa}")
        nid2len = parse_gfa_lengths(args.gfa)
        merge_lengths(node_to_counts, nid2len)

    output_file = args.output
    if not output_file:
        output_file = "reads_by_node.json" if args.output_format == "json" else "reads_by_node.pkl"

    print(f"\nSaving output to {output_file} as {args.output_format.upper()}...")

    if args.output_format == "json":
        # Convert int keys to strings for JSON compatibility
        json_data = {str(k): v for k, v in node_to_counts.items()}
        with open(output_file, "w") as f:
            json.dump(json_data, f, indent=2)
    else:
        with open(output_file, "wb") as f:
            pickle.dump(node_to_counts, f, protocol=pickle.HIGHEST_PROTOCOL)

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_count}")
    print(f"  Perfect read segments: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
    print(f"  Not-perfect reads segments: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
    print(f"  Nodes included: {len(node_to_counts)}")
    if args.gfa:
        n_with_len = sum(1 for v in node_to_counts.values() if v.get("length", 0) > 0)
        print(f"  Nodes with length>0: {n_with_len}")
    print(f"  Elapsed time: {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()

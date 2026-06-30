#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
import time
import concurrent.futures

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import vg_pb2
import json
import pickle
from typing import Dict


# ─────────────────────────────────────────────────────────────────────────────
# Low-level helpers

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


# ─────────────────────────────────────────────────────────────────────────────
# Core per-group processing

def process_group_serialized(args):
    group_number, messages, chrom_name = args
    total = 0
    perfect = 0
    not_perfect = 0
    # node_id -> {"perfect": X, "not_perfect": Y, "max_read_length": R, "max_cigar_length": C}
    node_counts: Dict[int, Dict[str, int]] = {}

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

            # Compute segment read-length on this node (sum of to_length), and CIGAR string
            seg_to_len_sum = 0
            cigar_parts = []
            for edit in mapping.edit:
                from_len = edit.from_length
                to_len = edit.to_length
                seq_len = len(edit.sequence)

                if from_len == to_len:
                    if seq_len == 0:
                        cigar_parts.append(f"{from_len}M")
                    else:
                        cigar_parts.append(f"{from_len}X")
                        is_perfect = False
                elif from_len > 0 and to_len == 0:
                    cigar_parts.append(f"{from_len}D")
                    is_perfect = False
                elif from_len == 0 and to_len > 0:
                    cigar_parts.append(f"{to_len}I")
                    is_perfect = False
                else:
                    # Any other case is non-perfect
                    is_perfect = False

                seg_to_len_sum += to_len

            cigar_len = len("".join(cigar_parts))  # count chars in CIGAR

            rec = node_counts.get(node_id)
            if rec is None:
                rec = {
                    "perfect": 0,
                    "not_perfect": 0,
                    "max_read_length": 0,
                    "max_cigar_length": 0,
                }
                node_counts[node_id] = rec

            if is_perfect:
                rec["perfect"] += 1
                perfect += 1
            else:
                rec["not_perfect"] += 1
                not_perfect += 1

            # Update maxima
            if seg_to_len_sum > rec["max_read_length"]:
                rec["max_read_length"] = seg_to_len_sum
            if cigar_len > rec["max_cigar_length"]:
                rec["max_cigar_length"] = cigar_len

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
# Main

def main():
    parser = argparse.ArgumentParser(
        description="Parse a GAM file and count perfect/non-perfect mappings per node; "
                    "also record per-node max_read_length (sum of to_length) and max_cigar_length."
    )
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=4, help="Number of worker processes (default: 4)")
    parser.add_argument("--max_pending", type=int, default=16, help="Max number of groups in flight (default: 16)")
    parser.add_argument("--chr", default="", help="Chromosome name to filter on (default: \"\")")
    parser.add_argument("--milestone", type=int, default=100_000_000, help="Number of reads between progress updates")
    parser.add_argument("--output_format", choices=["json", "pickle"], default="json", help="Output format (default: json)")
    parser.add_argument("--output", help="Output file name (overrides default)")
    args = parser.parse_args()
    start_time = time.perf_counter()

    total_count = 0
    perfect_count = 0
    not_perfect_count = 0
    milestone = args.milestone
    # node_id -> {"perfect": int, "not_perfect": int, "max_read_length": int, "max_cigar_length": int}
    node_to_counts: Dict[int, Dict[str, int]] = {}

    for group_number, node_reads, t, p, np in process_groups_pipeline(
            args.filename, args.threads, args.max_pending, chrom_name=args.chr):

        total_count += t
        perfect_count += p
        not_perfect_count += np

        # Merge per-node maps
        for node_id, rec in node_reads.items():
            out = node_to_counts.get(node_id)
            if out is None:
                node_to_counts[node_id] = rec
            else:
                out["perfect"] += rec["perfect"]
                out["not_perfect"] += rec["not_perfect"]
                if rec["max_read_length"] > out["max_read_length"]:
                    out["max_read_length"] = rec["max_read_length"]
                if rec["max_cigar_length"] > out["max_cigar_length"]:
                    out["max_cigar_length"] = rec["max_cigar_length"]

        if total_count >= milestone:
            elapsed = time.perf_counter() - start_time
            print("\nMilestone reached:")
            print(f"  Total reads processed: {total_count}")
            print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
            print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
            print(f"  Elapsed time: {elapsed:.2f} seconds.")
            milestone += args.milestone

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
    print(f"  Elapsed time: {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()

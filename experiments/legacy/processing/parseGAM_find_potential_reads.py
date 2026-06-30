#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import sys
import concurrent.futures
import vg_pb2
import json
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
    node_to_reads = {}

    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)

        if chrom_name and not is_on_chromosome(alignment, chrom_name):
            continue

        total += 1

        if alignment.identity == 1.0:
            perfect += 1
            continue

        not_perfect += 1
        aln_dict = MessageToDict(alignment)

        for mapping in alignment.path.mapping:
            if mapping.position.node_id:
                node_id = mapping.position.node_id
                if node_id not in node_to_reads:
                    node_to_reads[node_id] = []
                node_to_reads[node_id].append(aln_dict)

    return group_number, node_to_reads, total, perfect, not_perfect


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


def main():
    parser = argparse.ArgumentParser(description="Parse a GAM file, group by node ID, and filter non-perfect alignments.")
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=4, help="Number of worker processes (default: 4)")
    parser.add_argument("--max_pending", type=int, default=16, help="Max number of groups in flight (default: 16)")
    parser.add_argument("--chr", default="", help="Chromosome name to filter on (default: \"\")")
    parser.add_argument("--output_json", default="reads_by_node.json", help="Output JSON with node_id -> reads")
    parser.add_argument("--milestone", type=int, default=100_000_000,
                        help="Number of reads between progress updates (default: 100000000)")

    args = parser.parse_args()

    start_time = time.perf_counter()

    total_count = 0
    perfect_count = 0
    not_perfect_count = 0
    milestone = args.milestone
    node_to_reads_all = {}

    for group_number, node_reads, t, p, np in process_groups_pipeline(
            args.filename, args.threads, args.max_pending, chrom_name=args.chr):

        total_count += t
        perfect_count += p
        not_perfect_count += np

        for node_id, reads in node_reads.items():
            if node_id not in node_to_reads_all:
                node_to_reads_all[node_id] = []
            node_to_reads_all[node_id].extend(reads)

        if total_count >= milestone:
            elapsed = time.perf_counter() - start_time
            print("\nEarly Stop Summary:")
            print(f"  Total reads processed: {total_count}")
            print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
            print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
            print(f"Elapsed time: {elapsed:.2f} seconds.")
            milestone += args.milestone

    with open(args.output_json, "w") as out_f:
        json.dump(node_to_reads_all, out_f)

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_count}")
    print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
    print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
    print(f"Elapsed time: {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()
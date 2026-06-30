#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import sys
import concurrent.futures
import vg_pb2  # Import the generated protobuf module


def read_varint(stream):
    """Read a varint from the stream."""
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
    """Check if the file is gzipped by reading its magic number."""
    with open(filename, 'rb') as f:
        magic = f.read(2)
    return magic == b'\x1f\x8b'


def parse_gam_file_groups(filename, expected_tag="GAM"):
    """
    Generator that traverses the GAM file group by group.

    Each group starts with a varint count and a type tag. If the tag matches
    the expected_tag (default "GAM"), the generator yields a tuple:
    (group_number, list of raw message bytes for each alignment in the group).
    """
    if is_gzipped(filename):
        f = gzip.open(filename, 'rb')
    else:
        f = open(filename, 'rb')

    group_number = 0
    try:
        while True:
            try:
                group_count = read_varint(f)
            except EOFError:
                break  # End of file reached
            group_number += 1
            print(f"Reading Group {group_number}: {group_count} messages")

            if group_count == 0:
                continue

            try:
                tag_size = read_varint(f)
            except EOFError:
                print("Unexpected EOF when reading type tag size.")
                break
            tag_bytes = f.read(tag_size)
            if len(tag_bytes) != tag_size:
                print("Unexpected EOF when reading type tag bytes.")
                break
            tag = tag_bytes.decode("utf-8")
            print(f"  Type tag: {tag}")

            if tag != expected_tag:
                print("  Skipping group with unexpected tag.")
                # Skip remaining messages in this group.
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        break
                    f.seek(msg_size, os.SEEK_CUR)
                continue
            else:
                messages = []
                for _ in range(group_count - 1):
                    try:
                        msg_size = read_varint(f)
                    except EOFError:
                        print("Unexpected EOF when reading message size.")
                        break
                    msg_bytes = f.read(msg_size)
                    if len(msg_bytes) != msg_size:
                        print("Unexpected EOF when reading message bytes.")
                        break
                    messages.append(msg_bytes)
                yield (group_number, messages)
    finally:
        f.close()


def process_group(group_tuple):
    """
    Parse each raw message in the group into a vg_pb2.Alignment.
    Returns a tuple (group_number, list of parsed alignments).
    """
    group_number, messages = group_tuple
    alignments = []
    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)
        alignments.append(alignment)
    return (group_number, alignments)


def process_groups_pipeline(filename, threads, max_pending=10):
    """
    Read groups from the GAM file and process them concurrently as they are read.
    The max_pending parameter limits the number of groups in flight.
    Yields processed group results as soon as they are available.
    """
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        pending_futures = []
        # Iterate over groups one at a time
        for group in parse_gam_file_groups(filename):
            future = executor.submit(process_group, group)
            pending_futures.append(future)
            # If too many futures are pending, wait for at least one to finish.
            if len(pending_futures) >= max_pending:
                done, not_done = concurrent.futures.wait(
                    pending_futures, return_when=concurrent.futures.FIRST_COMPLETED)
                for completed in done:
                    yield completed.result()
                pending_futures = list(not_done)  # Convert set back to list
        # Yield any remaining futures
        for future in concurrent.futures.as_completed(pending_futures):
            yield future.result()


def main():
    parser = argparse.ArgumentParser(
        description="Multi-threaded parsing of a GAM file by group using vg_pb2.")
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of worker threads to use (default: 4)")
    parser.add_argument("--max_pending", type=int, default=10,
                        help="Maximum number of groups to process concurrently (default: 10)")
    args = parser.parse_args()

    start_time = time.perf_counter()
    for group_number, alignments in process_groups_pipeline(args.filename, args.threads, args.max_pending):
        print(f"Processed Group {group_number}: Parsed {len(alignments)} alignments")
        for alignment in alignments:
            print(f"  Alignment name: {alignment.name}")
    end_time = time.perf_counter()
    print(f"Elapsed time: {end_time - start_time:.6f} seconds")


if __name__ == "__main__":
    main()

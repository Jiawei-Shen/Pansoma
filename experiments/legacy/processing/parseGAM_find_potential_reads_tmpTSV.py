#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import vg_pb2
import sqlite3
import concurrent.futures


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
    records = []

    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)

        if chrom_name and not is_on_chromosome(alignment, chrom_name):
            continue

        total += 1
        # if alignment.identity == 1.0:
        #     perfect += 1
        #     continue
        #
        # not_perfect += 1
        read_seq = alignment.sequence
        read_qual = alignment.quality.decode('latin1')
        read_offset = 0
        mapping_quality = alignment.mapping_quality

        for mapping in alignment.path.mapping:
            if not mapping.position.node_id:
                continue

            node_id = mapping.position.node_id
            offset = mapping.position.offset if mapping.HasField("position") else 0

            node_seq = ""
            node_qual = ""

            for edit in mapping.edit:
                aligned_len = edit.from_length

                if aligned_len > 0:
                    seg_seq = read_seq[read_offset:read_offset + aligned_len]
                    seg_qual = read_qual[read_offset:read_offset + aligned_len]
                    if edit.sequence:
                        seg_seq = seg_seq.lower()
                    node_seq += seg_seq
                    node_qual += seg_qual
                    read_offset += aligned_len

                elif edit.sequence:
                    insert_len = len(edit.sequence)
                    seg_seq = read_seq[read_offset:read_offset + insert_len].lower()
                    seg_qual = read_qual[read_offset:read_offset + insert_len]
                    node_seq += seg_seq
                    node_qual += seg_qual
                    read_offset += insert_len

            records.append((
                node_id,
                offset,
                node_seq,
                node_qual,
                mapping_quality
            ))

    return group_number, records, total#, perfect, not_perfect

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


def setup_sqlite_db(db_path):
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute("PRAGMA temp_store = MEMORY")

    conn.execute("""
        CREATE TABLE IF NOT EXISTS alignments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            node_id INTEGER,
            offset INTEGER,
            sequence TEXT,
            quality TEXT,
            mapping_quality INTEGER
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_node_id ON alignments(node_id)")
    return conn


def main():
    parser = argparse.ArgumentParser(description="Parse a GAM file and create a queryable SQLite database indexed by node_id.")
    parser.add_argument("filename", help="Path to the GAM file")
    parser.add_argument("--threads", type=int, default=8, help="Number of worker processes")
    parser.add_argument("--max_pending", type=int, default=16, help="Max groups in parallel")
    parser.add_argument("--chr", default="", help="Chromosome name to filter on")
    parser.add_argument("--output_db", default="./tmp/pileup.sqlite", help="Output SQLite filename")
    parser.add_argument("--tmp_dir", default='./tmp', help="Temp directory (unused here, kept for compatibility)")
    parser.add_argument("--milestone", type=int, default=10_000_000, help="Progress update interval")

    args = parser.parse_args()
    os.makedirs(args.tmp_dir, exist_ok=True)

    db_path = os.path.join(".", args.output_db)
    conn = setup_sqlite_db(db_path)
    cursor = conn.cursor()

    start_time = time.perf_counter()
    total_count = 0
    # perfect_count = 0
    # not_perfect_count = 0
    milestone = args.milestone

    for group_number, records, t in process_groups_pipeline(
        args.filename, args.threads, args.max_pending, chrom_name=args.chr):

        total_count += t

        cursor.executemany("""
            INSERT INTO alignments (node_id, offset, sequence, quality, mapping_quality)
            VALUES (?, ?, ?, ?, ?)
        """, records)

        if total_count >= milestone:
            elapsed = time.perf_counter() - start_time
            print(f"\nReached {total_count} reads.")
            # print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}%)")
            # print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}%)")
            print(f"  Elapsed time: {elapsed:.2f} seconds.")
            milestone += args.milestone

    conn.commit()
    conn.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_count}")
    # print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}%)")
    # print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}%)")
    print(f"  Elapsed time: {elapsed:.2f} seconds.")
    print(f"\nSQLite database saved to: {db_path}")


if __name__ == "__main__":
    main()
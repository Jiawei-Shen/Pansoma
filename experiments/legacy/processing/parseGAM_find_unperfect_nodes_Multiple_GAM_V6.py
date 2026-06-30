#!/usr/bin/env python3
import argparse
import gzip
import os
import time
import concurrent.futures
import vg_pb2
import json
import pickle
from typing import Dict, Iterable, List, Tuple, Optional


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
    """
    Yield (group_number, messages_bytes_list) for a GAM file.

    Each 'group' starts with a varint group_count, followed by a 'GAM' tag
    and group_count-1 Alignment messages (varint size + bytes).
    """
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

            # Skip non-GAM groups by seeking past their payloads
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
    """
    Worker function for a single group of serialized Alignment messages.
    Returns: (group_number, node_counts, total, perfect, not_perfect)
    """
    group_number, messages, chrom_name, min_mapq = args
    total = 0
    perfect = 0
    not_perfect = 0
    # node_id -> {"perfect": X, "not_perfect": Y, "max_read_length": R, "max_cigar_length": C}
    node_counts: Dict[int, Dict[str, int]] = {}

    for msg_bytes in messages:
        alignment = vg_pb2.Alignment()
        alignment.ParseFromString(msg_bytes)

        # Skip low mapping quality reads
        if alignment.mapping_quality < min_mapq:
            continue

        if chrom_name and not is_on_chromosome(alignment, chrom_name):
            continue

        total += 1

        for mapping in alignment.path.mapping:
            if not mapping.position.node_id:
                continue

            node_id = mapping.position.node_id
            is_perfect = True

            # Compute segment read-length on this node (sum of to_length), and CIGAR string length
            seg_to_len_sum = 0
            cigar_len = 0
            for edit in mapping.edit:
                from_len = edit.from_length
                to_len = edit.to_length
                seq_len = len(edit.sequence)

                if from_len == to_len:
                    if seq_len == 0:
                        # e.g. 10M -> adds len("10M") to cigar length
                        cigar_len += len(f"{from_len}M")
                    else:
                        cigar_len += len(f"{from_len}X")
                        is_perfect = False
                elif from_len > 0 and to_len == 0:
                    cigar_len += len(f"{from_len}D")
                    is_perfect = False
                elif from_len == 0 and to_len > 0:
                    cigar_len += len(f"{to_len}I")
                    is_perfect = False
                else:
                    # Any other case is non-perfect
                    is_perfect = False
                    # still account for cigar text conservatively
                    # (use a generic op length to contribute to max length)
                    op_len = to_len if to_len > 0 else from_len
                    cigar_len += len(str(op_len)) + 1  # e.g. "5?" -> 2

                seg_to_len_sum += to_len

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


def process_groups_across_files(
    filenames: List[str],
    threads: int,
    max_pending: int,
    chrom_name: str,
    min_mapq: int,
):
    """
    Single ProcessPool shared across all input files.
    Yields tuples from completed futures as they finish.
    """
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        pending_futures = []
        for fname in filenames:
            for group_number, messages in parse_gam_file_groups(fname):
                args = (group_number, messages, chrom_name, min_mapq)
                pending_futures.append(executor.submit(process_group_serialized, args))

                # Backpressure to avoid unbounded memory
                if len(pending_futures) >= max_pending:
                    done, not_done = concurrent.futures.wait(
                        pending_futures,
                        return_when=concurrent.futures.FIRST_COMPLETED
                    )
                    for completed in done:
                        yield completed.result()
                    pending_futures = list(not_done)

        # Drain remaining
        for future in concurrent.futures.as_completed(pending_futures):
            yield future.result()


# ─────────────────────────────────────────────────────────────────────────────
# Argument helpers

def collect_input_files(positional: List[str], from_list: Optional[str]) -> List[str]:
    files = []
    if positional:
        files.extend(positional)
    if from_list:
        with open(from_list, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    files.append(line)
    # De-dup while preserving order
    seen = set()
    uniq = []
    for p in files:
        if p not in seen:
            uniq.append(p)
            seen.add(p)
    if not uniq:
        raise SystemExit("No input GAM files provided.")
    return uniq


# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parse one or more GAM files and count perfect/non-perfect mappings per node;\n"
            "also record per-node max_read_length (sum of to_length) and max_cigar_length.\n"
            "Outputs a single merged result across all inputs."
        )
    )
    parser.add_argument("input_gams", nargs="*", help="Path(s) to GAM or GAM.GZ files")
    parser.add_argument("--from_list", help="Text file with one GAM path per line (merged with positional files)")
    parser.add_argument("--threads", type=int, default=4, help="Number of worker processes (default: 4)")
    parser.add_argument("--max_pending", type=int, default=16, help="Max number of groups in flight (default: 16)")
    parser.add_argument("--chr", default="", help='Chromosome name to filter on (default: "")')
    parser.add_argument("--min_mapq", type=int, default=6, help="Minimum MAPQ to keep (default: 6; old code skipped <=5)")
    parser.add_argument("--milestone", type=int, default=100_000_000, help="Reads between progress updates (global)")
    parser.add_argument("--output_format", choices=["json", "pickle"], default="json", help="Output format (default: json)")
    parser.add_argument("--output", help="Output file name (overrides default)")
    parser.add_argument("--per_file_stats", action="store_true", help="Print per-file stats (approximate)")
    args = parser.parse_args()

    filenames = collect_input_files(args.input_gams, args.from_list)

    start_time = time.perf_counter()

    total_count = 0
    perfect_count = 0
    not_perfect_count = 0
    milestone_next = args.milestone
    # node_id -> {"perfect": int, "not_perfect": int, "max_read_length": int, "max_cigar_length": int}
    node_to_counts: Dict[int, Dict[str, int]] = {}

    # Optional: rough per-file counters (based on group processing order)
    per_file_stats = {}
    if args.per_file_stats:
        # map of index -> filename, and running tally when we switch files
        per_file_stats = {fn: {"total": 0, "perfect": 0, "not_perfect": 0} for fn in filenames}

    # We don't have explicit file boundaries in the executor results. If per-file
    # stats are important, consider running files sequentially. For now, we only
    # report a global summary by default. If --per_file_stats is requested, we
    # compute approximate per-file stats by processing files sequentially here.
    # To ensure correctness, we’ll process files sequentially if per_file_stats=True.
    if args.per_file_stats:
        for fn in filenames:
            local_total = local_perfect = local_not_perfect = 0
            for group_number, messages in parse_gam_file_groups(fn):
                res = process_group_serialized((group_number, messages, args.chr, args.min_mapq))
                _, node_reads, t, p, np = res

                total_count += t
                perfect_count += p
                not_perfect_count += np
                local_total += t
                local_perfect += p
                local_not_perfect += np

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

                if total_count >= milestone_next:
                    elapsed = time.perf_counter() - start_time
                    print("\nMilestone reached:")
                    print(f"  Total reads processed: {total_count}")
                    print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
                    print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
                    print(f"  Elapsed time: {elapsed:.2f} seconds.")
                    milestone_next += args.milestone

            per_file_stats[fn]["total"] = local_total
            per_file_stats[fn]["perfect"] = local_perfect
            per_file_stats[fn]["not_perfect"] = local_not_perfect
    else:
        # Fast path: use a shared pool across all files
        for group_number, node_reads, t, p, np in process_groups_across_files(
            filenames, args.threads, args.max_pending, chrom_name=args.chr, min_mapq=args.min_mapq
        ):
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

            if total_count >= milestone_next:
                elapsed = time.perf_counter() - start_time
                print("\nMilestone reached:")
                print(f"  Total reads processed: {total_count}")
                print(f"  Perfect reads: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
                print(f"  Not-perfect reads: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
                print(f"  Elapsed time: {elapsed:.2f} seconds.")
                milestone_next += args.milestone

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
    print("\nFinal Summary (merged across all inputs):")
    print(f"  Total reads processed: {total_count}")
    if total_count > 0:
        print(f"  Perfect read segments: {perfect_count} ({(perfect_count / total_count * 100):.2f}% of total)")
        print(f"  Not-perfect read segments: {not_perfect_count} ({(not_perfect_count / total_count * 100):.2f}% of total)")
    else:
        print("  Perfect read segments: 0 (0.00% of total)")
        print("  Not-perfect read segments: 0 (0.00% of total)")
    print(f"  Nodes included: {len(node_to_counts)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds.")

    if args.per_file_stats:
        print("\nPer-file summary (sequential mode for exact tallies):")
        for fn, st in per_file_stats.items():
            t = st["total"]
            p = st["perfect"]
            npv = st["not_perfect"]
            pct_p = (p / t * 100) if t else 0.0
            pct_np = (npv / t * 100) if t else 0.0
            print(f"  {fn}")
            print(f"    Total: {t}")
            print(f"    Perfect: {p} ({pct_p:.2f}%)")
            print(f"    Not-perfect: {npv} ({pct_np:.2f}%)")


if __name__ == "__main__":
    main()

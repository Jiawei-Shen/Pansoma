#!/usr/bin/env python3
import argparse
import gzip
import json
import pickle
import time
import gc
from collections import defaultdict
import vg_pb2

# ─────────────────────────────────────────────────────────────────────────────
# Lightweight segment container using __slots__ to cut per-object overhead
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'rq', 'strand')
    def __init__(self, offset, seq, bq, rq, strand):
        self.offset = offset
        self.seq    = seq
        self.bq     = bq
        self.rq     = rq
        self.strand = strand

# ─────────────────────────────────────────────────────────────────────────────
# Utility: Read varint from stream
def read_varint(stream):
    value = 0
    shift_amount = 0
    while True:
        byte_pair = stream.read(1)
        if not byte_pair:
            raise EOFError("EOF while reading varint")
        byte_value = byte_pair[0]
        value |= (byte_value & 0x7F) << shift_amount
        if not (byte_value & 0x80):
            return value
        shift_amount += 7

# Utility: Check if a file is gzipped
def file_is_gzip(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"

# Generator: Yield one GAM record at a time
def gam_record_iter(path, tag="GAM"):
    open_func = gzip.open if file_is_gzip(path) else open
    with open_func(path, "rb") as f:
        while True:
            try:
                group_count = read_varint(f)
            except EOFError:
                break
            if group_count == 0:
                continue
            try:
                tag_len = read_varint(f)
                group_tag = f.read(tag_len).decode()
            except (EOFError, UnicodeDecodeError):
                break
            if group_tag != tag:
                # skip the rest of this group
                for _ in range(group_count - 1):
                    skip_len = read_varint(f)
                    f.seek(skip_len, 1)
                continue
            # yield each message in the group
            for _ in range(group_count - 1):
                try:
                    msg_size = read_varint(f)
                    yield f.read(msg_size)
                except EOFError:
                    break

# ─────────────────────────────────────────────────────────────────────────────
# Process one alignment record, reusing a single Alignment instance
def process_alignment(raw_message, wanted_nodes, chrom_filter):
    segment_dict = {}
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw_message)

    # optional chromosome filter
    if chrom_filter and not any(pos.name == chrom_filter for pos in alignment.refpos):
        return segment_dict, 0

    read_sequence = alignment.sequence
    read_quality_bytes = alignment.quality
    mapping_quality = alignment.mapping_quality
    read_offset = 0
    reads = 1

    for mapping in alignment.path.mapping:
        node_id = mapping.position.node_id

        if node_id not in wanted_nodes:
            # advance offset past this mapping
            for edit in mapping.edit:
                read_offset += max(edit.from_length, len(edit.sequence))
            continue

        node_offset = mapping.position.offset
        strand_char = "-" if mapping.position.is_reverse else "+"
        sequence_parts = []
        quality_bytes = bytearray()

        # collect all matched/inserted fragments
        for edit in mapping.edit:
            if edit.from_length:
                frag = read_sequence[read_offset: read_offset + edit.from_length]
                sequence_parts.append(frag.lower() if edit.sequence else frag)
                quality_bytes.extend(read_quality_bytes[read_offset: read_offset + edit.from_length])
                read_offset += edit.from_length
            elif edit.sequence:
                ins_len = len(edit.sequence)
                frag = read_sequence[read_offset: read_offset + ins_len]
                sequence_parts.append(frag.lower())
                quality_bytes.extend(read_quality_bytes[read_offset: read_offset + ins_len])
                read_offset += ins_len

        # store as Slot-based object
        seg = Segment(
            node_offset,
            "".join(sequence_parts),
            bytes(quality_bytes),
            mapping_quality,
            strand_char
        )
        segment_dict.setdefault(node_id, []).append(seg)

    return segment_dict, reads

# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
def run_pipeline(gam_path, stats_path, output_prefix, output_format, milestone_step, chrom_filter):
    print(f"Loading stats: {stats_path}")
    with open(stats_path, "rb") as stats_file:
        stats_data = pickle.load(stats_file)

    # filter wanted nodes
    wanted_nodes = {
        int(node_id)
        for node_id, stat in stats_data.items()
        if stat["not_perfect"] > 1 and stat["not_perfect"] / (stat["perfect"] + stat["not_perfect"]) > 0.10
    }
    del stats_data
    gc.collect()

    print(f"Filtered {len(wanted_nodes)} nodes from stats file.")

    merged_segments = defaultdict(list)
    total_reads = 0
    next_milestone = milestone_step
    start_time = time.perf_counter()

    # # reuse one Alignment object
    # alignment = vg_pb2.Alignment()

    for raw_msg in gam_record_iter(gam_path):
        segment_dict, read_count = process_alignment(raw_msg, wanted_nodes, chrom_filter)
        total_reads += read_count

        for node_id, segs in segment_dict.items():
            merged_segments[node_id].extend(segs)

        if total_reads >= next_milestone:
            elapsed = time.perf_counter() - start_time
            print(f"{total_reads} reads processed | {elapsed:.1f} seconds")
            next_milestone += milestone_step

    # finalize in-memory dict
    merged_segments = dict(merged_segments)

    # write outputs
    if output_format in ("pkl", "both"):
        pickle_path = output_prefix + ".pkl"
        with open(pickle_path, "wb") as pkl_f:
            pickle.dump(merged_segments, pkl.HIGHEST_PROTOCOL)
        print(f"Wrote Pickle: {pickle_path}")

    if output_format in ("json", "both"):
        json_path = output_prefix + ".json"
        with open(json_path, "w") as j_f:
            json.dump({
                str(k): [
                    {
                        "offset": seg.offset,
                        "seq":    seg.seq,
                        "bq":     seg.bq.hex(),
                        "rq":     seg.rq,
                        "strand": seg.strand
                    } for seg in v
                ] for k, v in merged_segments.items()
            }, j_f)
        print(f"Wrote JSON: {json_path}")

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads}")
    print(f"  Nodes included: {len(merged_segments)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")

# ─────────────────────────────────────────────────────────────────────────────
# Command-line interface
def main():
    parser = argparse.ArgumentParser(description="Extract read-segments by node from a GAM file.")
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the node stats pickle file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--fmt", choices=["json", "pkl", "both"], default="json", help="Output format")
    parser.add_argument("--milestone", type=int, default=10_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        output_format=args.fmt,
        milestone_step=args.milestone,
        chrom_filter=args.chr
    )

if __name__ == "__main__":
    main()

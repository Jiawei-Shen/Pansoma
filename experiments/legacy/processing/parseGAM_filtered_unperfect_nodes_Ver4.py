#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
from collections import defaultdict
import vg_pb2


# ─────────────────────────────────────────────────────────────────────────────
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'cigar', 'rq', 'strand')
    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset = offset
        self.seq    = seq
        self.bq     = bq
        self.cigar  = cigar
        self.rq     = rq
        self.strand = strand


RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size


# ─────────────────────────────────────────────────────────────────────────────
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

def file_is_gzip(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"

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
                for _ in range(group_count - 1):
                    skip_len = read_varint(f)
                    f.seek(skip_len, 1)
                continue
            for _ in range(group_count - 1):
                try:
                    msg_size = read_varint(f)
                    yield f.read(msg_size)
                except EOFError:
                    break

# ─────────────────────────────────────────────────────────────────────────────
def build_cigar(mapping_edits):
    cigar_parts = []
    for edit in mapping_edits:
        from_len = edit.from_length
        to_len = edit.to_length
        edit_len = len(edit.sequence)

        if from_len == to_len:
            if edit_len == 0:
                cigar_parts.append(f"{from_len}M")
            else:
                cigar_parts.append(f"{from_len}X")
        elif from_len > 0 and to_len == 0:
            cigar_parts.append(f"{from_len}D")
        elif from_len == 0 and to_len > 0:
            cigar_parts.append(f"{to_len}I")
        else:
            raise ValueError(f"Unexpected edit: from_length={from_len}, sequence_length={to_len}")
    return "".join(cigar_parts)

# ─────────────────────────────────────────────────────────────────────────────
def process_alignment(raw_message, wanted_nodes, chrom_filter):
    segment_dict = {}
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw_message)

    # Filter out low mapping quality
    if alignment.mapping_quality <= 10:
        return segment_dict

    if chrom_filter and not any(pos.name == chrom_filter for pos in alignment.refpos):
        return segment_dict

    read_sequence = alignment.sequence
    read_quality = alignment.quality
    mapping_quality = alignment.mapping_quality
    read_offset = 0

    for mapping in alignment.path.mapping:
        node_id = mapping.position.node_id

        if node_id not in wanted_nodes:
            for edit in mapping.edit:
                # read_offset += max(edit.from_length, len(edit.sequence))
                read_offset += edit.to_length
            continue

        node_offset = mapping.position.offset
        strand_char = b"-" if mapping.position.is_reverse else b"+"

        sequence_parts = []
        quality_parts = bytearray()
        cigar_parts = []

        for edit in mapping.edit:
            from_len = edit.from_length
            to_len = edit.to_length
            edit_len = len(edit.sequence)

            # Build CIGAR part
            if from_len == to_len:
                if edit_len == 0:
                    cigar_parts.append(f"{from_len}M")
                else:
                    cigar_parts.append(f"{from_len}X")
            elif from_len > 0 and to_len == 0:
                cigar_parts.append(f"{from_len}D")
            elif from_len == 0 and to_len > 0:
                cigar_parts.append(f"{to_len}I")
            else:
                raise ValueError(f"Unexpected edit: from_length={from_len}, sequence_length={edit_len}")

            # Append sequence and quality
            # edit_length = max(from_len, edit_len)
            edit_length = to_len
            sequence_fragment = read_sequence[read_offset: read_offset + edit_length]
            quality_fragment = read_quality[read_offset: read_offset + edit_length]

            sequence_parts.append(sequence_fragment.upper())
            quality_parts.extend(quality_fragment)
            read_offset += edit_length

        # Construct final padded fields
        cigar_string = "".join(cigar_parts)
        seq_final = "".join(sequence_parts).encode().ljust(150, b'\x00')[:150]
        bq_final = bytes(quality_parts).ljust(150, b'\x00')[:150]
        cigar_bytes = cigar_string.encode().ljust(20, b'\x00')[:20]

        seg = Segment(
            offset=node_offset,
            seq=seq_final,
            bq=bq_final,
            cigar=cigar_bytes,
            rq=mapping_quality,
            strand=strand_char
        )
        segment_dict.setdefault(node_id, []).append(seg)

    return segment_dict

# ─────────────────────────────────────────────────────────────────────────────
def initialize_output_files(stats_path, output_prefix):
    with open(stats_path, "rb") as stats_file:
        stats_data = pickle.load(stats_file)

    block_infos = {}
    wanted_nodes = set()

    current_offset = 0
    total_nodes = 0

    for node_id_str, stat in stats_data.items():
        total_nodes += 1
        node_id = int(node_id_str)
        perfect = stat["perfect"]
        not_perfect = stat["not_perfect"]
        if not_perfect > 1 and not_perfect / (perfect + not_perfect) > 0.10:
            wanted_nodes.add(node_id)
            n_records = perfect + not_perfect

            block_infos[node_id] = {
                "offset": current_offset,
                "n_records": n_records,
                "current_pos": 0
            }
            current_offset += 4 + 4 + 2 + n_records * RECORD_SIZE

    print(f"Filtered {len(wanted_nodes)} nodes from {total_nodes} total nodes "
          f"({len(wanted_nodes) / total_nodes:.2%} selected).")
    del stats_data
    gc.collect()

    dat_path = output_prefix + ".dat"

    with open(dat_path, "wb") as f:
        f.write(b"MYFMT\1")
        f.write(struct.pack("<BBI16s", 0, 1, len(block_infos), b'\x00' * 16))

        blank_record = RECORD_STRUCT.pack(0, b'\x00'*150, b'\x00'*150, b'\x00'*20, 0, b'+')
        for node_id, info in block_infos.items():
            block_header = struct.pack("<I I H", node_id, info["n_records"], 0)
            f.write(block_header + blank_record * info["n_records"])

    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx_file:
        idx_file.write(struct.pack("<I", len(block_infos)))
        for node_id, info in block_infos.items():
            idx_file.write(struct.pack("<I Q I I H", node_id, info["offset"],
                                       4 + 4 + 2 + info["n_records"] * RECORD_SIZE, info["n_records"], 0))

    return block_infos, dat_path, wanted_nodes

# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter):
    print(f"Initializing output files...")
    block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)
    print(f"Output file created: {dat_path}")

    BUFFER_SEGMENTS = 200_000_000  # it takes about 240GB memory

    next_milestone = milestone_step
    total_reads = 0
    total_segments = 0
    start_time = time.perf_counter()

    dat_fh = open(dat_path, "r+b")
    segment_buffer = defaultdict(list)

    def flush_segment_buffer():
        nonlocal total_segments
        for node_id, segs in segment_buffer.items():
            if not segs:
                continue
            info = block_infos[node_id]
            base_offset = info["offset"] + 4 + 4 + 2

            batch_blob = bytearray()
            for seg in segs:
                batch_blob += RECORD_STRUCT.pack(seg.offset, seg.seq, seg.bq, seg.cigar, seg.rq, seg.strand)

            pos = base_offset + info["current_pos"] * RECORD_SIZE
            dat_fh.seek(pos)
            dat_fh.write(batch_blob)

            info["current_pos"] += len(segs)

        segment_buffer.clear()
        total_segments = 0

    for raw_msg in gam_record_iter(gam_path):
        segment_dict = process_alignment(raw_msg, wanted_nodes, chrom_filter)
        total_reads += 1

        for node_id, segs in segment_dict.items():
            segment_buffer[node_id].extend(segs)
            total_segments += len(segs)

        if total_segments >= BUFFER_SEGMENTS:
            flush_segment_buffer()

        if total_reads >= next_milestone:
            elapsed = time.perf_counter() - start_time
            print(f"{total_reads} reads processed | {elapsed:.1f} seconds")
            next_milestone += milestone_step

    flush_segment_buffer()
    dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="GAM segment extractor with CIGAR generation.")
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the node stats pickle file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr
    )

if __name__ == "__main__":
    main()

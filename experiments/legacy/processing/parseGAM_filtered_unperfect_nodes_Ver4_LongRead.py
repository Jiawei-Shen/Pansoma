#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
import os
from collections import defaultdict
import vg_pb2


# ─────────────────────────────────────────────────────────────────────────────
# Segment container (stores raw, unpadded bytes for seq/bq; padded at pack time)
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'cigar', 'rq', 'strand')
    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset = offset
        self.seq    = seq      # bytes (unpadded)
        self.bq     = bq       # bytes (unpadded)
        self.cigar  = cigar    # bytes (unpadded; will be padded to 30)
        self.rq     = rq       # int (MAPQ)
        self.strand = strand   # b'+' or b'-'


# ─────────────────────────────────────────────────────────────────────────────
# New variable-size record layout utilities

# Global header: magic + version block
GLOBAL_MAGIC = b"MYFMT\x01"                              # 6 bytes: "MYFMT" + 0x01
GLOBAL_VER_PACK = struct.Struct("<BBI16s")               # major, minor, block_count, reserved[16]
GLOBAL_MAJOR, GLOBAL_MINOR = 0, 2                        # bump minor for new layout
GLOBAL_HEADER_SIZE = len(GLOBAL_MAGIC) + GLOBAL_VER_PACK.size  # 6 + 22 = 28

# Per-node block header now includes node_length
# <I I H I>  -> node_id (u32), n_records (u32), flags (u16=0), node_length (u32)
BLOCK_HDR_PACK = struct.Struct("<I I H I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size  # 14

def make_record_struct(node_length: int) -> struct.Struct:
    """
    Build the struct for a single segment record for a given node_length.
    Layout: <h {L}s {L}s 30s h c
      - i16 offset
      - seq[L] bytes
      - bq[L] bytes
      - cigar[30] bytes (ASCII, null-padded)
      - i16 rq (MAPQ)
      - char strand ('+' / '-')
    """
    return struct.Struct(f"<h{node_length}s{node_length}s30shc")

def record_size(node_length: int) -> int:
    return make_record_struct(node_length).size


# ─────────────────────────────────────────────────────────────────────────────
# GAM parsing

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
# CIGAR builder (unchanged semantics; now stored in 30 bytes)
def build_cigar(mapping_edits):
    cigar_parts = []
    for edit in mapping_edits:
        from_len = edit.from_length
        to_len = edit.to_length
        edit_len = len(edit.sequence)

        if from_len == to_len:
            if edit_len == 0:
                cigar_parts.append(f"{from_len}M")  # match
            else:
                cigar_parts.append(f"{from_len}X")  # substitution
        elif from_len > 0 and to_len == 0:
            cigar_parts.append(f"{from_len}D")      # deletion
        elif from_len == 0 and to_len > 0:
            cigar_parts.append(f"{to_len}I")        # insertion
        else:
            raise ValueError(f"Unexpected edit: from_length={from_len}, to_length={to_len}")
    return "".join(cigar_parts)


# ─────────────────────────────────────────────────────────────────────────────
# Alignment → Segment
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

        # Advance read_offset even if node not wanted, to keep position in read correct
        if node_id not in wanted_nodes:
            for edit in mapping.edit:
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

            # CIGAR part
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

            # Append sequence/quality by to_len
            edit_length = to_len
            sequence_fragment = read_sequence[read_offset: read_offset + edit_length]
            quality_fragment  = read_quality[read_offset: read_offset + edit_length]
            sequence_parts.append(sequence_fragment.upper())
            quality_parts.extend(quality_fragment)
            read_offset += edit_length

        cigar_string = "".join(cigar_parts).encode()
        seq_bytes = "".join(sequence_parts).encode()
        bq_bytes  = bytes(quality_parts)

        seg = Segment(
            offset=node_offset,
            seq=seq_bytes,          # unpadded; pad to node_length on write
            bq=bq_bytes,            # unpadded
            cigar=cigar_string,     # will be padded/truncated to 30
            rq=mapping_quality,
            strand=strand_char
        )
        segment_dict.setdefault(node_id, []).append(seg)

    return segment_dict


# ─────────────────────────────────────────────────────────────────────────────
# Initialize outputs (now variable record sizes and header-offset fixed)
def initialize_output_files(stats_path, output_prefix, default_node_length=150):
    with open(stats_path, "rb") as stats_file:
        stats_data = pickle.load(stats_file)

    block_infos = {}
    wanted_nodes = set()

    # Start after global header to fix the header-offset bug
    current_offset = GLOBAL_HEADER_SIZE
    total_nodes = 0
    warned_default = False

    for node_id_key, stat in stats_data.items():
        total_nodes += 1
        # keys in PKL might be str or int; normalize
        node_id = int(node_id_key)

        perfect = int(stat.get("perfect", 0))
        not_perfect = int(stat.get("not_perfect", 0))
        node_len = int(stat.get("length", 0))

        if node_len <= 0:
            node_len = default_node_length
            if not warned_default:
                print(f"[warn] Missing node length in PKL for some nodes; using default={default_node_length}. "
                      f"Example node_id={node_id}")
                warned_default = True

        # your selection rule
        if (perfect + not_perfect) > 0 and not_perfect > 1 and not_perfect / (perfect + not_perfect) > 0.10:
            wanted_nodes.add(node_id)
            n_records = perfect + not_perfect

            rec_sz = record_size(node_len)
            blk_sz = BLOCK_HDR_SIZE + n_records * rec_sz

            block_infos[node_id] = {
                "offset": current_offset,   # absolute from file start
                "n_records": n_records,
                "current_pos": 0,
                "node_length": node_len,
                "record_size": rec_sz,
                "block_size": blk_sz,
            }
            current_offset += blk_sz

    print(f"Filtered {len(wanted_nodes)} nodes from {total_nodes} total nodes "
          f"({(len(wanted_nodes) / max(total_nodes,1)):.2%} selected).")
    del stats_data
    gc.collect()

    dat_path = output_prefix + ".dat"

    # Write .dat file: global header + blocks
    with open(dat_path, "wb") as f:
        # Global header
        f.write(GLOBAL_MAGIC)
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR, GLOBAL_MINOR, len(block_infos), b'\x00' * 16))

        # Blocks
        for node_id, info in block_infos.items():
            node_len = info["node_length"]
            n_records = info["n_records"]

            # Block header with node_length
            f.write(BLOCK_HDR_PACK.pack(node_id, n_records, 0, node_len))

            # Preallocate blank records for the block
            rec_pack = make_record_struct(node_len)
            blank = rec_pack.pack(
                0,
                b'\x00' * node_len,
                b'\x00' * node_len,
                b'\x00' * 30,
                0,
                b'+'
            )
            for _ in range(n_records):
                f.write(blank)

    # Write .idx: include node_length for convenience
    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx_file:
        # count
        idx_file.write(struct.pack("<I", len(block_infos)))
        # entries: node_id (u32), offset (u64), block_size (u32), n_records (u32), flags (u16), node_length (u32)
        for node_id, info in block_infos.items():
            idx_file.write(struct.pack(
                "<I Q I I H I",
                node_id,
                info["offset"],
                info["block_size"],
                info["n_records"],
                0,
                info["node_length"]
            ))

    return block_infos, dat_path, wanted_nodes


# ─────────────────────────────────────────────────────────────────────────────
# Reuse existing .dat/.idx instead of initializing
def load_existing_output_files(output_prefix):
    idx_path = output_prefix + ".idx"
    dat_path = output_prefix + ".dat"

    if not (os.path.exists(idx_path) and os.path.exists(dat_path)):
        raise FileNotFoundError(f"Expected existing files: {idx_path} and {dat_path}")

    # Parse .idx
    with open(idx_path, "rb") as f:
        raw = f.read(4)
        if len(raw) != 4:
            raise RuntimeError("Corrupt .idx: cannot read block count")
        (count,) = struct.unpack("<I", raw)
        f.seek(0, os.SEEK_END)
        remaining = f.tell() - 4
        if count <= 0 or remaining <= 0:
            raise RuntimeError("Empty or corrupt .idx")
        entry_size = remaining // count
        if entry_size not in (22, 26):  # 22B (old, no node_length) or 26B (new, with node_length)
            entry_size = 26
        f.seek(4)

        entries = []
        for _ in range(count):
            data = f.read(entry_size)
            if len(data) != entry_size:
                raise RuntimeError("Corrupt .idx: truncated entry")
            if entry_size == 26:
                node_id, offset, block_size, n_records, flags, node_len = struct.unpack("<I Q I I H I", data)
            else:
                node_id, offset, block_size, n_records, flags = struct.unpack("<I Q I I H", data)
                node_len = 0
            entries.append((node_id, offset, block_size, n_records, flags, node_len))

    # Verify .dat header & block_count; fill missing node_length from .dat
    with open(dat_path, "rb") as df:
        magic = df.read(len(GLOBAL_MAGIC))
        if magic != GLOBAL_MAGIC:
            raise RuntimeError("Invalid .dat magic/version")
        majors, minors, dat_count, _ = GLOBAL_VER_PACK.unpack(df.read(GLOBAL_VER_PACK.size))
        if dat_count != len(entries):
            print(f"[warn] .dat block_count ({dat_count}) != .idx count ({len(entries)})")

        block_infos = {}
        wanted_nodes = set()
        for node_id, offset, block_size, n_records, flags, node_len in entries:
            if node_len <= 0:
                df.seek(offset, os.SEEK_SET)
                hdr = df.read(BLOCK_HDR_SIZE)
                if len(hdr) != BLOCK_HDR_SIZE:
                    raise RuntimeError(f"Corrupt .dat: cannot read block header at {offset}")
                nid2, nrec2, flg2, node_len = BLOCK_HDR_PACK.unpack(hdr)
                if nid2 != node_id or nrec2 != n_records:
                    print(f"[warn] .dat/.idx mismatch for node {node_id} (idx n={n_records}, dat n={nrec2})")
            rec_sz = record_size(node_len)
            block_infos[node_id] = {
                "offset": offset,
                "n_records": n_records,
                "current_pos": 0,       # NOTE: assumes empty blocks; this does not resume in-place
                "node_length": node_len,
                "record_size": rec_sz,
                "block_size": block_size,
            }
            wanted_nodes.add(node_id)

    print(f"Reusing existing output with {len(block_infos)} node blocks from {idx_path}")
    return block_infos, dat_path, wanted_nodes


# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing):
    if use_existing:
        print("Reusing existing .dat/.idx...")
        block_infos, dat_path, wanted_nodes = load_existing_output_files(output_prefix)
    else:
        print(f"Initializing output files...")
        block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)
        print(f"Output file created: {dat_path}")

    BUFFER_SEGMENTS = 400_000_000  # number of segments buffered before flushing

    next_milestone = milestone_step
    total_reads = 0
    total_segments = 0
    start_time = time.perf_counter()

    dat_fh = open(dat_path, "r+b")
    segment_buffer = defaultdict(list)

    def flush_segment_buffer():
        nonlocal total_segments
        if not segment_buffer:
            return
        for node_id, segs in segment_buffer.items():
            if not segs:
                continue
            info = block_infos[node_id]
            base_offset = info["offset"] + BLOCK_HDR_SIZE  # start of records for this block
            node_len = info["node_length"]
            rec_pack = make_record_struct(node_len)

            # Build contiguous blob for this node batch
            batch = bytearray()
            for seg in segs:
                batch += rec_pack.pack(
                    int(seg.offset),
                    seg.seq.ljust(node_len, b'\x00')[:node_len],
                    seg.bq.ljust(node_len, b'\x00')[:node_len],
                    seg.cigar.ljust(30, b'\x00')[:30],
                    int(seg.rq),
                    seg.strand if seg.strand in (b'+', b'-') else b'+'
                )

            pos = base_offset + info["current_pos"] * info["record_size"]
            dat_fh.seek(pos, os.SEEK_SET)
            dat_fh.write(batch)

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
            segment_buffer = defaultdict(list)

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
    parser = argparse.ArgumentParser(
        description="GAM segment extractor with per-node variable record sizes and node_length from PKL."
    )
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the node stats pickle file (used unless --use-existing)")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    parser.add_argument("--use-existing", action="store_true",
                        help="Reuse existing initialized output (output_prefix.dat/.idx) instead of reinitializing")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr,
        use_existing=args.use_existing
    )

if __name__ == "__main__":
    main()

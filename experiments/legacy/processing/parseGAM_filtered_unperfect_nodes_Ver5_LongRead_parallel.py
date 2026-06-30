#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import vg_pb2


# ─────────────────────────────────────────────────────────────────────────────
# Segment container (raw, unpadded; pad at write time)
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'cigar', 'rq', 'strand')
    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset = offset         # i16 on-disk
        self.seq    = seq            # bytes (unpadded)
        self.bq     = bq             # bytes (unpadded)
        self.cigar  = cigar          # bytes (unpadded; will be padded to 30)
        self.rq     = rq             # i16 on-disk
        self.strand = strand         # b'+' or b'-'


# ─────────────────────────────────────────────────────────────────────────────
# On-disk layout helpers

GLOBAL_MAGIC = b"MYFMT\x01"                          # 6 bytes: "MYFMT" + version byte
GLOBAL_VER_PACK = struct.Struct("<BBI16s")           # major, minor, block_count, reserved[16]
GLOBAL_MAJOR, GLOBAL_MINOR = 0, 2                    # bump minor for new layout
GLOBAL_HEADER_SIZE = len(GLOBAL_MAGIC) + GLOBAL_VER_PACK.size  # 6 + 22 = 28

# Per-node block header (node_length added)
# <I I H I> -> node_id(u32), n_records(u32), flags(u16=0), node_length(u32)
BLOCK_HDR_PACK = struct.Struct("<I I H I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size  # 14

def make_record_struct(node_length: int) -> struct.Struct:
    """
    Per-record struct for given node_length:
      <h {L}s {L}s 30s h c
        i16 offset
        seq[L] bytes
        bq[L] bytes
        cigar[30] bytes (ASCII, null-padded)
        i16 rq
        char strand
    """
    return struct.Struct(f"<h{node_length}s{node_length}s30shc")

def record_size(node_length: int) -> int:
    return make_record_struct(node_length).size


# ─────────────────────────────────────────────────────────────────────────────
# GAM decoding

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
# CIGAR construction (stored in 30B field)
def build_cigar(mapping_edits):
    cigar_parts = []
    for edit in mapping_edits:
        from_len = edit.from_length
        to_len   = edit.to_length
        seq_len  = len(edit.sequence)
        if from_len == to_len:
            cigar_parts.append(f"{from_len}{'M' if seq_len == 0 else 'X'}")
        elif from_len > 0 and to_len == 0:
            cigar_parts.append(f"{from_len}D")
        elif from_len == 0 and to_len > 0:
            cigar_parts.append(f"{to_len}I")
        else:
            raise ValueError(f"Unexpected edit: from_length={from_len}, to_length={to_len}")
    return "".join(cigar_parts)


# ─────────────────────────────────────────────────────────────────────────────
# Alignment → Segment conversion
def process_alignment(raw_message, wanted_nodes, chrom_filter):
    segment_dict = {}
    aln = vg_pb2.Alignment()
    aln.ParseFromString(raw_message)

    # MAPQ filter
    if aln.mapping_quality <= 10:
        return segment_dict

    # Chrom filter (if provided)
    if chrom_filter and not any(pos.name == chrom_filter for pos in aln.refpos):
        return segment_dict

    read_seq = aln.sequence
    read_bq  = aln.quality
    mapq     = aln.mapping_quality
    read_off = 0

    for mapping in aln.path.mapping:
        node_id = mapping.position.node_id

        # Advance read offset even if we skip writing (keeps position in read consistent)
        if node_id not in wanted_nodes:
            for ed in mapping.edit:
                read_off += ed.to_length
            continue

        node_offset = mapping.position.offset
        strand_char = b"-" if mapping.position.is_reverse else b"+"

        seq_parts = []
        bq_parts  = bytearray()
        cig_parts = []

        for ed in mapping.edit:
            from_len = ed.from_length
            to_len   = ed.to_length
            seq_len  = len(ed.sequence)

            # CIGAR
            if from_len == to_len:
                cig_parts.append(f"{from_len}{'M' if seq_len == 0 else 'X'}")
            elif from_len > 0 and to_len == 0:
                cig_parts.append(f"{from_len}D")
            elif from_len == 0 and to_len > 0:
                cig_parts.append(f"{to_len}I")
            else:
                raise ValueError(f"Unexpected edit: from_length={from_len}, sequence_length={seq_len}")

            # Sequence/quality slices (by to_len)
            seg_len = to_len
            if seg_len > 0:
                seq_fragment = read_seq[read_off: read_off + seg_len]
                bq_fragment  = read_bq[read_off: read_off + seg_len]
                seq_parts.append(seq_fragment.upper())
                bq_parts.extend(bq_fragment)
            read_off += seg_len

        seg = Segment(
            offset=node_offset,
            seq="".join(seq_parts).encode(),
            bq=bytes(bq_parts),
            cigar="".join(cig_parts).encode(),
            rq=mapq,
            strand=strand_char
        )
        segment_dict.setdefault(node_id, []).append(seg)

    return segment_dict


# ─────────────────────────────────────────────────────────────────────────────
# Initialize outputs (variable record sizes; header-offset fixed)
def initialize_output_files(stats_path, output_prefix, default_node_length=150):
    with open(stats_path, "rb") as stats_file:
        stats_data = pickle.load(stats_file)

    block_infos = {}
    wanted_nodes = set()
    total_nodes = 0
    warned_default = False

    # First block begins after global header
    current_offset = GLOBAL_HEADER_SIZE

    for node_id_key, stat in stats_data.items():
        total_nodes += 1
        node_id = int(node_id_key)

        perfect      = int(stat.get("perfect", 0))
        not_perfect  = int(stat.get("not_perfect", 0))
        node_len     = int(stat.get("length", 0))  # read node_length from PKL

        if node_len <= 0:
            node_len = default_node_length
            if not warned_default:
                print(f"[warn] Missing node length in PKL for some nodes; using default={default_node_length}. "
                      f"Example node_id={node_id}")
                warned_default = True

        # selection rule
        total_rec = perfect + not_perfect
        if total_rec > 0 and not_perfect > 1 and (not_perfect / total_rec) > 0.10:
            wanted_nodes.add(node_id)
            n_records = total_rec
            rec_sz    = record_size(node_len)
            blk_sz    = BLOCK_HDR_SIZE + n_records * rec_sz

            block_infos[node_id] = {
                "offset": current_offset,   # absolute offset from file start
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
    with open(dat_path, "wb") as f:
        # Global header
        f.write(GLOBAL_MAGIC)
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR, GLOBAL_MINOR, len(block_infos), b'\x00' * 16))

        # Blocks with preallocated blank records
        for node_id, info in block_infos.items():
            L = info["node_length"]
            N = info["n_records"]
            f.write(BLOCK_HDR_PACK.pack(node_id, N, 0, L))

            rec_pack = make_record_struct(L)
            blank = rec_pack.pack(0, b'\x00'*L, b'\x00'*L, b'\x00'*30, 0, b'+')
            # Write N blank records (avoid giant concatenation)
            for _ in range(N):
                f.write(blank)

    # Build .idx (include node_length for convenience)
    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx_file:
        idx_file.write(struct.pack("<I", len(block_infos)))
        for node_id, info in block_infos.items():
            # node_id(u32), offset(u64), block_size(u32), n_records(u32), flags(u16), node_length(u32)
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


# NEW: load existing .dat/.idx instead of initializing
def load_existing_output_files(output_prefix):
    idx_path = output_prefix + ".idx"
    dat_path = output_prefix + ".dat"

    if not (os.path.exists(idx_path) and os.path.exists(dat_path)):
        raise FileNotFoundError(f"Expected existing files: {idx_path} and {dat_path}")

    # Read .idx header
    with open(idx_path, "rb") as f:
        raw = f.read(4)
        if len(raw) != 4:
            raise RuntimeError("Corrupt .idx: cannot read block count")
        (count,) = struct.unpack("<I", raw)
        # Determine per-entry size to support older/newer variants (22 vs 26 bytes)
        f.seek(0, os.SEEK_END)
        total = f.tell()
        remaining = total - 4
        if count <= 0:
            raise RuntimeError("Empty .idx (no blocks)")
        entry_size = remaining // count
        if entry_size not in (22, 26):
            # Default to 26 and hope; otherwise fail
            entry_size = 26
        f.seek(4)

        entries = []
        for _ in range(count):
            if entry_size == 26:
                data = f.read(26)
                if len(data) != 26: raise RuntimeError("Corrupt .idx: truncated entry")
                node_id, offset, block_size, n_records, flags, node_len = struct.unpack("<I Q I I H I", data)
            else:
                data = f.read(22)
                if len(data) != 22: raise RuntimeError("Corrupt .idx: truncated entry")
                node_id, offset, block_size, n_records, flags = struct.unpack("<I Q I I H", data)
                node_len = 0  # unknown; read from .dat
            entries.append((node_id, offset, block_size, n_records, flags, node_len))

    # Verify .dat global header and block count
    with open(dat_path, "rb") as df:
        magic = df.read(len(GLOBAL_MAGIC))
        if magic != GLOBAL_MAGIC:
            raise RuntimeError("Invalid .dat magic/version")
        majors, minors, dat_count, _ = GLOBAL_VER_PACK.unpack(df.read(GLOBAL_VER_PACK.size))
        if dat_count != len(entries):
            print(f"[warn] .dat block_count ({dat_count}) != .idx count ({len(entries)})")

        # Build block_infos, pulling node_length from .dat when missing
        block_infos = {}
        wanted_nodes = set()
        for node_id, offset, block_size, n_records, flags, node_len in entries:
            if node_len <= 0:
                df.seek(offset, os.SEEK_SET)
                hdr = df.read(BLOCK_HDR_SIZE)
                if len(hdr) != BLOCK_HDR_SIZE:
                    raise RuntimeError(f"Corrupt .dat: cannot read block header at {offset}")
                nid2, nrec2, flg2, node_len = BLOCK_HDR_PACK.unpack(hdr)
                # light validation
                if nid2 != node_id or nrec2 != n_records:
                    print(f"[warn] .dat/.idx mismatch for node {node_id} (idx n={n_records}, dat n={nrec2})")

            rec_sz = record_size(node_len)
            block_infos[node_id] = {
                "offset": offset,
                "n_records": n_records,
                "current_pos": 0,          # assume unused/empty file
                "node_length": node_len,
                "record_size": rec_sz,
                "block_size": block_size,
            }
            wanted_nodes.add(node_id)

    print(f"Reusing existing output with {len(block_infos)} node blocks from {idx_path}")
    return block_infos, dat_path, wanted_nodes


# ─────────────────────────────────────────────────────────────────────────────
# Faster flusher: pack_into + pwrite + file-order + optional concurrency
ZEROS = b"\x00" * (1024 * 1024)  # 1 MB shared pad buffer

def _ensure_zeros(n):
    global ZEROS
    if len(ZEROS) < n:
        new_len = max(n, len(ZEROS) * 2)
        ZEROS = b"\x00" * new_len

def make_flusher(block_infos, dat_fh, io_workers=8, fsync_every=0, verify_bounds=True):
    """
    Returns a flush() closure that writes all buffered segments efficiently.
    - io_workers: number of concurrent pwrite workers (1 = synchronous).
    - fsync_every: call os.fsync(fd) every N flushes (0 = never).
    - verify_bounds: extra assertions to catch offset math mistakes.
    """
    fd = dat_fh.fileno()
    supports_pwrite = hasattr(os, "pwrite")
    flush_counter = {"n": 0}  # mutable box for closure

    # cache for struct objects per node length
    pack_cache = {}
    def get_pack(L):
        st = pack_cache.get(L)
        if st is None:
            st = make_record_struct(L)
            pack_cache[L] = st
        return st

    def flush(segment_buffer):
        if not segment_buffer:
            return

        # Order nodes by on-disk offset for locality
        items = sorted(segment_buffer.items(), key=lambda kv: block_infos[kv[0]]["offset"])

        write_tasks = []  # (base_pos, buf, node_id, batch_n)
        for node_id, segs in items:
            if not segs:
                continue
            info = block_infos[node_id]
            L        = info["node_length"]
            rec_size = info["record_size"]
            rec_pack = get_pack(L)

            base_pos   = info["offset"] + BLOCK_HDR_SIZE + info["current_pos"] * rec_size
            batch_n    = len(segs)
            total_size = batch_n * rec_size

            # Bound checks (avoid silent corruption)
            if verify_bounds:
                assert info["current_pos"] + batch_n <= info["n_records"], \
                    f"node {node_id}: write past reserved records ({info['current_pos']} + {batch_n} > {info['n_records']})"
                block_end = info["offset"] + info["block_size"]
                assert base_pos + total_size <= block_end, \
                    f"node {node_id}: write exceeds block (end={base_pos+total_size} > {block_end})"

            # Pre-allocate batch buffer and pack in-place
            buf = bytearray(total_size)
            _ensure_zeros(max(L, 30))

            def pad_exact(b: bytes, n: int) -> bytes:
                lb = len(b)
                if lb >= n:
                    return b[:n]
                return b + ZEROS[:n - lb]

            off = 0
            for seg in segs:
                seq_f = pad_exact(seg.seq, L)
                bq_f  = pad_exact(seg.bq,  L)
                cg_f  = pad_exact(seg.cigar, 30)

                rec_pack.pack_into(
                    buf, off,
                    int(seg.offset),
                    seq_f, bq_f, cg_f,
                    int(seg.rq),
                    seg.strand if seg.strand in (b'+', b'-') else b'+'
                )
                off += rec_size

            write_tasks.append((base_pos, buf, node_id, batch_n))

        # Perform writes
        if supports_pwrite and (io_workers or 0) > 1 and len(write_tasks) > 1:
            workers = min(io_workers, len(write_tasks))
            with ThreadPoolExecutor(max_workers=workers) as ex:
                futs = [ex.submit(os.pwrite, fd, buf, base_pos) for base_pos, buf, _, _ in write_tasks]
                for fut in as_completed(futs):
                    _ = fut.result()
        else:
            # Fallback: single-threaded writes (seek+write if pwrite not available)
            for base_pos, buf, _, _ in write_tasks:
                if supports_pwrite:
                    os.pwrite(fd, buf, base_pos)
                else:
                    dat_fh.seek(base_pos, os.SEEK_SET)
                    dat_fh.write(buf)

        # Update write positions after success
        for _, _, node_id, batch_n in write_tasks:
            block_infos[node_id]["current_pos"] += batch_n

        # Optional durability checkpoint
        flush_counter["n"] += 1
        if fsync_every and (flush_counter["n"] % fsync_every == 0):
            os.fsync(fd)

        # Clear the caller's buffer
        segment_buffer.clear()

    return flush


# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter,
                 buffer_segments, io_workers, fsync_every, verify_bounds, use_existing):
    if use_existing:
        print("Reusing existing .dat/.idx...")
        block_infos, dat_path, wanted_nodes = load_existing_output_files(output_prefix)
    else:
        print("Initializing output files...")
        block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)
        print(f"Output file created: {dat_path}")

    next_milestone = milestone_step
    total_reads = 0
    total_segments = 0
    start_time = time.perf_counter()

    dat_fh = open(dat_path, "r+b")
    segment_buffer = defaultdict(list)
    flush = make_flusher(
        block_infos, dat_fh,
        io_workers=io_workers,
        fsync_every=fsync_every,
        verify_bounds=verify_bounds
    )

    for raw_msg in gam_record_iter(gam_path):
        segment_dict = process_alignment(raw_msg, wanted_nodes, chrom_filter)
        total_reads += 1

        for node_id, segs in segment_dict.items():
            segment_buffer[node_id].extend(segs)
            total_segments += len(segs)

        # Flush when buffered segment count exceeds threshold
        if total_segments >= buffer_segments:
            flush(segment_buffer)
            segment_buffer = defaultdict(list)
            total_segments = 0

        if total_reads >= next_milestone:
            elapsed = time.perf_counter() - start_time
            print(f"{total_reads} reads processed | {elapsed:.1f} seconds")
            next_milestone += milestone_step

    # Final flush
    if segment_buffer:
        flush(segment_buffer)

    dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")


# ─────────────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description="GAM segment extractor (variable per-node record size, node_length from PKL, fast flush)."
    )
    p.add_argument("gam_path", help="Path to the GAM file")
    p.add_argument("stats_pickle", help="PKL with per-node {'perfect','not_perfect','length'} (length preferred)")
    p.add_argument("output_prefix", help="Prefix for output files (.dat/.idx)")
    p.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval in reads")
    p.add_argument("--chr", default="", help="Optional chromosome filter (matches Alignment.refpos.name)")
    p.add_argument("--buffer-segments", type=int, default=400_000_000,
                   help="Flush when buffered segments reach this count (default: 100,000,000)")
    p.add_argument("--io-workers", type=int, default=8,
                   help="Concurrent pwrite workers (1 disables threading)")
    p.add_argument("--fsync-every", type=int, default=0,
                   help="Call fsync() after every N flushes (0 disables)")
    p.add_argument("--no-verify-bounds", action="store_true",
                   help="Disable extra assertions on write bounds (faster, less safe)")
    p.add_argument("--use-existing", action="store_true",
                   help="Reuse existing initialized output (output_prefix.dat/.idx) instead of reinitializing")
    args = p.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr,
        buffer_segments=args.buffer_segments,
        io_workers=max(1, args.io_workers),
        fsync_every=max(0, args.fsync_every),
        verify_bounds=not args.no_verify_bounds,
        use_existing=args.use_existing
    )


if __name__ == "__main__":
    main()

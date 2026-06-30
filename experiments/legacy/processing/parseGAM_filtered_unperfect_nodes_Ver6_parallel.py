#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
import os
from collections import defaultdict
import vg_pb2  # Assumes you have the compiled protobuf Python file
import multiprocessing as mp

# ─────────────────────────────────────────────────────────────────────────────
# Segment container
class Segment:
    __slots__ = ('offset', 'seq', 'bq', 'cigar', 'rq', 'strand')
    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset = offset
        self.seq    = seq      # bytes (unpadded)
        self.bq     = bq       # bytes (unpadded)
        self.cigar  = cigar    # bytes (unpadded)
        self.rq     = rq       # int (MAPQ)
        self.strand = strand   # b'+' or b'-'

# ─────────────────────────────────────────────────────────────────────────────
# File layout

# Global header: magic + version block
GLOBAL_MAGIC = b"MYFMT\x01"                                   # 6 bytes: "MYFMT" + 0x01
GLOBAL_VER_PACK = struct.Struct("<BBI16s")                    # major, minor, block_count, reserved[16]
GLOBAL_MAJOR, GLOBAL_MINOR = 0, 5                             # minor=5: no node_length; per-block maxima used
GLOBAL_HEADER_SIZE = len(GLOBAL_MAGIC) + GLOBAL_VER_PACK.size # 6 + 22 = 28

# Per-node block header (NO node_length)
# <I I H I I> -> node_id (u32), n_records (u32), flags (u16),
#                max_read_length (u32), max_cigar_length (u32)
BLOCK_HDR_PACK = struct.Struct("<I I H I I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size  # 18

def make_record_struct(max_read_len: int, max_cigar_len: int) -> struct.Struct:
    """
    Per-record layout:
      <h {R}s {R}s {C}s h c
        - i16 offset
        - seq[R] bytes
        - bq[R] bytes
        - cigar[C] bytes (ASCII)
        - i16 rq (MAPQ)
        - char strand ('+' / '-')
    """
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")

def record_size(max_read_len: int, max_cigar_len: int) -> int:
    return make_record_struct(max_read_len, max_cigar_len).size

# ─────────────────────────────────────────────────────────────────────────────
# GAM parsing helpers

def read_varint(stream):
    value = 0
    shift_amount = 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError("EOF while reading varint")
        v = b[0]
        value |= (v & 0x7F) << shift_amount
        if not (v & 0x80):
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
            if group_count == 0: continue
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
# Alignment → Segment
def process_alignment(raw_message, wanted_nodes, chrom_filter):
    segment_dict = {}
    aln = vg_pb2.Alignment()
    aln.ParseFromString(raw_message)

    if aln.mapping_quality <= 10:
        return segment_dict
    if chrom_filter and not any(pos.name == chrom_filter for pos in aln.refpos):
        return segment_dict

    read_sequence = aln.sequence
    read_quality  = aln.quality
    mapq = aln.mapping_quality
    read_offset = 0

    for mapping in aln.path.mapping:
        nid = mapping.position.node_id

        if nid not in wanted_nodes:
            for e in mapping.edit:
                read_offset += e.to_length
            continue

        node_offset = mapping.position.offset
        strand_char = b"-" if mapping.position.is_reverse else b"+"
        seq_parts, bq_parts, cigar_parts = [], bytearray(), []

        for e in mapping.edit:
            fL, tL = e.from_length, e.to_length
            if fL == tL: cigar_parts.append(f"{fL}M" if not e.sequence else f"{fL}X")
            elif fL > 0 and tL == 0: cigar_parts.append(f"{fL}D")
            elif fL == 0 and tL > 0: cigar_parts.append(f"{tL}I")
            else: raise ValueError(f"Unexpected edit: from_length={fL}, to_length={tL}")

            if tL > 0:
                seq_parts.append(read_sequence[read_offset : read_offset + tL].upper())
                bq_parts.extend(read_quality[read_offset : read_offset + tL])
                read_offset += tL

        seg = Segment(offset=node_offset, seq="".join(seq_parts).encode(), bq=bytes(bq_parts),
                      cigar="".join(cigar_parts).encode(), rq=mapq, strand=strand_char)
        segment_dict.setdefault(nid, []).append(seg)

    return segment_dict

# ─────────────────────────────────────────────────────────────────────────────
# File Initialization (Serial)
def initialize_output_files(stats_path, output_prefix):
    with open(stats_path, "rb") as fh:
        stats_data = pickle.load(fh)

    wanted_nodes, node_counts, maxima = set(), {}, {}
    for node_id_key, stat in stats_data.items():
        nid = int(node_id_key)
        perfect = int(stat.get("perfect", 0))
        not_perfect = int(stat.get("not_perfect", 0))
        if (perfect + not_perfect) > 0 and not_perfect > 1 and not_perfect / (perfect + not_perfect) > 0.05:
            wanted_nodes.add(nid)
            node_counts[nid] = perfect + not_perfect
            R = int(stat.get("max_read_length", 1) or 1)
            C = int(stat.get("max_cigar_length", 1) or 1)
            maxima[nid] = (max(1, R), max(1, C))

    print(f"Filtered {len(wanted_nodes)} nodes from {len(stats_data)} total.")
    del stats_data
    gc.collect()

    block_infos = {}
    current_offset = GLOBAL_HEADER_SIZE
    for nid in sorted(list(wanted_nodes)): # Sort for deterministic layout
        nrec = node_counts[nid]
        R, C = maxima[nid]
        rec_sz = record_size(R, C)
        blk_sz = BLOCK_HDR_SIZE + nrec * rec_sz
        block_infos[nid] = {"offset": current_offset, "n_records": nrec, "current_pos": 0,
                            "max_read_len": R, "max_cigar_len": C, "record_size": rec_sz, "block_size": blk_sz}
        current_offset += blk_sz

    dat_path = output_prefix + ".dat"
    with open(dat_path, "wb") as f:
        f.write(GLOBAL_MAGIC)
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR, GLOBAL_MINOR, len(block_infos), b'\x00' * 16))
        for nid, info in block_infos.items():
            f.write(BLOCK_HDR_PACK.pack(nid, info["n_records"], 0, info["max_read_len"], info["max_cigar_len"]))
        # Pre-allocate the file by seeking to the end and writing a null byte
        if current_offset > GLOBAL_HEADER_SIZE:
            f.seek(current_offset - 1)
            f.write(b'\x00')

    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx:
        idx.write(struct.pack("<I", len(block_infos)))
        for nid, info in block_infos.items():
            idx.write(struct.pack("<I Q I I H I I", nid, info["offset"], info["block_size"],
                                  info["n_records"], 0, info["max_read_len"], info["max_cigar_len"]))

    return block_infos, dat_path, wanted_nodes

def load_existing_output_files(output_prefix):
    # This function remains unchanged from your script
    # ... (implementation from the question)
    pass # Placeholder for brevity

# ─────────────────────────────────────────────────────────────────────────────
# WORKER GLOBALS & FUNCTIONS (Must be at the top level for pickling)

# Globals to be populated in each worker process by the initializer
g_fd = -1
g_block_infos = None
g_locks = None
g_wanted_nodes = None
g_chrom_filter = None

def init_worker(fd, block_infos, locks, wanted_nodes, chrom_filter):
    """Initializes globals for each worker process in the pool."""
    global g_fd, g_block_infos, g_locks, g_wanted_nodes, g_chrom_filter
    g_fd = fd
    g_block_infos = block_infos
    g_locks = locks
    g_wanted_nodes = wanted_nodes
    g_chrom_filter = chrom_filter

def process_and_write_worker(raw_msg):
    """
    A single worker function that processes a GAM record and writes the output.
    This is the core of the streaming pipeline.
    """
    segs_by_node = process_alignment(raw_msg, g_wanted_nodes, g_chrom_filter)
    if not segs_by_node:
        return 0

    segments_written = 0
    for nid, segs in segs_by_node.items():
        if not segs: continue

        # CRITICAL SECTION: Atomically claim a write position for this node.
        with g_locks[nid]:
            info = g_block_infos[nid]
            current_write_idx = info["current_pos"]
            n_segs_to_write = len(segs)
            info["current_pos"] = current_write_idx + n_segs_to_write
            g_block_infos[nid] = info # Write back the modified dictionary

        # The lock is now released. We can proceed with serialization and I/O.
        R, C = info["max_read_len"], info["max_cigar_len"]
        rec_pack = make_record_struct(R, C)
        rec_sz = info["record_size"]
        buf = bytearray(rec_sz * n_segs_to_write)
        off = 0
        for s in segs:
            rec_pack.pack_into(buf, off, int(s.offset), s.seq.ljust(R, b'\x00')[:R], s.bq.ljust(R, b'\x00')[:R], s.cigar.ljust(C, b'\x00')[:C], int(s.rq), s.strand if s.strand in (b'+', b'-') else b'+')
            off += rec_sz

        base_offset = info["offset"] + BLOCK_HDR_SIZE
        write_pos = base_offset + current_write_idx * rec_sz
        os.pwrite(g_fd, buf, write_pos)
        segments_written += n_segs_to_write

    return segments_written

# ─────────────────────────────────────────────────────────────────────────────
# Main Pipeline
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing, num_workers):
    if not hasattr(os, 'pwrite'):
        raise NotImplementedError("This parallel script requires os.pwrite(), which is not available on this OS (e.g., Windows).")

    if use_existing:
        print("Reusing existing .dat/.idx...")
        block_infos, dat_path, wanted_nodes = load_existing_output_files(output_prefix)
    else:
        print("Initializing output files from PKL maxima...")
        block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)

    start_time = time.perf_counter()
    dat_fh = open(dat_path, "r+b")
    dat_fd = dat_fh.fileno()

    # ---- SETUP SHARED STATE USING A MANAGER ----
    print("Setting up shared memory manager for locks and state...")
    manager = mp.Manager()
    mgr_block_infos = manager.dict(block_infos)
    mgr_locks = manager.dict({nid: manager.Lock() for nid in wanted_nodes})
    print("Manager setup complete.")

    total_reads = 0
    next_milestone = milestone_step
    init_args = (dat_fd, mgr_block_infos, mgr_locks, wanted_nodes, chrom_filter)

    with mp.Pool(processes=num_workers, initializer=init_worker, initargs=init_args) as pool:
        gam_iterator = gam_record_iter(gam_path)
        # Use imap_unordered to process records as they come in, improving throughput.
        # chunksize helps reduce the overhead of sending tasks to workers.
        for _ in pool.imap_unordered(process_and_write_worker, gam_iterator, chunksize=250):
            total_reads += 1
            if total_reads >= next_milestone:
                elapsed = time.perf_counter() - start_time
                print(f"{total_reads:,} reads processed | {elapsed:.1f} seconds")
                next_milestone += milestone_step

    dat_fh.close()
    manager.shutdown()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads:,}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")

# ─────────────────────────────────────────────────────────────────────────────
# Entry Point
def main():
    parser = argparse.ArgumentParser(
        description="Streaming parallel GAM segment extractor with per-node locking.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the stats pickle with node maxima")
    parser.add_argument("output_prefix", help="Prefix for output files (.dat, .idx)")
    parser.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval (number of reads)")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    parser.add_argument("--use-existing", action="store_true", help="Reuse existing initialized output files")
    parser.add_argument("--workers", type=int, default=min(4, os.cpu_count()), help="Number of worker processes to use")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr,
        use_existing=args.use_existing,
        num_workers=args.workers
    )

if __name__ == "__main__":
    main()
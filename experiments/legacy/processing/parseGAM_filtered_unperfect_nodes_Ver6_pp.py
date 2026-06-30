#!/usr/bin/env python3
# Option B (mmap version): parent packs buffers; workers write to mmap at fixed offsets.

import argparse
import gzip
import pickle
import struct
import time
import gc
import os
import mmap
from collections import defaultdict
import multiprocessing as mp
import vg_pb2  # compiled protobuf for VG Alignment

# ─────────────────────────────────────────────────────────────────────────────
# Segment container
class Segment:
    __slots__ = ("offset", "seq", "bq", "cigar", "rq", "strand")

    def __init__(self, offset, seq, bq, cigar, rq, strand):
        self.offset = offset     # int
        self.seq = seq           # bytes (unpadded)
        self.bq = bq             # bytes (unpadded)
        self.cigar = cigar       # bytes (unpadded)
        self.rq = rq             # int (MAPQ)
        self.strand = strand     # b'+' or b'-'


# ─────────────────────────────────────────────────────────────────────────────
# File layout
GLOBAL_MAGIC = b"MYFMT\x01"
GLOBAL_VER_PACK = struct.Struct("<BBI16s")  # major, minor, nblocks, reserved16
GLOBAL_MAJOR, GLOBAL_MINOR = 0, 5
GLOBAL_HEADER_SIZE = len(GLOBAL_MAGIC) + GLOBAL_VER_PACK.size
BLOCK_HDR_PACK = struct.Struct("<I I H I I")  # node_id, n_records, flags, maxR, maxC
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size


def make_record_struct(max_read_len: int, max_cigar_len: int) -> struct.Struct:
    # <h seq[R] bq[R] cigar[C] h c>
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")


def record_size(max_read_len: int, max_cigar_len: int) -> int:
    return make_record_struct(max_read_len, max_cigar_len).size


# ─────────────────────────────────────────────────────────────────────────────
# GAM parsing
def read_varint(stream):
    value, shift_amount = 0, 0
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
            if group_count == 0:
                continue
            try:
                tag_len = read_varint(f)
                group_tag = f.read(tag_len).decode()
            except (EOFError, UnicodeDecodeError):
                break
            if group_tag != tag:
                for _ in range(group_count - 1):
                    try:
                        f.seek(read_varint(f), 1)
                    except EOFError:
                        break
                continue
            for _ in range(group_count - 1):
                try:
                    size = read_varint(f)
                    yield f.read(size)
                except EOFError:
                    break


def process_alignment(raw_message, wanted_nodes, chrom_filter):
    """
    Returns: dict[node_id] -> [Segment, ...] for segments belonging to wanted_nodes
    """
    out = {}
    aln = vg_pb2.Alignment()
    aln.ParseFromString(raw_message)

    if aln.mapping_quality <= 10:
        return out
    if chrom_filter and not any(getattr(rp, "name", "") == chrom_filter for rp in aln.refpos):
        return out

    seq = aln.sequence
    bqs = aln.quality
    mapq = aln.mapping_quality
    read_off = 0

    for mapping in aln.path.mapping:
        nid = mapping.position.node_id
        if nid not in wanted_nodes:
            for e in mapping.edit:
                read_off += e.to_length
            continue

        node_off = mapping.position.offset
        strand = b"-" if mapping.position.is_reverse else b"+"
        seq_parts = []
        bq_parts = bytearray()
        cigar_parts = []

        for e in mapping.edit:
            fL, tL = e.from_length, e.to_length
            if fL == tL:
                cigar_parts.append(f"{fL}M" if not e.sequence else f"{fL}X")
            elif fL > 0 and tL == 0:
                cigar_parts.append(f"{fL}D")
            elif fL == 0 and tL > 0:
                cigar_parts.append(f"{tL}I")

            if tL > 0:
                seq_parts.append(seq[read_off : read_off + tL].upper())
                bq_parts.extend(bqs[read_off : read_off + tL])
                read_off += tL

        seg = Segment(
            offset=node_off,
            seq=("".join(seq_parts)).encode(),
            bq=bytes(bq_parts),
            cigar=("".join(cigar_parts)).encode(),
            rq=mapq,
            strand=strand,
        )
        out.setdefault(nid, []).append(seg)

    return out


# ─────────────────────────────────────────────────────────────────────────────
# Output initialization
def initialize_output_files(stats_path, output_prefix):
    with open(stats_path, "rb") as fh:
        stats = pickle.load(fh)

    wanted_nodes = set()
    node_counts = {}
    maxima = {}

    for node_id_key, stat in stats.items():
        nid = int(node_id_key)
        perfect = int(stat.get("perfect", 0))
        not_perfect = int(stat.get("not_perfect", 0))
        total = perfect + not_perfect
        if total > 0 and not_perfect > 1 and not_perfect / total > 0.05:
            wanted_nodes.add(nid)
            node_counts[nid] = total
            R = int(stat.get("max_read_length", 1) or 1)
            C = int(stat.get("max_cigar_length", 1) or 1)
            maxima[nid] = (max(1, R), max(1, C))

    print(f"Filtered {len(wanted_nodes)} nodes from {len(stats)} total.")
    del stats
    gc.collect()

    block_infos = {}
    current_offset = GLOBAL_HEADER_SIZE
    for nid in sorted(wanted_nodes):  # deterministic block order
        nrec = node_counts[nid]
        R, C = maxima[nid]
        rec_sz = record_size(R, C)
        blk_sz = BLOCK_HDR_SIZE + nrec * rec_sz
        block_infos[nid] = {
            "offset": current_offset,
            "n_records": nrec,
            "current_pos": 0,
            "max_read_len": R,
            "max_cigar_len": C,
            "record_size": rec_sz,
            "block_size": blk_sz,
        }
        current_offset += blk_sz

    # Create and preallocate .dat
    dat_path = output_prefix + ".dat"
    with open(dat_path, "wb") as f:
        # 1) write global header
        f.write(GLOBAL_MAGIC)
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR, GLOBAL_MINOR, len(block_infos), b'\x00' * 16))

        # 2) pre-allocate whole file
        if current_offset > GLOBAL_HEADER_SIZE:
            f.seek(current_offset - 1)
            f.write(b'\x00')

        # 3) write each block header at its declared offset
        for nid, info in block_infos.items():
            f.seek(info["offset"], os.SEEK_SET)
            f.write(BLOCK_HDR_PACK.pack(
                nid,
                info["n_records"],
                0,  # flags
                info["max_read_len"],
                info["max_cigar_len"],
            ))

    # Write .idx
    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx:
        idx.write(struct.pack("<I", len(block_infos)))
        for nid, info in block_infos.items():
            idx.write(
                struct.pack(
                    "<I Q I I H I I",
                    nid,
                    info["offset"],
                    info["block_size"],
                    info["n_records"],
                    0,
                    info["max_read_len"],
                    info["max_cigar_len"],
                )
            )

    total_size = current_offset  # final byte length of the .dat file
    return block_infos, dat_path, wanted_nodes, total_size


def load_existing_output_files(output_prefix):
    raise NotImplementedError("--use-existing loader not implemented in this example.")


# ─────────────────────────────────────────────────────────────────────────────
# Worker: mmap writer
g_map = None  # type: mmap.mmap
g_map_len = 0
g_dat_path = None


def init_flush_worker_map(dat_path: str, map_len: int):
    """Each worker opens and mmaps the file for write access."""
    global g_map, g_map_len, g_dat_path
    g_dat_path = dat_path
    g_map_len = map_len
    # r+b required; size must be already preallocated
    fh = open(dat_path, "r+b", buffering=0)
    # ACCESS_WRITE maps file for read/write and keeps changes visible to other processes
    g_map = mmap.mmap(fh.fileno(), map_len, access=mmap.ACCESS_WRITE)
    # We intentionally keep fh open (mapped file keeps FD ref); OS will close on process exit.


def flush_worker_map(job):
    """job: (write_pos:int, buf:bytes) → write into mmap slice."""
    write_pos, buf = job
    end_pos = write_pos + len(buf)
    if end_pos > g_map_len:
        raise RuntimeError(f"mmap write out of bounds: [{write_pos}, {end_pos}) > {g_map_len}")
    g_map[write_pos:end_pos] = buf
    # Rely on OS page cache; explicit flush not needed per write. Parent will flush at end.


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing, num_flush_workers):
    if use_existing:
        block_infos, dat_path, wanted_nodes, total_size = load_existing_output_files(output_prefix)
    else:
        block_infos, dat_path, wanted_nodes, total_size = initialize_output_files(stats_path, output_prefix)

    BUFFER_SEGMENTS = 10_000_000  # total segments batched before a flush
    next_milestone = milestone_step
    total_reads = 0
    total_segments = 0
    start_time = time.perf_counter()

    # Parent opens file too (so we can flush at end)
    dat_fh = open(dat_path, "r+b", buffering=0)

    # Prefer "fork" on Linux; fall back if unavailable (spawn will still work because
    # each worker mmaps the path independently).
    try:
        mpc = mp.get_context("fork")
    except ValueError:
        mpc = mp.get_context()

    pool = mpc.Pool(
        processes=num_flush_workers,
        initializer=init_flush_worker_map,
        initargs=(dat_path, total_size),
    )

    segment_buffer = defaultdict(list)

    def flush_segment_buffer_parallel():
        nonlocal total_segments
        if not segment_buffer:
            return

        n_nodes = len(segment_buffer)
        n_segs = sum(len(v) for v in segment_buffer.values())
        print(f"[flush] start: nodes={n_nodes:,}, segs={n_segs:,}")
        t0 = time.perf_counter()

        jobs = []
        planned_bytes = 0

        # Pop items until buffer is empty
        while segment_buffer:
            nid, segs = segment_buffer.popitem()
            if not segs:
                continue
            info = block_infos.get(nid)
            if not info:
                continue

            if info["current_pos"] + len(segs) > info["n_records"]:
                raise RuntimeError(
                    f"Block overflow for node {nid}: "
                    f"{info['current_pos']} + {len(segs)} > {info['n_records']}"
                )

            base_offset = info["offset"] + BLOCK_HDR_SIZE
            R, C = info["max_read_len"], info["max_cigar_len"]
            rec_pack = make_record_struct(R, C)
            rec_sz = rec_pack.size

            buf = bytearray(rec_sz * len(segs))
            off = 0
            for s in segs:
                rec_pack.pack_into(
                    buf, off,
                    int(s.offset),
                    s.seq.ljust(R, b"\x00")[:R],
                    s.bq.ljust(R, b"\x00")[:R],
                    s.cigar.ljust(C, b"\x00")[:C],
                    int(s.rq),
                    s.strand if s.strand in (b"+", b"-") else b"+",
                )
                off += rec_sz

            write_pos = base_offset + info["current_pos"] * rec_sz
            info["current_pos"] += len(segs)

            planned_bytes += len(buf)
            jobs.append((write_pos, bytes(buf)))

        t_pack = time.perf_counter()
        if jobs:
            print(f"[flush] dispatching {len(jobs):,} jobs, total≈{planned_bytes / 1e6:.2f} MB")
            pool.map(flush_worker_map, jobs, chunksize=32)

        t_end = time.perf_counter()
        print(
            f"[flush] done: segs={n_segs:,}, "
            f"pack={t_pack - t0:.3f}s, io={t_end - t_pack:.3f}s, total={t_end - t0:.3f}s"
        )

        total_segments = 0

    try:
        for raw_msg in gam_record_iter(gam_path):
            segs_by_node = process_alignment(raw_msg, wanted_nodes, chrom_filter)
            total_reads += 1

            for nid, segs in segs_by_node.items():
                segment_buffer[nid].extend(segs)
                total_segments += len(segs)

            if total_segments >= BUFFER_SEGMENTS:
                flush_segment_buffer_parallel()

            if total_reads >= next_milestone:
                elapsed = time.perf_counter() - start_time
                print(f"{total_reads:,} reads processed | {elapsed:.1f} s")
                next_milestone += milestone_step

        # final drain
        flush_segment_buffer_parallel()
    finally:
        pool.close()
        pool.join()

        # Parent-side final flush to ensure all dirty pages are persisted
        try:
            # Map briefly and flush; or fsync the file descriptor.
            with mmap.mmap(dat_fh.fileno(), length=0, access=mmap.ACCESS_READ) as m:
                m.flush()
        except Exception:
            try:
                os.fsync(dat_fh.fileno())
            except Exception:
                pass

        dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads:,}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
def main():
    p = argparse.ArgumentParser(
        description="GAM segment extractor with parallel mmap writes (parent packs).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("gam_path", help="Path to the GAM file")
    p.add_argument("stats_pickle", help="Path to the stats pickle with node maxima")
    p.add_argument("output_prefix", help="Prefix for output files (.dat, .idx)")
    p.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval")
    p.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    p.add_argument("--use-existing", action="store_true", help="Reuse existing initialized output files")
    p.add_argument("--flush-workers", type=int, default=4, help="Workers for I/O (mmap slice writes)")
    args = p.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr,
        use_existing=args.use_existing,
        num_flush_workers=args.flush_workers,
    )


if __name__ == "__main__":
    main()

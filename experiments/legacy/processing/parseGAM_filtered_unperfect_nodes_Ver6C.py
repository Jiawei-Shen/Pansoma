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
    if max_read_len <= 0 or max_cigar_len <= 0:
        raise ValueError("max_read_len and max_cigar_len must be > 0")
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
# CIGAR builder (same semantics; not used by writer but handy for debug)
def build_cigar(mapping_edits):
    parts = []
    for e in mapping_edits:
        fL, tL, sL = e.from_length, e.to_length, len(e.sequence)
        if fL == tL:
            parts.append(f"{fL}M" if sL == 0 else f"{fL}X")
        elif fL > tL:
            parts.append(f"{fL - tL}D")
        elif fL < tL:
            parts.append(f"{tL - fL}I")
        else:
            raise ValueError(f"Unexpected edit: from_length={fL}, to_length={tL}")
    return "".join(parts)

# ─────────────────────────────────────────────────────────────────────────────
# Alignment → Segment (packing happens later using per-block maxima)
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

        # advance whether wanted or not
        if nid not in wanted_nodes:
            for e in mapping.edit:
                read_offset += e.to_length
            continue

        node_offset = mapping.position.offset
        strand_char = b"-" if mapping.position.is_reverse else b"+"

        seq_parts = []
        bq_parts  = bytearray()
        cigar_parts = []

        for e in mapping.edit:
            fL, tL, sL = e.from_length, e.to_length, len(e.sequence)

            if fL == tL:
                cigar_parts.append(f"{fL}M" if sL == 0 else f"{fL}X")
            elif fL > 0 and tL == 0:
                cigar_parts.append(f"{fL}D")
            elif fL == 0 and tL > 0:
                cigar_parts.append(f"{tL}I")
            else:
                raise ValueError(f"Unexpected edit: from_length={fL}, sequence_length={sL}")

            if tL:
                seq_frag = read_sequence[read_offset: read_offset + tL]
                bq_frag  = read_quality[read_offset: read_offset + tL]
                seq_parts.append(seq_frag.upper())
                bq_parts.extend(bq_frag)
                read_offset += tL

        seg = Segment(
            offset=node_offset,
            seq="".join(seq_parts).encode(),
            bq=bytes(bq_parts),
            cigar="".join(cigar_parts).encode(),
            rq=mapq,
            strand=strand_char
        )
        segment_dict.setdefault(nid, []).append(seg)

    return segment_dict

# ─────────────────────────────────────────────────────────────────────────────
# Initialize outputs (read per-node maxima directly from stats PKL)
def initialize_output_files(stats_path, output_prefix):
    with open(stats_path, "rb") as fh:
        stats_data = pickle.load(fh)

    # Expecting: stats_data[node_id] has keys:
    #   perfect, not_perfect, max_read_length, max_cigar_length
    wanted_nodes = set()
    node_counts = {}
    maxima = {}
    total_nodes = 0

    for node_id_key, stat in stats_data.items():
        total_nodes += 1
        nid = int(node_id_key)
        perfect = int(stat.get("perfect", 0))
        not_perfect = int(stat.get("not_perfect", 0))
        # Selection rule (keep your original threshold 0.06 here)
        if (perfect + not_perfect) > 0 and not_perfect > 1 and not_perfect / (perfect + not_perfect) > 0.05:
            wanted_nodes.add(nid)
            node_counts[nid] = perfect + not_perfect
            # Pull maxima directly from PKL (fallback to 1 if missing/bad)
            R = int(stat.get("max_read_length", 1) or 1)
            C = int(stat.get("max_cigar_length", 1) or 1)
            if R <= 0: R = 1
            if C <= 0: C = 1
            maxima[nid] = (R, C)

    print(f"Filtered {len(wanted_nodes)} nodes from {total_nodes} total nodes "
          f"({(len(wanted_nodes) / max(total_nodes,1)):.2%} selected).")
    del stats_data
    gc.collect()

    # lay out blocks using maxima from PKL
    block_infos = {}
    current_offset = GLOBAL_HEADER_SIZE
    for nid in wanted_nodes:
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

    dat_path = output_prefix + ".dat"

    # write .dat with sparse preallocation (no per-record blank writes)
    with open(dat_path, "wb") as f:
        f.write(GLOBAL_MAGIC)
        f.write(GLOBAL_VER_PACK.pack(GLOBAL_MAJOR, GLOBAL_MINOR, len(block_infos), b'\x00' * 16))

        # Keep deterministic order for writing: rely on insertion order from block_infos
        for nid, info in block_infos.items():
            nrec = info["n_records"]
            R = info["max_read_len"]
            C = info["max_cigar_len"]

            # block header
            f.write(BLOCK_HDR_PACK.pack(nid, nrec, 0, R, C))

            # sparse region for nrec records
            rec_sz = info["record_size"]
            if nrec > 0:
                # Seek from current position (right after header) to last byte of the block payload
                f.seek(nrec * rec_sz - 1, os.SEEK_CUR)
                f.write(b"\x00")
            # else: empty block payload

    # write .idx
    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx:
        idx.write(struct.pack("<I", len(block_infos)))
        # node_id (u32), offset (u64), block_size (u32), n_records (u32),
        # flags (u16), max_read_len (u32), max_cigar_len (u32)  → 30 bytes per entry
        for nid, info in block_infos.items():
            idx.write(struct.pack(
                "<I Q I I H I I",
                nid,
                info["offset"],
                info["block_size"],
                info["n_records"],
                0,
                info["max_read_len"],
                info["max_cigar_len"]
            ))

    return block_infos, dat_path, wanted_nodes

# ─────────────────────────────────────────────────────────────────────────────
# Reuse existing .dat/.idx (new format only)
def load_existing_output_files(output_prefix):
    idx_path = output_prefix + ".idx"
    dat_path = output_prefix + ".dat"
    if not (os.path.exists(idx_path) and os.path.exists(dat_path)):
        raise FileNotFoundError(f"Expected existing files: {idx_path} and {dat_path}")

    with open(idx_path, "rb") as f:
        raw = f.read(4)
        if len(raw) != 4:
            raise RuntimeError("Corrupt .idx: cannot read block count")
        (count,) = struct.unpack("<I", raw)
        entries = []
        for _ in range(count):
            data = f.read(30)  # fixed new entry size
            if len(data) != 30:
                raise RuntimeError("Corrupt .idx: truncated entry (expect 30 bytes/entry)")
            nid, offset, block_size, n_records, flags, max_read_len, max_cigar_len = struct.unpack("<I Q I I H I I", data)
            entries.append((nid, offset, block_size, n_records, flags, max_read_len, max_cigar_len))

    with open(dat_path, "rb") as df:
        magic = df.read(len(GLOBAL_MAGIC))
        if magic != GLOBAL_MAGIC:
            raise RuntimeError("Invalid .dat magic/version")
        majors, minors, dat_count, _ = GLOBAL_VER_PACK.unpack(df.read(GLOBAL_VER_PACK.size))
        if dat_count != len(entries):
            print(f"[warn] .dat block_count ({dat_count}) != .idx count ({len(entries)})")

        block_infos = {}
        wanted_nodes = set()
        for nid, offset, block_size, n_records, flags, R, C in entries:
            df.seek(offset, os.SEEK_SET)
            hdr = df.read(BLOCK_HDR_SIZE)
            if len(hdr) != BLOCK_HDR_SIZE:
                raise RuntimeError(f"Corrupt .dat: cannot read block header at {offset}")
            nid2, nrec2, flg2, R2, C2 = struct.unpack("<I I H I I", hdr)
            if (nid2 != nid) or (nrec2 != n_records):
                print(f"[warn] .dat/.idx mismatch for node {nid} (idx n={n_records}, dat n={nrec2})")
            # prefer values from dat
            R = R2 or R or 1
            C = C2 or C or 1
            rec_sz = record_size(R, C)
            block_infos[nid] = {
                "offset": offset,
                "n_records": n_records,
                "current_pos": 0,
                "max_read_len": R,
                "max_cigar_len": C,
                "record_size": rec_sz,
                "block_size": block_size,
            }
            wanted_nodes.add(nid)

    print(f"Reusing existing output with {len(block_infos)} node blocks from {idx_path}")
    return block_infos, dat_path, wanted_nodes

# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter, use_existing):
    if use_existing:
        print("Reusing existing .dat/.idx...")
        block_infos, dat_path, wanted_nodes = load_existing_output_files(output_prefix)
    else:
        print("Initializing output files from PKL maxima...")
        block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)
        print(f"Output file created: {dat_path}")

    next_milestone = milestone_step
    total_reads = 0
    start_time = time.perf_counter()

    # Open for in-place updates (do NOT truncate)
    dat_fh = open(dat_path, "r+b")

    # # cache of per-node packers to avoid rebuilding structs
    # packers = {}

    def write_segments_now(nid, segs):
        """Write a small batch (this read’s segments for a node) immediately."""
        if not segs:
            return

        info = block_infos[nid]
        R, C   = info["max_read_len"], info["max_cigar_len"]
        rec_sz = info["record_size"]

        # bounds check: don’t overrun block
        if info["current_pos"] + len(segs) > info["n_records"]:
            raise RuntimeError(
                f"Block overflow for node {nid}: current_pos={info['current_pos']}, "
                f"len(segs)={len(segs)}, n_records={info['n_records']}"
            )

        rec_pack = make_record_struct(R, C)
        # rec_pack = packers.get(nid)
        # if rec_pack is None:
        #     rec_pack = make_record_struct(R, C)
        #     packers[nid] = rec_pack

        # Build one small contiguous batch and write once
        n = len(segs)
        batch = bytearray(n * rec_sz)
        off = 0
        for s in segs:
            batch[off:off+rec_sz] = rec_pack.pack(
                int(s.offset),
                s.seq.ljust(R, b'\x00')[:R],
                s.bq.ljust(R, b'\x00')[:R],
                s.cigar.ljust(C, b'\x00')[:C],
                int(s.rq),
                s.strand if s.strand in (b'+', b'-') else b'+'
            )
            off += rec_sz

        base_offset = info["offset"] + BLOCK_HDR_SIZE
        pos = base_offset + info["current_pos"] * rec_sz
        dat_fh.seek(pos, os.SEEK_SET)
        dat_fh.write(batch)
        info["current_pos"] += n

    # Main loop: streaming; no global buffer
    for raw_msg in gam_record_iter(gam_path):
        segs_by_node = process_alignment(raw_msg, wanted_nodes, chrom_filter)
        total_reads += 1

        for nid, segs in segs_by_node.items():
            write_segments_now(nid, segs)

        if total_reads >= next_milestone:
            elapsed = time.perf_counter() - start_time
            print(f"{total_reads} reads processed | {elapsed:.1f} seconds")
            next_milestone += milestone_step

    dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Total reads processed: {total_reads}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="GAM segment extractor with per-block maxima read directly from stats PKL (streaming writes, sparse preallocation)."
    )
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the stats pickle (must contain max_read_length/max_cigar_length per node)")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--milestone", type=int, default=1_000_000, help="Progress report interval")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    parser.add_argument("--use-existing", action="store_true",
                        help="Reuse existing initialized output (output_prefix.dat/.idx) in NEW format")
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

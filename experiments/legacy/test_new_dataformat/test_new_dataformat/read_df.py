#!/usr/bin/env python3
import struct
import argparse
import sys
import os

# ─────────────────────────────────────────────────────────────────────────────
# Block headers (new + old)
# Latest (per-block maxima): <I I H I I>  (node_id, n_records, flags, max_read_len, max_cigar_len) → 18B
BLOCK_HDR_PACK_LATEST = struct.Struct("<I I H I I")   # 18 bytes
# Old (padded) node_length: <I I H 2x I> → 16B
BLOCK_HDR_PACK_PADDED = struct.Struct("<I I H 2x I")  # 16 bytes
# Old (compact) node_length: <I I H I> → 14B
BLOCK_HDR_PACK_14B    = struct.Struct("<I I H I")     # 14 bytes

def make_record_struct_old(node_length: int) -> struct.Struct:
    """
    OLD format: per-read record sized by node_length:
      <h {L}s {L}s {L}s h c>
        - i16 offset
        - seq[L] bytes
        - bq[L] bytes
        - cigar[L] bytes
        - i16 rq (MAPQ)
        - char strand
    """
    return struct.Struct(f"<h{node_length}s{node_length}s{node_length}shc")

def make_record_struct_latest(max_read_len: int, max_cigar_len: int) -> struct.Struct:
    """
    LATEST format: per-read record sized by per-block maxima:
      <h {R}s {R}s {C}s h c>
        - i16 offset
        - seq[R] bytes
        - bq[R] bytes
        - cigar[C] bytes
        - i16 rq (MAPQ)
        - char strand
    """
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")

# ─────────────────────────────────────────────────────────────────────────────
def load_index(idx_path):
    """
    Supports:
      • Latest fixed-size entry: 30B <I Q I I H I I>  (adds max_read_len, max_cigar_len)
      • Newer fixed-size entry: 26B <I Q I I H I>    (has node_length)
      • Older fixed-size entry: 22B <I Q I I H>      (no length fields)
      • Legacy variable-size with metadata_len per entry: <I Q I I H> + metadata_len bytes

    Returns: { node_id: {"start": offset, "size": block_size, "n_records": n_records,
                         "flags": flags, "node_length": node_len_or_None,
                         "max_read_len": R_or_None, "max_cigar_len": C_or_None } }
    """
    node_index = {}

    with open(idx_path, "rb") as f:
        header = f.read(4)
        if len(header) < 4:
            raise ValueError(f"Index file too small: {idx_path}")
        blocks_num, = struct.unpack("<I", header)

        file_size = os.fstat(f.fileno()).st_size
        remaining = file_size - 4

        def read_fixed_30():
            rec = f.read(30)
            if len(rec) != 30:
                raise EOFError("Unexpected EOF in 30B index entry")
            node_id, start, size, nrec, flags, max_read_len, max_cigar_len = struct.unpack("<I Q I I H I I", rec)
            return node_id, start, size, nrec, flags, None, max_read_len, max_cigar_len

        def read_fixed_26():
            rec = f.read(26)
            if len(rec) != 26:
                raise EOFError("Unexpected EOF in 26B index entry")
            node_id, start, size, nrec, flags, node_len = struct.unpack("<I Q I I H I", rec)
            return node_id, start, size, nrec, flags, node_len, None, None

        def read_fixed_22():
            rec = f.read(22)
            if len(rec) != 22:
                raise EOFError("Unexpected EOF in 22B index entry")
            node_id, start, size, nrec, flags = struct.unpack("<I Q I I H", rec)
            return node_id, start, size, nrec, flags, None, None, None

        def read_legacy_with_meta():
            hdr = f.read(22)
            if len(hdr) != 22:
                raise EOFError("Unexpected EOF in legacy index header")
            node_id, start, size, nrec, meta_len = struct.unpack("<I Q I I H", hdr)
            if meta_len:
                skipped = f.read(meta_len)
                if len(skipped) != meta_len:
                    raise EOFError("Unexpected EOF while skipping metadata")
            return node_id, start, size, nrec, 0, None, None, None

        # Pick strategy by exact division (prefer latest 30B)
        if blocks_num > 0 and remaining % blocks_num == 0:
            entry_size = remaining // blocks_num
            if entry_size == 30:
                strategy = "fixed30"
            elif entry_size == 26:
                strategy = "fixed26"
            elif entry_size == 22:
                strategy = "fixed22"
            else:
                strategy = "legacy"
        else:
            strategy = "legacy"

        for _ in range(blocks_num):
            if strategy == "fixed30":
                node_id, start, size, nrec, flags, node_len, R, C = read_fixed_30()
            elif strategy == "fixed26":
                node_id, start, size, nrec, flags, node_len, R, C = read_fixed_26()
            elif strategy == "fixed22":
                node_id, start, size, nrec, flags, node_len, R, C = read_fixed_22()
            else:
                node_id, start, size, nrec, flags, node_len, R, C = read_legacy_with_meta()

            node_index[node_id] = {
                "start": start,
                "size": size,
                "n_records": nrec,
                "flags": flags,
                "node_length": node_len,   # may be None
                "max_read_len": R,         # may be None
                "max_cigar_len": C,        # may be None
            }

    return node_index

# ─────────────────────────────────────────────────────────────────────────────
def _read_block_header(dat_file, block_start):
    """
    Read .dat block header at block_start.

    Tries (in order):
      1) LATEST 18B: <I I H I I> → returns ('latest', (R,C))
      2) OLD PADDED 16B: <I I H 2x I> → returns ('old', L)
      3) OLD 14B: <I I H I> → returns ('old', L)

    Returns: (node_id, n_records, flags, lengths, header_size_bytes, fmt)
      where lengths = (R, C) if fmt=='latest' else L (int)
    """
    # Try latest 18B
    dat_file.seek(block_start, os.SEEK_SET)
    hdr = dat_file.read(BLOCK_HDR_PACK_LATEST.size)
    if len(hdr) == BLOCK_HDR_PACK_LATEST.size:
        nid, nrec, flags, R, C = BLOCK_HDR_PACK_LATEST.unpack(hdr)
        if (1 <= R <= 1_000_000) and (1 <= C <= 1_000_000) and (0 <= nrec < 10_000_000):
            return nid, nrec, flags, (R, C), BLOCK_HDR_PACK_LATEST.size, 'latest'

    # Try old 16B padded
    dat_file.seek(block_start, os.SEEK_SET)
    hdr = dat_file.read(BLOCK_HDR_PACK_PADDED.size)
    if len(hdr) == BLOCK_HDR_PACK_PADDED.size:
        nid, nrec, flags, L = BLOCK_HDR_PACK_PADDED.unpack(hdr)
        if (1 <= L <= 1_000_000) and (0 <= nrec < 10_000_000):
            return nid, nrec, flags, L, BLOCK_HDR_PACK_PADDED.size, 'old'

    # Fallback: old 14B
    dat_file.seek(block_start, os.SEEK_SET)
    hdr = dat_file.read(BLOCK_HDR_PACK_14B.size)
    if len(hdr) != BLOCK_HDR_PACK_14B.size:
        raise RuntimeError(f"Cannot read block header at offset {block_start}")
    nid, nrec, flags, L = BLOCK_HDR_PACK_14B.unpack(hdr)
    if not (1 <= L <= 1_000_000 and 0 <= nrec < 10_000_000):
        raise RuntimeError(f"Suspicious block header at {block_start}: node_len={L}, nrec={nrec}")
    return nid, nrec, flags, L, BLOCK_HDR_PACK_14B.size, 'old'

# ─────────────────────────────────────────────────────────────────────────────
def read_segments(dat_path, node_index, node_id):
    if node_id not in node_index:
        raise ValueError(f"Node ID {node_id} not found in index")

    info = node_index[node_id]
    block_start = info["start"]

    records = []

    with open(dat_path, "rb") as f:
        # Read header from .dat to get lengths and true n_records
        nid, n_records, flags, lengths, hdr_size, fmt = _read_block_header(f, block_start)
        if nid != node_id:
            # Tolerate mismatch but continue
            pass

        if fmt == 'latest':
            R, C = lengths
            rec_struct = make_record_struct_latest(int(R), int(C))
        else:
            L = int(lengths)
            rec_struct = make_record_struct_old(L)

        rec_size = rec_struct.size

        # Jump to records start
        block_content_start = block_start + hdr_size
        f.seek(block_content_start, os.SEEK_SET)
        for _ in range(n_records):
            raw = f.read(rec_size)
            # print(raw)
            if len(raw) < rec_size:
                break
            offset, seq, bq, cigar_bytes, rq, strand = rec_struct.unpack(raw)

            record = {
                "offset": int(offset),
                "sequence": seq.rstrip(b'\x00').decode('ascii', errors='ignore'),
                "base_quality": list(bq.rstrip(b'\x00')),
                "cigar": cigar_bytes.rstrip(b'\x00').decode('ascii', errors='ignore'),
                "mapping_quality": int(rq),
                "strand": (strand.decode('ascii', errors='ignore')
                           if isinstance(strand, (bytes, bytearray)) else chr(strand)),
            }
            # print(record)
            records.append(record)

    return records

# ─────────────────────────────────────────────────────────────────────────────
def print_node_ids(node_index, n):
    node_ids = list(node_index.keys())
    node_ids.sort()
    print(f"Total nodes: {len(node_ids)}")
    print("Listing first {} node IDs:".format(n if n != -1 else len(node_ids)))
    for nid in node_ids[:n if n != -1 else None]:
        print(nid)

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Read segments or print node IDs from binary .dat/.idx files (supports latest and older formats)."
    )
    parser.add_argument("dat_path", help="Path to .dat file")
    parser.add_argument("idx_path", help="Path to .idx file")
    parser.add_argument("--node_id", type=int, help="Node ID to fetch")
    parser.add_argument("--topn", type=int, default=20, help="Show top N records (default 20, use -1 to show all)")
    parser.add_argument("--print-nodes", type=int, help="Print N node IDs from the index and exit")
    args = parser.parse_args()

    # Load index regardless of operation
    node_index = load_index(args.idx_path)

    # If printing node IDs only
    if args.print_nodes is not None:
        print_node_ids(node_index, args.print_nodes)
        return

    # node_id is required if not using --print-nodes
    if args.node_id is None:
        print("Error: --node_id is required unless --print-nodes is specified.", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    # Proceed to read segments
    records = read_segments(args.dat_path, node_index, args.node_id)
    print(f"Node {args.node_id} has {len(records)} records.")

    if args.topn == -1:
        for r in records:
            print(r)
    else:
        for r in records[:args.topn]:
            print(r)

if __name__ == "__main__":
    main()

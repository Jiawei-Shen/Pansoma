#!/usr/bin/env python3
import argparse
import gzip
import pickle
import struct
import time
import gc
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import itertools

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

# Globals in worker processes
WANTED_NODES = None
CHROM_FILTER = None

def _worker_init(wanted_nodes, chrom_filter):
    global WANTED_NODES, CHROM_FILTER
    WANTED_NODES = set(wanted_nodes)
    CHROM_FILTER = chrom_filter

# ─────────────────────────────────────────────────────────────────────────────
def read_varint(stream):
    value = 0
    shift_amount = 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError("EOF while reading varint")
        byte_value = b[0]
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
                # Stream-safe skip (works for gzip and plain files)
                for _ in range(group_count - 1):
                    skip_len = read_varint(f)
                    _ = f.read(skip_len)
                continue
            for _ in range(group_count - 1):
                try:
                    msg_size = read_varint(f)
                    yield f.read(msg_size)
                except EOFError:
                    break

# ─────────────────────────────────────────────────────────────────────────────
def _build_cigar(mapping_edits):
    parts = []
    for edit in mapping_edits:
        from_len = edit.from_length
        to_len = edit.to_length
        edit_len = len(edit.sequence)
        if from_len == to_len:
            parts.append(f"{from_len}M" if edit_len == 0 else f"{from_len}X")
        elif from_len > 0 and to_len == 0:
            parts.append(f"{from_len}D")
        elif from_len == 0 and to_len > 0:
            parts.append(f"{to_len}I")
        else:
            raise ValueError(f"Unexpected edit: from_length={from_len}, sequence_length={to_len}")
    return "".join(parts)

def _process_alignment_to_packed(raw_message):
    out = {}
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw_message)

    if alignment.mapping_quality <= 10:
        return out

    if CHROM_FILTER and not any(pos.name == CHROM_FILTER for pos in alignment.refpos):
        return out

    read_sequence = alignment.sequence
    read_quality = alignment.quality
    mapping_quality = alignment.mapping_quality
    read_offset = 0

    for mapping in alignment.path.mapping:
        node_id = mapping.position.node_id

        if node_id not in WANTED_NODES:
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

            if from_len == to_len:
                cigar_parts.append(f"{from_len}M" if edit_len == 0 else f"{from_len}X")
            elif from_len > 0 and to_len == 0:
                cigar_parts.append(f"{from_len}D")
            elif from_len == 0 and to_len > 0:
                cigar_parts.append(f"{to_len}I")
            else:
                raise ValueError(f"Unexpected edit: from_length={from_len}, sequence_length={edit_len}")

            # edit_length = max(from_len, edit_len)
            edit_length = to_len
            sequence_fragment = read_sequence[read_offset: read_offset + edit_length]
            quality_fragment = read_quality[read_offset: read_offset + edit_length]
            sequence_parts.append(sequence_fragment.upper())
            quality_parts.extend(quality_fragment)
            read_offset += edit_length

        cigar_string = "".join(cigar_parts)
        seq_final = "".join(sequence_parts).encode().ljust(150, b'\x00')[:150]
        bq_final = bytes(quality_parts).ljust(150, b'\x00')[:150]
        cigar_bytes = cigar_string.encode().ljust(20, b'\x00')[:20]

        rec = RECORD_STRUCT.pack(
            node_offset,
            seq_final,
            bq_final,
            cigar_bytes,
            mapping_quality,
            strand_char
        )

        buf = out.get(node_id)
        if buf is None:
            buf = bytearray()
            out[node_id] = buf
        buf += rec

    return out

def _process_batch(raw_messages):
    combined = {}
    for msg in raw_messages:
        d = _process_alignment_to_packed(msg)
        for nid, blob in d.items():
            if nid in combined:
                combined[nid] += blob
            else:
                combined[nid] = bytearray(blob)
    return combined, len(raw_messages)

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
            block_infos[node_id] = {"offset": current_offset, "n_records": n_records, "current_pos": 0}
            current_offset += 4 + 4 + 2 + n_records * RECORD_SIZE

    print(f"Filtered {len(wanted_nodes)} nodes from {total_nodes} total nodes "
          f"({len(wanted_nodes) / total_nodes:.2%} selected).")
    del stats_data
    gc.collect()

    dat_path = output_prefix + ".dat"
    with open(dat_path, "wb") as f:
        f.write(b"MYFMT\1")
        f.write(struct.pack("<BBI16s", 0, 1, len(block_infos), b'\x00' * 16))
        blank = RECORD_STRUCT.pack(0, b'\x00'*150, b'\x00'*150, b'\x00'*20, 0, b'+')
        for node_id, info in block_infos.items():
            header = struct.pack("<I I H", node_id, info["n_records"], 0)
            f.write(header + blank * info["n_records"])

    idx_path = output_prefix + ".idx"
    with open(idx_path, "wb") as idx_file:
        idx_file.write(struct.pack("<I", len(block_infos)))
        for node_id, info in block_infos.items():
            idx_file.write(struct.pack(
                "<I Q I I H",
                node_id,
                info["offset"],
                4 + 4 + 2 + info["n_records"] * RECORD_SIZE,
                info["n_records"],
                0
            ))

    return block_infos, dat_path, wanted_nodes

# ─────────────────────────────────────────────────────────────────────────────
def _batched(iterable, n):
    it = iter(iterable)
    while True:
        batch = list(itertools.islice(it, n))
        if not batch:
            return
        yield batch

# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(gam_path, stats_path, output_prefix, milestone_step, chrom_filter,
                 num_workers=4, batch_size=25000, buffer_segments=40_000_000, max_inflight=None):
    print("Initializing output files...")
    block_infos, dat_path, wanted_nodes = initialize_output_files(stats_path, output_prefix)
    print(f"Output file created: {dat_path}")

    processed_reads = 0
    total_dispatched = 0
    next_milestone = milestone_step
    start_time = time.perf_counter()

    dat_fh = open(dat_path, "r+b")
    node_buffers = defaultdict(bytearray)
    buffered_segments = 0

    def flush_node_buffers():
        nonlocal buffered_segments
        for node_id, blob in node_buffers.items():
            if not blob:
                continue
            info = block_infos.get(node_id)
            if info is None:
                continue
            base = info["offset"] + 4 + 4 + 2
            n_new = len(blob) // RECORD_SIZE
            # Optional safety guard against overflow:
            if info["current_pos"] + n_new > info["n_records"]:
                raise RuntimeError(f"Node {node_id} overflow: {info['current_pos']} + {n_new} > {info['n_records']}")
            pos = base + info["current_pos"] * RECORD_SIZE
            dat_fh.seek(pos)
            dat_fh.write(blob)
            info["current_pos"] += n_new
        node_buffers.clear()
        buffered_segments = 0

    if max_inflight is None:
        max_inflight = max(1, 2 * num_workers)

    with ProcessPoolExecutor(max_workers=num_workers,
                             initializer=_worker_init,
                             initargs=(wanted_nodes, chrom_filter)) as pool:
        futures = []

        for batch in _batched(gam_record_iter(gam_path), batch_size):
            # throttle in-flight tasks
            while len(futures) >= max_inflight:
                drained = False
                for fu in list(futures):
                    if fu.done():
                        futures.remove(fu)
                        (result, n_done) = fu.result()
                        for nid, blob in result.items():
                            node_buffers[nid] += blob
                            buffered_segments += (len(blob) // RECORD_SIZE)
                        processed_reads += n_done
                        drained = True
                        # processed-progress
                        if processed_reads >= next_milestone:
                            elapsed = time.perf_counter() - start_time
                            print(f"{processed_reads} reads processed | {elapsed:.1f} seconds")
                            next_milestone += milestone_step
                        break
                if not drained:
                    time.sleep(0.01)

            futures.append(pool.submit(_process_batch, batch))
            total_dispatched += len(batch)

            # opportunistic drain
            for fu in list(futures):
                if fu.done():
                    futures.remove(fu)
                    (result, n_done) = fu.result()
                    for nid, blob in result.items():
                        node_buffers[nid] += blob
                        buffered_segments += (len(blob) // RECORD_SIZE)
                    processed_reads += n_done
                    if processed_reads >= next_milestone:
                        elapsed = time.perf_counter() - start_time
                        print(f"{processed_reads} reads processed | {elapsed:.1f} seconds")
                        next_milestone += milestone_step

            if buffered_segments >= buffer_segments:
                flush_node_buffers()

        # final drain
        for fu in as_completed(futures):
            (result, n_done) = fu.result()
            for nid, blob in result.items():
                node_buffers[nid] += blob
                buffered_segments += (len(blob) // RECORD_SIZE)
            processed_reads += n_done
            if processed_reads >= next_milestone:
                elapsed = time.perf_counter() - start_time
                print(f"{processed_reads} reads processed | {elapsed:.1f} seconds")
                next_milestone += milestone_step

    flush_node_buffers()
    dat_fh.close()

    elapsed = time.perf_counter() - start_time
    print("\nFinal Summary:")
    print(f"  Reads processed: {processed_reads}")
    print(f"  Reads dispatched: {total_dispatched}")
    print(f"  Nodes included: {len(block_infos)}")
    print(f"  Elapsed time: {elapsed:.2f} seconds")

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="GAM segment extractor with CIGAR generation (parallel, processed-progress).")
    parser.add_argument("gam_path", help="Path to the GAM file")
    parser.add_argument("stats_pickle", help="Path to the node stats pickle file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--milestone", type=int, default=10_000_000, help="Progress report interval (#reads processed)")
    parser.add_argument("--chr", default="", help="Optional chromosome name to filter on")
    parser.add_argument("--workers", type=int, default=4, help="Number of worker processes")
    parser.add_argument("--batch", type=int, default=25000, help="GAM messages per worker task")
    parser.add_argument("--buffer_segments", type=int, default=10_000_000,
                        help="Flush threshold (#records) for writer buffers")
    parser.add_argument("--max_inflight", type=int, default=None,
                        help="Max number of in-flight worker tasks (default = 2*workers)")
    args = parser.parse_args()

    run_pipeline(
        gam_path=args.gam_path,
        stats_path=args.stats_pickle,
        output_prefix=args.output_prefix,
        milestone_step=args.milestone,
        chrom_filter=args.chr if args.chr else None,
        num_workers=args.workers,
        batch_size=args.batch,
        buffer_segments=args.buffer_segments,
        max_inflight=args.max_inflight
    )

if __name__ == "__main__":
    main()

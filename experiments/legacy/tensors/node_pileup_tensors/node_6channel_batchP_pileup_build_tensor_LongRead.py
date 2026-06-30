#!/usr/bin/env python3
"""
Build tensors for .idx/.dat nodes in memory-safe waves.

Key memory fixes:
  • Stream the .idx file into waves; don't build an all-nodes dict.
  • For each wave, stream merged_json and pick only needed node records.
  • Publish only wave-scoped GLOBAL_NODE_SEQS / GLOBAL_NODE_AF_BINS before forking.

Channel 6 policy (AF):
  - Integer bins 0..7 per base.
  - If merged.json AF is digit string '0'..'7' per base -> use digits as bins.
  - If AF is a float list -> bin into 8 levels:
      0–1e-6→0, 1e-6–1e-5→1, 1e-5–1e-4→2, 1e-4–1e-3→3,
      1e-3–1e-2→4, 1e-2–0.1→5, 0.1–0.5→6, 0.5–1.0→7
  - If a node has NO AF -> channel 6 is 0 everywhere, except the variant column is 3.

Usage example:
  python build_all_nodes_from_idx_waves.py DAT.dat IDX.idx OUT_DIR merged.json[.gz] \
    --num_workers 32 --chunksize 512 --wave_size 100000 --variant_type snp
"""

import argparse
import gzip
import json
import math
import os
import re
import struct
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Optional, Tuple, Iterable, Set

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Logging

def ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def log(msg: str) -> None:
    print(f"[{ts()}] {msg}", file=sys.stderr, flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# Formats / constants

RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size

# NEW FORMAT ADDITIONS (kept minimal; rest of code unchanged)
# Per-node .dat block header: <I I H I> = node_id (u32), n_records (u32), flags (u16), node_length (u32)
BLOCK_HDR_PACK = struct.Struct("<I I H I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size  # 14 bytes

def make_record_struct(node_length: int) -> struct.Struct:
    # Per-node variable-size record
    # <h {L}s {L}s 30s h c>
    return struct.Struct(f"<h{node_length}s{node_length}s30shc")

BASE_TO_INDEX = {'A': 20, 'C': 30, 'G': 50, 'T': 70, 'N': 10, '*': 90, '_PADDING_': 0}
PADDING_BASE_INDEX = 0

CIGAR_OP_TO_INDEX = {'M': 10, 'N': 20, 'S': 30, 'I': 40, 'D': 50, 'H': 60, 'P': 70, '=': 80, 'X': 90, '_PADDING_': 0}
CIGAR_PADDING_INDEX = 0

INDEX_TO_BASE_FOR_VIEW = {20: 'A', 30: 'C', 50: 'G', 70: 'T', 10: 'N', 90: '*', 0: '0'}

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
DEFAULT_QUALITY_PADDING = 0
DEFAULT_MAPPING_QUALITY_PADDING = -1
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1

# Worker globals (wave-scoped)
worker_dat_file = None
worker_base_output_dir = None
worker_need_view = False
GLOBAL_NODE_SEQS: Dict[int, str] = {}
GLOBAL_NODE_AF_BINS: Dict[int, Optional[np.ndarray]] = {}  # None -> missing AF

# ─────────────────────────────────────────────────────────────────────────────
# Helpers

def reverse_complement(sequence: str) -> str:
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def calculate_window_start(variant_pos: int, window_size: int) -> int:
    center_index = window_size // 2
    return variant_pos - center_index

def decode_cigar_to_int_ops(cigar_string: str) -> List[Tuple[int, str]]:
    if not cigar_string or cigar_string == '*':
        return []
    try:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse CIGAR '{cigar_string}': {e}\n")
        return []

# ─────────────────────────────────────────────────────────────────────────────
# Stream the IDX in waves (no giant dict)

def iter_idx_waves(idx_path: str, wave_size: int) -> Iterable[List[Tuple[int, int, int]]]:
    """Yield lists of (node_id, offset, n_records) of size up to wave_size.

    Supports old 22B entries (<I Q I I H>) and new 26B entries (<I Q I I H I>).
    """
    with open(idx_path, 'rb') as f:
        header = f.read(4)
        if len(header) < 4:
            raise RuntimeError("Index file too small")
        n = struct.unpack('<I', header)[0]

        # Determine entry size from remaining bytes
        f.seek(0, os.SEEK_END)
        remaining = f.tell() - 4
        f.seek(4, os.SEEK_SET)
        entry_size = remaining // max(n, 1)
        if entry_size == 26:
            unpack_fmt = '<I Q I I H I'
        elif entry_size == 22:
            unpack_fmt = '<I Q I I H'
        else:
            # default to new format
            entry_size = 26
            unpack_fmt = '<I Q I I H I'

        wave: List[Tuple[int, int, int]] = []
        for i in range(n):
            rec = f.read(entry_size)
            if len(rec) < entry_size:
                log(f"Index ended prematurely at record {i+1}")
                break
            if entry_size == 26:
                node_id, offset, _block_size, n_records, _flags, _node_len = struct.unpack(unpack_fmt, rec)
            else:
                node_id, offset, _block_size, n_records, _flags = struct.unpack(unpack_fmt, rec)
            wave.append((node_id, offset, n_records))
            if len(wave) == wave_size:
                wave.sort(key=lambda t: t[1])  # by dat offset
                yield wave
                wave = []
        if wave:
            wave.sort(key=lambda t: t[1])
            yield wave

# ─────────────────────────────────────────────────────────────────────────────
# Streaming JSON reader: yield node dicts without loading whole file

def _open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def iter_nodes_from_merged(merged_path: str) -> Iterable[dict]:
    """
    Stream node objects from merged JSON that is either:
      • a top-level list: [ {node}, {node}, ... ]
      • or a dict wrapper: {"nodes":[ {node}, ... ], ...}
    We do a lightweight brace-tracking parse to extract each { ... } inside the first '[' array.
    """
    with _open_maybe_gzip(merged_path, "rt") as f:
        buf = []
        in_array = False
        depth = 0
        in_str = False
        esc = False

        # Scan to the first '['
        while True:
            ch = f.read(1)
            if not ch:
                return
            if ch == '[':
                in_array = True
                break

        # Now extract JSON objects delimited by balanced '{...}' up to matching ']'
        while True:
            ch = f.read(1)
            if not ch:
                break

            if not in_array:
                continue

            # Skip whitespace and commas until an object starts
            if depth == 0:
                if ch.isspace() or ch == ',':
                    continue
                if ch == ']':
                    break
                if ch != '{':
                    # Unexpected token; keep scanning
                    continue
                # Start of object
                buf = ['{']
                depth = 1
                in_str = False
                esc = False
                continue

            # We're inside an object
            buf.append(ch)

            if in_str:
                if esc:
                    esc = False
                elif ch == '\\':
                    esc = True
                elif ch == '"':
                    in_str = False
            else:
                if ch == '"':
                    in_str = True
                elif ch == '{':
                    depth += 1
                elif ch == '}':
                    depth -= 1
                    if depth == 0:
                        # end of object
                        obj_str = ''.join(buf)
                        try:
                            yield json.loads(obj_str)
                        except Exception as e:
                            sys.stderr.write(f"Warning: failed to parse node object: {e}\n")
                        buf = []
    # done

# ─────────────────────────────────────────────────────────────────────────────
# AF binning

def af_float_to_bin(x: float) -> int:
    try:
        x = float(x)
    except Exception:
        return 0
    if x <= 0.0:          return 0
    if x < 1e-6:          return 0
    if x < 1e-5:          return 1
    if x < 1e-4:          return 2
    if x < 1e-3:          return 3
    if x < 1e-2:          return 4
    if x < 0.1:           return 5
    if x < 0.5:           return 6
    return 7

def normalize_af_to_bins_compact(af_field, expected_len: int) -> Optional[np.ndarray]:
    """
    Return np.uint8 array length==expected_len with values 0..7, or None if AF truly missing.
    Accept list[float] or digit string "0".."7" per base.
    """
    if af_field is None:
        return None
    if isinstance(af_field, str):
        if not af_field:
            return None
        # Map ascii '0'..'7' -> 0..7; other chars -> 0
        a = np.frombuffer(af_field.encode('ascii'), dtype=np.uint8)
        a = np.where((a >= ord('0')) & (a <= ord('7')), a - ord('0'), 0).astype(np.uint8)
        if a.size < expected_len:
            out = np.zeros(expected_len, dtype=np.uint8)
            out[:a.size] = a
            return out
        return a[:expected_len].copy()
    if isinstance(af_field, list):
        # bin floats
        out = np.zeros(expected_len, dtype=np.uint8)
        upto = min(expected_len, len(af_field))
        # vectorizing via list-comp for clarity
        out[:upto] = np.fromiter((af_float_to_bin(v) for v in af_field[:upto]), count=upto, dtype=np.uint8)
        return out
    return None

# ─────────────────────────────────────────────────────────────────────────────
# Load only the current wave's sequences/AF from merged JSON (streaming)

def load_wave_seqs_and_af_bins_stream(merged_path: str, wanted_ids: Set[int]) -> Tuple[Dict[int, str], Dict[int, Optional[np.ndarray]], int]:
    seqs: Dict[int, str] = {}
    af_bins: Dict[int, Optional[np.ndarray]] = {}
    missing = 0
    wanted = set(wanted_ids)
    found = 0

    for rec in iter_nodes_from_merged(merged_path):
        if not isinstance(rec, dict):
            continue
        if 'sequence' not in rec:
            continue
        nid_raw = rec.get('node_id')
        if nid_raw is None:
            # allow records that only have sequence but not node_id? skip
            continue
        try:
            nid = int(str(nid_raw))
        except Exception:
            continue
        if nid not in wanted:
            continue
        seq = str(rec.get('sequence', '')).upper()
        if not seq:
            missing += 1
            continue
        seqs[nid] = seq
        af_bins[nid] = normalize_af_to_bins_compact(rec.get('genomead_af'), len(seq))
        found += 1
        if found == len(wanted):
            break

    return seqs, af_bins, missing

# ─────────────────────────────────────────────────────────────────────────────
# Variant helpers

def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_quality_values, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_allele_for_indel=None):
    current_node_pos = read_offset_on_node
    current_read_pos = 0
    for length, op in read_cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            if current_node_pos <= target_node_pos < current_node_pos + length:
                offset_in_block = target_node_pos - current_node_pos
                read_idx = current_read_pos + offset_in_block
                if read_idx < len(read_sequence):
                    allele = read_sequence[read_idx].upper()
                    quality = read_quality_values[read_idx] if read_idx < len(read_quality_values) else 0
                    if expected_var_type in ('I', 'D'):
                        return "REF_STATE_FOR_INDEL", quality
                    return allele, quality
                return None, None
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                if current_read_pos + length <= len(read_sequence):
                    qualities = read_quality_values[current_read_pos: current_read_pos + length]
                    mean_quality = sum(qualities) / len(qualities) if qualities else 0.0
                    return read_sequence[current_read_pos: current_read_pos + length].upper(), mean_quality
                return None, None
            current_read_pos += length
        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I':
                    return "OTHER_FOR_INDEL", None
                if expected_var_type == 'D':
                    if 0 <= current_node_pos < len(node_sequence) and current_node_pos + length <= len(node_sequence):
                        deleted_seq_in_ref_context = node_sequence[current_node_pos: current_node_pos + length]
                        if deleted_seq_in_ref_context == expected_ref_allele_for_indel:
                            return "*", None
                        else:
                            return "OTHER_FOR_INDEL", None
                    return "OTHER_FOR_INDEL", None
                return "*", None
            current_node_pos += length
        elif op == 'S':
            current_read_pos += length
        elif op == 'N':
            current_node_pos += length
        if current_node_pos > target_node_pos + 1 and not (
            expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):
            break
    return None, None

def detect_variants_from_cigar(offset_on_node, cigar_ops_decoded, read_sequence, node_sequence):
    variants = []
    node_pos, read_pos = offset_on_node, 0
    node_seq_len, read_seq_len = len(node_sequence), len(read_sequence)
    for length, op in cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            for i in range(length):
                cur_node_p, cur_read_p = node_pos + i, read_pos + i
                if cur_node_p < node_seq_len and cur_read_p < read_seq_len:
                    node_base = node_sequence[cur_node_p].upper()
                    read_base = read_sequence[cur_read_p].upper()
                    if node_base != read_base and op != '=':
                        variants.append((cur_node_p, 'X', read_base, node_base))
                else:
                    break
            node_pos += length
            read_pos += length
        elif op == 'I':
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < node_seq_len else "*"
            variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor))
            read_pos += length
        elif op == 'D':
            deleted_sequence_from_ref = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= node_seq_len else ""
            if deleted_sequence_from_ref:
                variants.append((node_pos, 'D', "*", deleted_sequence_from_ref))
            node_pos += length
        elif op == 'S':
            read_pos += length
        elif op == 'N':
            node_pos += length
    return variants

def get_read_representation_in_window_for_view(segment_cigar_ops, segment_offset_on_node, segment_read_sequence,
                                               window_start_node, window_size, node_len):
    window_chars = [' '] * window_size
    node_pos, read_pos = segment_offset_on_node, 0
    read_seq_len = len(segment_read_sequence)
    for L, op in segment_cigar_ops:
        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = node_pos + i, read_pos + i
                win_idx = n_aln - window_start_node
                if 0 <= win_idx < window_size and r_aln < read_seq_len:
                    window_chars[win_idx] = segment_read_sequence[r_aln].upper()
            node_pos += L
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                win_idx = (node_pos + i) - window_start_node
                if 0 <= win_idx < window_size:
                    window_chars[win_idx] = '*'
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L
        if node_pos >= window_start_node + window_size and read_pos > 0:
            break
        if read_pos >= read_seq_len:
            break
    return window_chars

def get_read_tensor_rows_in_window(segment_cigar_ops, segment_offset_on_node,
                                   segment_read_sequence, segment_quality_values,
                                   mapping_quality,
                                   window_start_node, tensor_win_size, node_len):
    bases = [PADDING_BASE_INDEX] * tensor_win_size
    quals = [DEFAULT_QUALITY_PADDING] * tensor_win_size
    mapqs = [DEFAULT_MAPPING_QUALITY_PADDING] * tensor_win_size
    cigar_ops_indices = [CIGAR_PADDING_INDEX] * tensor_win_size
    node_pos, read_pos = segment_offset_on_node, 0
    read_seq_len = len(segment_read_sequence)
    qual_len = len(segment_quality_values)
    for L, op in segment_cigar_ops:
        op_idx = CIGAR_OP_TO_INDEX.get(op, CIGAR_PADDING_INDEX)
        if node_pos >= window_start_node + tensor_win_size and read_pos > 0:
            break
        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = node_pos + i, read_pos + i
                win_idx = n_aln - window_start_node
                if 0 <= win_idx < tensor_win_size:
                    cigar_ops_indices[win_idx] = op_idx
                    mapqs[win_idx] = mapping_quality
                    if r_aln < read_seq_len:
                        base_char = segment_read_sequence[r_aln].upper()
                        bases[win_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                        if r_aln < qual_len:
                            quals[win_idx] = segment_quality_values[r_aln]
            node_pos += L
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                win_idx = (node_pos + i) - window_start_node
                if 0 <= win_idx < tensor_win_size:
                    cigar_ops_indices[win_idx] = op_idx
                    mapqs[win_idx] = mapping_quality
                    bases[win_idx] = BASE_TO_INDEX['*']
                    quals[win_idx] = DEFAULT_QUALITY_PADDING
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L
        if read_pos >= read_seq_len:
            break
    return bases, quals, mapqs, cigar_ops_indices

# ─────────────────────────────────────────────────────────────────────────────
# Worker init & main worker

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker, need_view_flag):
    global worker_dat_file, worker_base_output_dir, worker_need_view
    worker_dat_file = open(dat_file_path_for_worker, 'rb')
    worker_base_output_dir = base_output_dir_for_worker
    worker_need_view = bool(need_view_flag)

def process_single_node_for_pileup(task_args):
    (node_id, dat_file_offset, n_records,
     min_af_threshold, min_variants_threshold, min_allele_bq_threshold,
     variant_type_to_process) = task_args

    global worker_dat_file, worker_base_output_dir, worker_need_view
    global GLOBAL_NODE_SEQS, GLOBAL_NODE_AF_BINS

    tensor_files_generated_for_node = 0
    if worker_dat_file is None or worker_base_output_dir is None:
        return node_id, None, tensor_files_generated_for_node

    node_sequence = GLOBAL_NODE_SEQS.get(node_id, "")
    if not node_sequence:
        return node_id, {}, tensor_files_generated_for_node
    node_len = len(node_sequence)

    af_bins = GLOBAL_NODE_AF_BINS.get(node_id, None)

    aligned_read_segments = []
    try:
        # NEW FORMAT: read block header to get node_length; then build per-node record struct
        worker_dat_file.seek(dat_file_offset, os.SEEK_SET)
        block_hdr = worker_dat_file.read(BLOCK_HDR_SIZE)
        if len(block_hdr) != BLOCK_HDR_SIZE:
            raise RuntimeError(f"Cannot read block header at offset {dat_file_offset} for node {node_id}")
        nid2, nrec2, _flags, node_length_from_dat = BLOCK_HDR_PACK.unpack(block_hdr)
        # not fatal if mismatch; continue using idx-provided counts
        rec_struct = make_record_struct(int(node_length_from_dat))
        rec_size = rec_struct.size

        # Move to first record (after block header) and read this block's payload
        worker_dat_file.seek(dat_file_offset + BLOCK_HDR_SIZE, os.SEEK_SET)
        bulk = worker_dat_file.read(n_records * rec_size)
        if len(bulk) < n_records * rec_size:
            n_records = len(bulk) // rec_size
            bulk = bulk[: n_records * rec_size]

        for (off_from_file, raw_seq, raw_qual, raw_cigar, mapq_val, strand_byte) in rec_struct.iter_unpack(bulk):
            if mapq_val < 10:
                continue
            try:
                seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                qual_values = list(raw_qual.rstrip(b'\0'))
                cigar_str = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand_char = strand_byte.decode('ascii') if isinstance(strand_byte, (bytes, bytearray)) else chr(strand_byte)
            except UnicodeDecodeError:
                continue
            if not seq or len(seq) != len(qual_values):
                continue
            cigar_ops = decode_cigar_to_int_ops(cigar_str)
            if not cigar_ops and cigar_str != '*':
                continue

            cur_seq = seq
            cur_qual = qual_values
            cur_ops = cigar_ops
            cur_off = off_from_file
            if strand_char == '-':
                cur_seq = reverse_complement(seq)
                cur_qual = qual_values[::-1]
                cur_ops = [op for op in reversed(cigar_ops)] if cigar_ops else []
                # aln_span = len(cur_seq)
                aln_span = len(cur_seq)

                # Adjust span: +I, -D (reference-consuming vs non-consuming edits)
                for L, op in cigar_ops:
                    if op == 'I':
                        aln_span -= L
                    elif op == 'D':
                        aln_span += L
                cur_off = node_len - aln_span - off_from_file
                if cur_off < 0:
                    continue

            aligned_read_segments.append({
                "offset_on_node": cur_off,
                "read_sequence": cur_seq,
                "processed_quality_values": cur_qual,
                "cigar_ops": cur_ops,
                "original_cigar_str": cigar_str,
                "strand": strand_char,
                "mapping_quality": mapq_val
            })
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()} Node {node_id}]: {e}\n")
        return node_id, None, tensor_files_generated_for_node

    if not aligned_read_segments:
        return node_id, {}, tensor_files_generated_for_node

    print(node_id, node_len , len(aligned_read_segments))
    # guard pathological nodes
    if len(aligned_read_segments) > 10000:
        return node_id, None, tensor_files_generated_for_node

    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    view_oriented_variant_data = {} if worker_need_view else None
    variant_headers_for_summary = []

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
        if variant_type_to_process == 'snp' and v_type != 'X':
            continue
        if variant_type_to_process == 'indel' and v_type not in ('I', 'D'):
            continue

        alt_allele_count = ref_allele_count = other_allele_count = locus_coverage = 0
        alt_allele_bq = []

        expected_ref = v_ref_from_cigar
        expected_alt = v_alt_from_cigar
        ref_for_indel_ctx = None

        if v_type == 'D':
            expected_alt = "*"
            if 0 <= v_pos < node_len:
                expected_ref = node_sequence[v_pos]
            ref_for_indel_ctx = v_ref_from_cigar
        elif v_type == 'I':
            expected_ref = node_sequence[v_pos] if 0 <= v_pos < node_len else "*"
            ref_for_indel_ctx = expected_ref

        for seg in aligned_read_segments:
            allele, bq = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["processed_quality_values"], seg["cigar_ops"],
                v_pos, node_sequence, v_type, ref_for_indel_ctx)
            if allele is not None:
                locus_coverage += 1
                if allele == expected_alt:
                    alt_allele_count += 1
                    if bq is not None:
                        alt_allele_bq.append(bq)
                elif allele == expected_ref or (v_type in ('I', 'D') and allele == "REF_STATE_FOR_INDEL"):
                    ref_allele_count += 1
                else:
                    other_allele_count += 1

        if alt_allele_count < min_variants_threshold:
            continue

        if v_type == 'X':
            tmp_af = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
            if tmp_af < min_af_threshold:
                continue
            mean_bq_tmp = sum(alt_allele_bq) / len(alt_allele_bq) if alt_allele_bq else 0.0
            if mean_bq_tmp < min_allele_bq_threshold:
                continue

        current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
        mean_alt_bq = sum(alt_allele_bq) / len(alt_allele_bq) if alt_allele_bq else 0.0

        variant_key_string = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

        # Reference row
        ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i, pos in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
            if 0 <= pos < node_len:
                ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[pos].upper(), BASE_TO_INDEX['N'])

        # AF bins row
        af_row_bins = np.zeros(TENSOR_WINDOW_SIZE, dtype=np.uint8)
        if af_bins is not None:
            for i, pos in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                if 0 <= pos < node_len and pos < af_bins.shape[0]:
                    af_row_bins[i] = af_bins[pos]
        else:
            variant_window_index = v_pos - window_start_pos
            if 0 <= variant_window_index < TENSOR_WINDOW_SIZE:
                af_row_bins[variant_window_index] = 3

        # Allocate channels
        H, W = 1 + TENSOR_MAX_READ_ROWS, TENSOR_WINDOW_SIZE
        ch1 = np.full((H, W), PADDING_BASE_INDEX, dtype=np.int8)
        ch2 = np.full((H, W), DEFAULT_QUALITY_PADDING, dtype=np.int8)
        ch3 = np.full((H, W), MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
        ch4 = np.full((H, W), DEFAULT_MAPPING_QUALITY_PADDING, dtype=np.int8)
        ch5 = np.full((H, W), CIGAR_PADDING_INDEX, dtype=np.int8)
        ch6 = np.zeros((H, W), dtype=np.int8)

        ch1[0, :] = np.asarray(ref_base_indices_row, dtype=np.int8)
        ch6[0, :] = af_row_bins.astype(np.int8, copy=False)

        reads_added = 0
        variant_window_index = v_pos - window_start_pos
        for seg in aligned_read_segments:
            if reads_added >= TENSOR_MAX_READ_ROWS:
                break
            base_row, qual_row, mapq_row, cigar_row = get_read_tensor_rows_in_window(
                seg["cigar_ops"], seg["offset_on_node"],
                seg["read_sequence"], seg["processed_quality_values"],
                max(0, min(int(seg["mapping_quality"]), 127)),
                window_start_pos, TENSOR_WINDOW_SIZE, node_len)

            if any(b != PADDING_BASE_INDEX for b in base_row):
                r = 1 + reads_added
                ch1[r, :] = np.asarray(base_row, dtype=np.int8)
                ch2[r, :] = np.asarray(qual_row, dtype=np.int8)
                ref_row = ch1[0, :].astype(np.int16)
                read_row = ch1[r, :].astype(np.int16)
                flags = np.full(W, MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
                mask_valid = (read_row != PADDING_BASE_INDEX) & (ref_row != PADDING_BASE_INDEX)
                flags[mask_valid] = (read_row[mask_valid] != ref_row[mask_valid]).astype(np.int8)
                if 0 <= variant_window_index < W and mask_valid[variant_window_index] and flags[variant_window_index] == 1:
                    flags[variant_window_index] = 5
                ch3[r, :] = flags
                ch4[r, :] = np.asarray(mapq_row, dtype=np.int8)
                ch5[r, :] = np.asarray(cigar_row, dtype=np.int8)
                ch6[r, :] = ch6[0, :]
                reads_added += 1

        tensor = np.stack([ch1, ch2, ch3, ch4, ch5, ch6], axis=0)
        node_dir = os.path.join(worker_base_output_dir, str(node_id))
        os.makedirs(node_dir, exist_ok=True)
        npy_name = f"{variant_key_string}.npy"
        np.save(os.path.join(node_dir, npy_name), tensor)

        variant_headers_for_summary.append({
            "variant_key": variant_key_string, "tensor_file": npy_name,
            "alt_allele_count": alt_allele_count,
            "ref_allele_count_at_locus": ref_allele_count,
            "other_allele_count_at_locus": other_allele_count,
            "coverage_at_locus": locus_coverage,
            "alt_allele_frequency": round(current_alt_freq, 4),
            "mean_alt_allele_base_quality": round(mean_alt_bq, 2)
        })
        tensor_files_generated_for_node += 1

    if variant_headers_for_summary:
        with open(os.path.join(worker_base_output_dir, str(node_id), "variant_summary.json"), 'w') as f:
            json.dump({"node_id": node_id, "node_length": len(node_sequence),
                       "variants_passing_af_filter": variant_headers_for_summary}, f, indent=2)

    return node_id, (view_oriented_variant_data or {}), tensor_files_generated_for_node

# ─────────────────────────────────────────────────────────────────────────────
# Optional pileup viewer

def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if max_variants_to_display == 0 or not node_data_for_display_view:
        return
    print(f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))
    for i, variant_key in enumerate(sorted_variant_keys):
        if i >= max_variants_to_display:
            print(f"\n  ... ({len(sorted_variant_keys) - i} more variants not shown)")
            break
        v_pos, v_type = int(variant_key.split('_')[0]), variant_key.split('_')[1]
        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)
        print(f"\n--- Variant: {variant_key} ---")
        ref_chars = []
        for j in range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE):
            ref_chars.append(full_node_sequence[j] if 0 <= j < len(full_node_sequence) else '0')
        print(f"  Node Ref: {''.join(ref_chars)}")
        marker_pos = v_pos - window_start_pos
        marker_line = [' '] * TENSOR_WINDOW_SIZE
        if 0 <= marker_pos < TENSOR_WINDOW_SIZE:
            marker_line[marker_pos] = '^'
        print(f"  Marker  : {''.join(marker_line)}")

# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    ap = argparse.ArgumentParser(
        description="Build tensors for .idx nodes in waves; wave-scoped seq/AF loading to control memory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("dat", help=".dat alignment file")
    ap.add_argument("idx", help=".idx index file")
    ap.add_argument("output", help="Base output directory")
    ap.add_argument("merged_json", help="Merged node JSON with sequence and AF (list or {'nodes':[...]}); .gz ok")
    ap.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Worker processes")
    ap.add_argument("--chunksize", type=int, default=512, help="executor.map chunksize")
    ap.add_argument("--wave_size", type=int, default=100000, help="Max nodes to submit per wave")
    ap.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                    help="Print pileups for top N variants per node (-1 for all)")
    ap.add_argument("--max_view_reads", type=int, default=20, help="Reads per pileup in view mode")
    ap.add_argument("--min_af", type=float, default=0.1, help="Min AF for SNP variant to pass")
    ap.add_argument("--min_variants", type=int, default=3, help="Min alt count to pass")
    ap.add_argument("--min_allele_bq", type=float, default=10.0, help="Min mean base quality for alt")
    ap.add_argument("--variant_type", type=str, default='all', choices=['snp', 'indel', 'all'],
                    help="Which variants to emit")
    args = ap.parse_args()

    for p in (args.dat, args.idx, args.merged_json):
        if not os.path.isfile(p):
            sys.exit(f"Error: file not found: {p}")
    os.makedirs(args.output, exist_ok=True)

    need_view = args.view is not None
    wave_size = max(1, int(args.wave_size))

    log(f"Streaming IDX in waves of {wave_size:,} and loading merged JSON per-wave.")

    total_tasks = 0
    total_processed = 0
    total_tensors = 0
    wave_num = 0

    for wave in iter_idx_waves(args.idx, wave_size):
        wave_num += 1
        wave_ids = [nid for (nid, _, _) in wave]
        wave_id_set = set(wave_ids)

        # Load only this wave's seq/AF
        t0 = time.time()
        seqs_map, af_bins_map, missing_seq_recs = load_wave_seqs_and_af_bins_stream(args.merged_json, wave_id_set)
        load_dt = time.time() - t0

        # Build tasks for nodes present in merged_json
        wave_tasks = []
        for nid, off, nrec in wave:
            if nid in seqs_map:
                wave_tasks.append((nid, off, nrec, args.min_af, args.min_variants, args.min_allele_bq, args.variant_type))
        if not wave_tasks:
            log(f"[Wave {wave_num}] No nodes in this wave had sequences in merged_json; skipping.")
            continue

        # Publish wave-scoped globals BEFORE forking
        GLOBAL_NODE_SEQS.clear()
        GLOBAL_NODE_AF_BINS.clear()
        GLOBAL_NODE_SEQS.update(seqs_map)
        GLOBAL_NODE_AF_BINS.update(af_bins_map)

        log(f"[Wave {wave_num}] Submitting {len(wave_tasks):,} nodes "
            f"(wave JSON load {load_dt:.2f}s; missing-in-records={missing_seq_recs})")

        processed = 0
        tensors = 0
        batch_nodes = 0
        batch_tensors = 0
        t_batch = time.time()

        with ProcessPoolExecutor(max_workers=args.num_workers,
                                 initializer=init_worker,
                                 initargs=(args.dat, args.output, need_view)) as ex:
            for node_id, view_data, tensor_count in ex.map(
                    process_single_node_for_pileup, wave_tasks, chunksize=max(1, args.chunksize)):
                total_processed += 1
                processed += 1
                batch_nodes += 1
                total_tensors += tensor_count
                tensors += tensor_count
                batch_tensors += tensor_count

                if need_view and view_data:
                    display_pileup_data(view_data, str(node_id), GLOBAL_NODE_SEQS.get(node_id, ""),
                                        args.max_view_reads,
                                        args.view if args.view != -1 else float('inf'))

                if batch_nodes >= 100000 or processed == len(wave_tasks):
                    dt = time.time() - t_batch
                    rate = batch_nodes / dt if dt > 0 else 0.0
                    log(f"[Wave {wave_num}] processed {processed:,}/{len(wave_tasks):,}  "
                        f"batch_nodes={batch_nodes:,} tensors_in_batch={batch_tensors:,}  "
                        f"in {dt:.2f}s ({rate:.1f} nodes/s)")
                    batch_nodes = 0
                    batch_tensors = 0
                    t_batch = time.time()

        log(f"[Wave {wave_num}] DONE — nodes:{processed:,}, tensors:{tensors:,}")

        # Free wave data
        GLOBAL_NODE_SEQS.clear()
        GLOBAL_NODE_AF_BINS.clear()
        import gc; gc.collect()

        total_tasks += len(wave_tasks)

    log(f"ALL WAVES DONE. Nodes processed: {total_processed:,}/{total_tasks:,}; total tensors: {total_tensors:,}")

if __name__ == "__main__":
    main()

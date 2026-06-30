#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np
from collections import defaultdict
import re
from concurrent.futures import ProcessPoolExecutor

# ─────────────────────────────────────────────────────────────────────────────
# Global/format structs

# Global .dat header (writer used "MYFMT\x01" + <BBI16s>)
DAT_GLOBAL_MAGIC = b"MYFMT\x01"
GLOBAL_VER_PACK = struct.Struct("<BBI16s")  # major, minor, block_count, reserved[16]
# Old block header (node_length):
OLD_BLOCK_HDR_PACK = struct.Struct("<I I H I")       # nid, nrec, flags, node_length
OLD_BLOCK_HDR_SIZE = OLD_BLOCK_HDR_PACK.size         # 14
# New block header (per-block maxima):
NEW_BLOCK_HDR_PACK = struct.Struct("<I I H I I")     # nid, nrec, flags, max_read_len, max_cigar_len
NEW_BLOCK_HDR_SIZE = NEW_BLOCK_HDR_PACK.size         # 18

# (Kept for reference; we build per-node record structs dynamically.)
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size

BASE_TO_INDEX = {
    'A': 20, 'C': 30, 'G': 50, 'T': 70,
    'N': 5,      # CHANGED
    '*': 1,      # CHANGED
    'I': 90,     # ADDED: insertion placeholder token
    '_PADDING_': 0
}
PADDING_BASE_INDEX = 0

CIGAR_OP_TO_INDEX = {
    'M': 10, 'N': 20, 'S': 30, 'I': 40, 'D': 50, 'H': 60, 'P': 70, '=': 80, 'X': 90,
    '_PADDING_': 0
}
CIGAR_PADDING_INDEX = 0

INDEX_TO_BASE_FOR_VIEW = {
    20: 'A', 30: 'C', 50: 'G', 70: 'T',
    5: 'N',   # CHANGED to reflect 'N' -> 5
    1: '*',   # CHANGED to reflect '*' -> 1
    0: ' '
}

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
DEFAULT_QUALITY_PADDING = 0
DEFAULT_MAPPING_QUALITY_PADDING = -1
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1

# NEW: drop INDELs longer than this (bp)
MAX_INDEL_LEN = 250

# Worker globals
worker_dat_file = None
worker_base_output_dir = None
worker_need_view = False
worker_is_new_format = None  # True=new (R,C), False=old (node_length)

# Read-only globals (broadcast to workers)
GLOBAL_NODE_SEQS = {}
GLOBAL_NODE_AF = {}

# ─────────────────────────────────────────────────────────────────────────────
# Format helpers

def read_dat_global_header(dat_path):
    """Returns (major, minor, block_count). Raises on error."""
    with open(dat_path, "rb") as f:
        magic = f.read(len(DAT_GLOBAL_MAGIC))
        if magic != DAT_GLOBAL_MAGIC:
            raise RuntimeError(f"Invalid .dat magic in {dat_path!r}: {magic!r}")
        major, minor, block_count, _ = GLOBAL_VER_PACK.unpack(f.read(GLOBAL_VER_PACK.size))
        return major, minor, block_count

def make_record_struct_old(node_length: int) -> struct.Struct:
    # Old format: seq[L], bq[L], cigar[L]
    return struct.Struct(f"<h{node_length}s{node_length}s{node_length}shc")

def make_record_struct_new(max_read_len: int, max_cigar_len: int) -> struct.Struct:
    # New format: seq[R], bq[R], cigar[C]
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")

# ─────────────────────────────────────────────────────────────────────────────
# Utility

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

def canonical_variant_key(v_pos, v_type, v_ref, v_alt, node_seq):
    node_len = len(node_seq)

    if v_type == 'I':
        anchor_pos = v_pos
        anchor_base = (v_ref if v_ref and v_ref != "*"
                       else (node_seq[anchor_pos].upper() if 0 <= anchor_pos < node_len else "N"))
        inserted = v_alt or ""
        return f"{anchor_pos}_I_{anchor_base}_{anchor_base + inserted}"

    if v_type == 'D':
        deleted = v_ref or ""
        if v_pos <= 0:
            # no known left anchor — use '*' for ALT and only the deleted bases for REF
            return f"0_D_{deleted}_{v_alt}"
        else:
            anchor_pos = v_pos - 1
            anchor_base = node_seq[anchor_pos].upper() if 0 <= anchor_pos < node_len else "N"
            return f"{anchor_pos}_D_{anchor_base + deleted}_{anchor_base}"

    return f"{v_pos}_{v_type}_{v_ref}_{v_alt}"

def calculate_window_start(variant_pos, window_size):
    return variant_pos - (window_size // 2)

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def load_full_idx_data(idx_path):
    """
    Accepts .idx entry sizes:
      - 30B new: <I Q I I H I I> (nid, offset, block_size, n_records, flags, R, C)
      - 26B old+: <I Q I I H I>   (nid, offset, block_size, n_records, flags, node_length)
      - 22B old : <I Q I I H>
    Returns dict: nid -> (offset, n_records)
    """
    idx_data_map = {}
    print(f"Loading full index data from {idx_path}...")
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                sys.stderr.write(f"Error: Index file {idx_path} too small.\n")
                return None
            (num_nodes_in_idx,) = struct.unpack('<I', f.read(4))
            print(f"  Index entries reported: {num_nodes_in_idx}")
            if num_nodes_in_idx == 0:
                return idx_data_map

            remaining = file_size - 4
            per = remaining // num_nodes_in_idx
            if per == 30:
                entry_size, unpack_fmt = 30, '<I Q I I H I I'
            elif per == 26:
                entry_size, unpack_fmt = 26, '<I Q I I H I'
            elif per == 22:
                entry_size, unpack_fmt = 22, '<I Q I I H'
            else:
                # Heuristic: try 30 first
                entry_size, unpack_fmt = 30, '<I Q I I H I I'

            for i in range(num_nodes_in_idx):
                rec = f.read(entry_size)
                if len(rec) != entry_size:
                    sys.stderr.write(f"Error: truncated idx at entry {i+1}.\n")
                    break
                up = struct.unpack(unpack_fmt, rec)
                if entry_size == 30:
                    node_id_from_idx, offset, _block_size, n_records, _flags, _R, _C = up
                elif entry_size == 26:
                    node_id_from_idx, offset, _block_size, n_records, _flags, _L = up
                else:
                    node_id_from_idx, offset, _block_size, n_records, _flags = up
                idx_data_map[node_id_from_idx] = (offset, n_records)
            print(f"  Loaded {len(idx_data_map)} index entries.")
        return idx_data_map
    except FileNotFoundError:
        sys.stderr.write(f"Error: Index file not found {idx_path}\n")
        return None
    except Exception as e:
        sys.stderr.write(f"Error parsing idx: {e}\n")
        return None

def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*':
        return []
    try:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse CIGAR '{cigar_string}': {e}\n")
        return []

def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_quality_values, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_allele_for_indel=None):
    current_node_pos = read_offset_on_node
    current_read_pos = 0

    saw_insertion_anchor = False
    bq_at_anchor = None

    for length, op in read_cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            if current_node_pos <= target_node_pos < current_node_pos + length:
                offset_in_block = target_node_pos - current_node_pos
                read_idx = current_read_pos + offset_in_block
                if read_idx < len(read_sequence):
                    allele = read_sequence[read_idx].upper()
                    quality = read_quality_values[read_idx] if read_idx < len(read_quality_values) else 0
                    if expected_var_type == 'D':
                        return "REF_STATE_FOR_INDEL", quality
                    if expected_var_type == 'I':
                        saw_insertion_anchor = True
                        bq_at_anchor = quality
                    else:
                        return allele, quality
                else:
                    if expected_var_type == 'I':
                        saw_insertion_anchor = True
                        bq_at_anchor = None
                    else:
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

    if expected_var_type == 'I' and saw_insertion_anchor:
        return "REF_STATE_FOR_INDEL", bq_at_anchor

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
                    # NEW: ignore positions where either side is N
                    if node_base != read_base and node_base != 'N' and read_base != 'N' and op != '=':
                        variants.append((cur_node_p, 'X', read_base, node_base))
                else:
                    break
            node_pos += length
            read_pos += length

        elif op == 'I':
            # NEW: drop long insertions
            if length > MAX_INDEL_LEN:
                read_pos += length
                continue
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            # NEW: ignore insertions containing N
            if inserted_sequence and 'N' not in inserted_sequence:
                ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
                ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < node_seq_len else "*"
                # Also ignore if anchor base is N
                if ref_base_at_anchor != 'N':
                    variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor))
            read_pos += length

        elif op == 'D':
            # NEW: drop long deletions
            if length > MAX_INDEL_LEN:
                node_pos += length
                continue
            if node_pos + length <= node_seq_len:
                deleted_sequence_from_ref = node_sequence[node_pos: node_pos + length].upper()
                # NEW: ignore deletions containing N
                if deleted_sequence_from_ref and 'N' not in deleted_sequence_from_ref:
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
                n_aln = node_pos + i
                win_idx = n_aln - window_start_node
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
                n_aln = node_pos + i
                win_idx = n_aln - window_start_node
                if 0 <= win_idx < tensor_win_size:
                    cigar_ops_indices[win_idx] = op_idx
                    mapqs[win_idx] = mapping_quality
                    bases[win_idx] = BASE_TO_INDEX['*']
                    quals[win_idx] = DEFAULT_QUALITY_PADDING
            node_pos += L
        elif op == 'I':
            # Represent insertion at its left anchor (node_pos - 1)
            anchor_node_pos = node_pos - 1
            anchor_idx = anchor_node_pos - window_start_node
            if 0 <= anchor_idx < tensor_win_size:
                cigar_ops_indices[anchor_idx] = op_idx          # mark CIGAR 'I' at anchor
                mapqs[anchor_idx] = mapping_quality             # place MAPQ at anchor
                bases[anchor_idx] = BASE_TO_INDEX['I']          # use 'I' token (90)

                # average BQ across the inserted run
                if L > 0 and read_pos + L <= qual_len:
                    ins_quals = segment_quality_values[read_pos: read_pos + L]
                    avg_bq = int(round(sum(ins_quals) / L)) if ins_quals else DEFAULT_QUALITY_PADDING
                else:
                    avg_bq = DEFAULT_QUALITY_PADDING
                quals[anchor_idx] = avg_bq

            # Insertions consume read only (no node columns)
            read_pos += L
        elif op == 'S':
            read_pos += L

        if read_pos >= read_seq_len:
            break
    return bases, quals, mapqs, cigar_ops_indices

# ─────────────────────────────────────────────────────────────────────────────
# Worker init & core

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker, need_view_flag, is_new_format_flag):
    global worker_dat_file, worker_base_output_dir, worker_need_view, worker_is_new_format
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
        worker_need_view = bool(need_view_flag)
        worker_is_new_format = bool(is_new_format_flag)
    except FileNotFoundError:
        sys.stderr.write(f"Error [Worker {os.getpid()}]: DAT file not found.\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()}] opening DAT file: {e}\n")
        sys.exit(1)

# ───── NEW: helper for row ordering ALT → OTHER → REF → UNINFORMATIVE ─────
def _row_group_key_for_variant(vt, expected_alt, expected_ref, allele_observed):
    """
    Return a sort key so similar alleles cluster together in tensor rows.
    Ordering: ALT first, then OTHER alleles, then REF, then uninformative.
    """
    if allele_observed is None:
        return (9,)   # uninformative goes last

    if vt == 'X':  # SNV
        if allele_observed == expected_alt:
            return (0,)  # ALT first
        if allele_observed == expected_ref:
            return (2,)  # REF after OTHER
        base_rank = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}.get(allele_observed, 5)
        return (1, base_rank, allele_observed)  # OTHER second

    if vt == 'I':  # Insertion
        if allele_observed == expected_alt:
            return (0, -len(allele_observed), allele_observed)
        if allele_observed == "REF_STATE_FOR_INDEL":
            return (2,)
        if allele_observed == "OTHER_FOR_INDEL":
            return (1,)
        if isinstance(allele_observed, str):
            return (1, -len(allele_observed), allele_observed)
        return (8,)

    if vt == 'D':  # Deletion
        if allele_observed == "*":
            return (0,)
        if allele_observed == "OTHER_FOR_INDEL":
            return (1,)
        if allele_observed == "REF_STATE_FOR_INDEL":
            return (2,)
        return (8,)

    return (9,)

def process_single_node_for_pileup(task_args):
    (node_id, dat_file_offset, n_records,
     min_af_threshold, min_variants_threshold, min_allele_bq_threshold,
     variant_type_to_process, min_depth_threshold) = task_args

    global worker_dat_file, worker_base_output_dir, worker_need_view, worker_is_new_format
    global GLOBAL_NODE_SEQS, GLOBAL_NODE_AF

    tensor_files_generated_for_node = 0
    if worker_dat_file is None or worker_base_output_dir is None:
        return node_id, None, tensor_files_generated_for_node

    node_sequence = GLOBAL_NODE_SEQS.get(node_id, "")
    if not node_sequence:
        return node_id, {}, tensor_files_generated_for_node
    genomead_af_list = GLOBAL_NODE_AF.get(node_id, [])
    node_len = len(node_sequence)
    aligned_read_segments = []

    try:
        # Read per-block header depending on format
        worker_dat_file.seek(dat_file_offset, os.SEEK_SET)
        if worker_is_new_format:
            hdr_bytes = worker_dat_file.read(NEW_BLOCK_HDR_SIZE)
            if len(hdr_bytes) != NEW_BLOCK_HDR_SIZE:
                raise RuntimeError(f"Cannot read new block header at offset {dat_file_offset} for node {node_id}")
            nid2, nrec2, _flags, R, C = NEW_BLOCK_HDR_PACK.unpack(hdr_bytes)
            if nrec2 != n_records or nid2 != node_id:
                # Not fatal
                pass
            rec_struct = make_record_struct_new(int(R or 1), int(C or 1))
            rec_size = rec_struct.size
        else:
            hdr_bytes = worker_dat_file.read(OLD_BLOCK_HDR_SIZE)
            if len(hdr_bytes) != OLD_BLOCK_HDR_SIZE:
                raise RuntimeError(f"Cannot read old block header at offset {dat_file_offset} for node {node_id}")
            nid2, nrec2, _flags, L = OLD_BLOCK_HDR_PACK.unpack(hdr_bytes)
            if nrec2 != n_records or nid2 != node_id:
                pass
            rec_struct = make_record_struct_old(int(L or 1))
            rec_size = rec_struct.size

        records_start = dat_file_offset + (NEW_BLOCK_HDR_SIZE if worker_is_new_format else OLD_BLOCK_HDR_SIZE)
        worker_dat_file.seek(records_start, os.SEEK_SET)

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
                cigar_str_original = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand_char = strand_byte.decode('ascii') if isinstance(strand_byte, (bytes, bytearray)) else chr(strand_byte)
            except UnicodeDecodeError:
                continue
            if not cigar_str_original or len(seq) != len(qual_values):
                continue

            original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
            if not original_decoded_cigar_ops and cigar_str_original != '*':
                continue

            current_read_sequence = seq
            current_quality_values = qual_values
            current_decoded_cigar_ops = original_decoded_cigar_ops
            current_offset_on_node = off_from_file

            if strand_char == '-':
                current_read_sequence = reverse_complement(seq)
                current_quality_values = qual_values[::-1]
                current_decoded_cigar_ops = list(reversed(original_decoded_cigar_ops)) if original_decoded_cigar_ops else []
                # span on node
                alignment_span_on_node = len(current_read_sequence)
                for Lx, opx in original_decoded_cigar_ops:
                    if opx == 'I':
                        alignment_span_on_node -= Lx
                    elif opx == 'D':
                        alignment_span_on_node += Lx
                current_offset_on_node = node_len - alignment_span_on_node - off_from_file
                if current_offset_on_node < 0:
                    continue

            aligned_read_segments.append({
                "offset_on_node": current_offset_on_node,
                "read_sequence": current_read_sequence,
                "processed_quality_values": current_quality_values,
                "cigar_ops": current_decoded_cigar_ops,
                "original_cigar_str": cigar_str_original,
                "strand": strand_char,
                "mapping_quality": mapq_val
            })
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()} for Node {node_id}]: {e}\n")
        return node_id, None, tensor_files_generated_for_node

    if not aligned_read_segments:
        return node_id, {}, tensor_files_generated_for_node

    # NEW: if too many reads, keep only top 400 by "non-M length" (prioritize mismatches/indels)
    if len(aligned_read_segments) > TENSOR_MAX_READ_ROWS:
        def non_m_length(seg):
            # Sum lengths for CIGAR ops that are not 'M'
            return sum(L for L, op in seg["cigar_ops"] if op != 'M')
        aligned_read_segments.sort(key=non_m_length, reverse=True)
        aligned_read_segments = aligned_read_segments[: (TENSOR_MAX_READ_ROWS * 2)]  # 200 * 2 = 400

    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    view_oriented_variant_data = {} if worker_need_view else None
    variant_headers_for_summary = []

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
        vt = v_type
        if (variant_type_to_process == 'snp' and vt != 'X') or \
           (variant_type_to_process == 'indel' and vt not in ('I', 'D')):
            continue

        alt_allele_count = ref_allele_count = other_allele_count = locus_coverage = 0
        alt_allele_base_qualities = []

        expected_ref_for_af = v_ref_from_cigar
        expected_alt_for_af = v_alt_from_cigar
        ref_allele_for_indel_context = None

        if vt == 'D':
            expected_alt_for_af = "*"
            if 0 <= v_pos < len(node_sequence):
                expected_ref_for_af = node_sequence[v_pos]
            ref_allele_for_indel_context = v_ref_from_cigar
        elif vt == 'I':
            expected_ref_for_af = node_sequence[v_pos] if 0 <= v_pos < len(node_sequence) else "*"
            ref_allele_for_indel_context = expected_ref_for_af

        # Collect allele observations and coverage
        for seg in aligned_read_segments:
            allele_observed, bq = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["processed_quality_values"], seg["cigar_ops"],
                v_pos, node_sequence, vt, ref_allele_for_indel_context)

            if allele_observed is not None:
                locus_coverage += 1
                if allele_observed == expected_alt_for_af:
                    alt_allele_count += 1
                    if bq is not None:
                        alt_allele_base_qualities.append(bq)
                elif allele_observed == expected_ref_for_af or (
                        vt in ('I', 'D') and allele_observed == "REF_STATE_FOR_INDEL"):
                    ref_allele_count += 1
                else:
                    other_allele_count += 1

        # NEW: depth gate — skip variants below min_depth
        if locus_coverage < min_depth_threshold:
            continue

        if alt_allele_count < min_variants_threshold:
            continue
        tmp_af = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
        if tmp_af < min_af_threshold:
            continue
        if vt == 'X':
            tmp_bq = sum(alt_allele_base_qualities) / len(alt_allele_base_qualities) if alt_allele_base_qualities else 0.0
            if tmp_bq < min_allele_bq_threshold:
                continue

        current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
        mean_alt_bq = sum(alt_allele_base_qualities) / len(alt_allele_base_qualities) if alt_allele_base_qualities else 0.0

        variant_key_string = canonical_variant_key(v_pos, vt, v_ref_from_cigar, v_alt_from_cigar, node_sequence)

        window_center_pos = v_pos + 1 if vt == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

        if worker_need_view:
            pileup_data_for_view_json = []
            for read_segment_idx, seg_data in enumerate(aligned_read_segments[: TENSOR_MAX_READ_ROWS + 50]):
                row_chars_for_view = get_read_representation_in_window_for_view(
                    seg_data["cigar_ops"], seg_data["offset_on_node"], seg_data["read_sequence"],
                    window_start_pos, TENSOR_WINDOW_SIZE, len(node_sequence))
                if any(char != ' ' for char in row_chars_for_view):
                    bases_for_view = [
                        (PADDING_BASE_INDEX if char == ' ' else BASE_TO_INDEX.get(char.upper(), BASE_TO_INDEX['N']))
                        for char in row_chars_for_view
                    ]
                    pileup_data_for_view_json.append({
                        "bases": bases_for_view,
                        "offset": seg_data["offset_on_node"],
                        "strand": seg_data["strand"],
                        "cigar": seg_data["original_cigar_str"]
                    })
            view_oriented_variant_data[variant_key_string] = {
                "pileup_reads_data": pileup_data_for_view_json[:TENSOR_MAX_READ_ROWS],
                "alt_allele_count": alt_allele_count, "ref_allele_count_at_locus": ref_allele_count,
                "other_allele_count_at_locus": other_allele_count, "coverage_at_locus": locus_coverage,
                "alt_allele_frequency": round(current_alt_freq, 4),
                "mean_alt_allele_base_quality": round(mean_alt_bq, 2)
            }

        # Build channels
        ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
            if 0 <= node_pos_in_window < len(node_sequence):
                ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[node_pos_in_window].upper(), BASE_TO_INDEX['N'])

        genomead_af_row = [0] * TENSOR_WINDOW_SIZE
        if genomead_af_list:
            for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                if 0 <= node_pos_in_window < len(node_sequence):
                    v = genomead_af_list[node_pos_in_window]
                    if isinstance(v, str):
                        if len(v) == 1 and '0' <= v <= '7':
                            b = ord(v) - 48
                        else:
                            try:
                                b = af_float_to_bin(float(v))
                            except Exception:
                                b = 0
                    else:
                        b = af_float_to_bin(v)
                    genomead_af_row[i] = int(b)

        H = 1 + TENSOR_MAX_READ_ROWS
        W = TENSOR_WINDOW_SIZE
        ch1 = np.full((H, W), PADDING_BASE_INDEX, dtype=np.int8)
        ch2 = np.full((H, W), DEFAULT_QUALITY_PADDING, dtype=np.int8)
        ch3 = np.full((H, W), MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
        ch4 = np.full((H, W), DEFAULT_MAPPING_QUALITY_PADDING, dtype=np.int8)
        ch5 = np.full((H, W), CIGAR_PADDING_INDEX, dtype=np.int8)
        ch6 = np.zeros((H, W), dtype=np.int8)

        ch1[0, :] = np.asarray(ref_base_indices_row, dtype=np.int8)
        ch2[0, :] = DEFAULT_QUALITY_PADDING
        ch3[0, :] = MISMATCH_COMPARISON_PADDING_VALUE
        ch4[0, :] = DEFAULT_MAPPING_QUALITY_PADDING
        ch5[0, :] = CIGAR_PADDING_INDEX
        ch6[0, :] = np.asarray(genomead_af_row, dtype=np.int8)

        # ───── NEW: group reads by allele → ALT, OTHER, REF, UNINFORMATIVE ─────
        grouped_reads = []
        for seg_data in aligned_read_segments:
            allele_observed, _bq_tmp = get_allele_from_read_at_node_pos(
                seg_data["offset_on_node"], seg_data["read_sequence"],
                seg_data["processed_quality_values"], seg_data["cigar_ops"],
                v_pos, node_sequence, vt,
                (expected_ref_for_af if vt in ('I','D') else None)
            )
            sort_key = _row_group_key_for_variant(vt, expected_alt_for_af, expected_ref_for_af, allele_observed)
            grouped_reads.append((sort_key, seg_data))
        grouped_reads.sort(key=lambda x: x[0])

        reads_added = 0
        variant_window_index = v_pos - window_start_pos
        for _key, seg_data in grouped_reads:
            if reads_added >= TENSOR_MAX_READ_ROWS:
                break
            mapq = max(0, min(int(seg_data["mapping_quality"]), 127))
            base_idx_row, quality_score_row, mapq_row, cigar_op_row = get_read_tensor_rows_in_window(
                seg_data["cigar_ops"], seg_data["offset_on_node"],
                seg_data["read_sequence"], seg_data["processed_quality_values"],
                mapq,
                window_start_pos, TENSOR_WINDOW_SIZE, len(node_sequence))

            if any(b != PADDING_BASE_INDEX for b in base_idx_row):
                r = 1 + reads_added
                ch1[r, :] = np.asarray(base_idx_row, dtype=np.int8)
                ch2[r, :] = np.asarray(quality_score_row, dtype=np.int8)
                ref_row = ch1[0, :].astype(np.int16)
                read_row = ch1[r, :].astype(np.int16)
                flags = np.full(W, MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
                mask_valid = (read_row != PADDING_BASE_INDEX) & (ref_row != PADDING_BASE_INDEX)
                flags[mask_valid] = (read_row[mask_valid] != ref_row[mask_valid]).astype(np.int8)
                if 0 <= variant_window_index < W and mask_valid[variant_window_index] and flags[variant_window_index] == 1:
                    flags[variant_window_index] = 5  # highlight focal mismatch
                ch3[r, :] = flags
                ch4[r, :] = np.asarray(mapq_row, dtype=np.int8)
                ch5[r, :] = np.asarray(cigar_op_row, dtype=np.int8)
                ch6[r, :] = ch6[0, :]
                reads_added += 1

        tensor = np.stack([ch1, ch2, ch3, ch4, ch5, ch6], axis=0)  # (6, H, W), int8

        node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
        os.makedirs(node_specific_output_dir, exist_ok=True)

        tensor_filename_npy = f"{variant_key_string}.npy"
        tensor_filepath_npy = os.path.join(node_specific_output_dir, tensor_filename_npy)
        try:
            np.save(tensor_filepath_npy, tensor)
            variant_headers_for_summary.append({
                "variant_key": variant_key_string, "tensor_file": tensor_filename_npy,
                "alt_allele_count": alt_allele_count, "ref_allele_count_at_locus": ref_allele_count,
                "other_allele_count_at_locus": other_allele_count, "coverage_at_locus": locus_coverage,
                "alt_allele_frequency": round(current_alt_freq, 4),
                "mean_alt_allele_base_quality": round(mean_alt_bq, 2)
            })
            tensor_files_generated_for_node += 1
        except Exception as e:
            sys.stderr.write(f"Error saving tensor for {variant_key_string}: {e}\n")

    if variant_headers_for_summary:
        summary_path = os.path.join(worker_base_output_dir, str(node_id), "variant_summary.json")
        with open(summary_path, 'w') as f:
            json.dump({"node_id": node_id, "node_length": len(node_sequence),
                       "variants_passing_af_filter": variant_headers_for_summary}, f, indent=2)

    return node_id, (view_oriented_variant_data or {}), tensor_files_generated_for_node

# ─────────────────────────────────────────────────────────────────────────────
# Pileup viewing (unchanged)

def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if max_variants_to_display == 0:
        return
    if not node_data_for_display_view:
        print(f"Info: No pileup data for node {node_id_str_for_display} (no variants met all filters).")
        return

    print(f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))

    for i, variant_key in enumerate(sorted_variant_keys):
        if i >= max_variants_to_display:
            print(f"\n  ... ({len(sorted_variant_keys) - i} more variants not shown due to --view limit)")
            break

        variant_data = node_data_for_display_view[variant_key]
        v_pos, v_type = int(variant_key.split('_')[0]), variant_key.split('_')[1]

        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

        print(f"\n--- Variant: {variant_key} ---")

        ref_chars = []
        for j in range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE):
            if 0 <= j < len(full_node_sequence):
                ref_chars.append(full_node_sequence[j])
            else:
                ref_chars.append('0')
        print(f"  Node Ref: {''.join(ref_chars)}")

        marker_pos_in_window = v_pos - window_start_pos
        marker_line = [' '] * TENSOR_WINDOW_SIZE
        if 0 <= marker_pos_in_window < TENSOR_WINDOW_SIZE:
            marker_line[marker_pos_in_window] = '^'
        print(f"  Marker  : {''.join(marker_line)}")

        for j, read_entry in enumerate(variant_data.get("pileup_reads_data", [])):
            if j >= max_reads_to_display_per_variant:
                print(f"  ... ({len(variant_data.get('pileup_reads_data', [])) - j} more reads not shown)")
                break
            bases_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in read_entry["bases"]])
            print(f"  Read {j + 1:3d}: {bases_str} (CIGAR:{read_entry['cigar']}, STRAND:{read_entry.get('strand', '?')})")

        print(
            f"  Alt Count: {variant_data.get('alt_allele_count', 'N/A')}, Ref Count: {variant_data.get('ref_allele_count_at_locus', 'N/A')}, Coverage: {variant_data.get('coverage_at_locus', 'N/A')}")
        print(
            f"  Alt Freq: {variant_data.get('alt_allele_frequency', 0.0):.4f}, Mean Alt BQ: {variant_data.get('mean_alt_allele_base_quality', 0.0):.2f}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered tensors from alignment data (supports old & new .dat/.idx formats).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat alignment file")
    parser.add_argument("idx", help=".idx index file")
    parser.add_argument("output", help="Base output directory")
    parser.add_argument("candidate_variants_json",
                        help="JSON with nodes and sequences to process.")
    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Number of worker processes")
    parser.add_argument("--chunksize", type=int, default=128, help="Task chunksize for executor.map()")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                        help="Print pileups for top N variants per node (-1 for all)")
    parser.add_argument("--max_view_reads", type=int, default=20, help="Max reads to show per pileup in view mode")
    parser.add_argument("--min_af", type=float, default=0.1, help="Minimum allele frequency")
    parser.add_argument("--min_variants", type=int, default=3, help="Minimum ALT-supporting reads")
    parser.add_argument("--min_allele_bq", type=float, default=10.0, help="Minimum mean BQ for ALT")
    parser.add_argument("--variant_type", type=str, default='all', choices=['snp', 'indel', 'all'],
                        help="Which variants to output tensors for")
    parser.add_argument("--milestone", type=int, default=100000, help="The milestone to print logs")
    parser.add_argument("--min_depth", type=int, default=20,
                        help="Minimum locus coverage required to emit a tensor")
    args = parser.parse_args()

    if not all([os.path.isfile(args.dat), os.path.isfile(args.idx), os.path.isfile(args.candidate_variants_json)]):
        sys.exit("Error: One or more input files (dat, idx, or json) were not found.")
    os.makedirs(args.output, exist_ok=True)

    need_view = (args.view is not None and args.view != 0)

    # Detect .dat format once (share to workers)
    try:
        major, minor, block_count = read_dat_global_header(args.dat)
        # New writer used minor >= 5 when switching to per-block (R,C)
        is_new_format = (minor >= 5)
        print(f".dat header -> version {major}.{minor}, blocks={block_count} → "
              f"{'NEW' if is_new_format else 'OLD'} format")
    except Exception as e:
        sys.exit(f"Error reading .dat global header: {e}")

    # Load candidate nodes (sequence + AF)
    node_sequences = {}
    node_af_data = {}
    node_ids_to_process = set()
    print(f"Loading nodes from {args.candidate_variants_json}...")
    try:
        with open(args.candidate_variants_json, 'r') as f:
            data = json.load(f)
            for node_obj in data.get('nodes', []):
                node_id_str = node_obj.get('node_id')
                sequence = node_obj.get('sequence')
                af_list = node_obj.get('genomead_af', [])
                if node_id_str and sequence:
                    try:
                        node_id_int = int(node_id_str)
                        node_sequences[node_id_int] = sequence.upper()
                        node_af_data[node_id_int] = af_list
                        node_ids_to_process.add(node_id_int)
                    except ValueError:
                        pass
    except Exception as e:
        sys.exit(f"Error reading or parsing JSON file: {e}")

    print(f"Found {len(node_sequences)} nodes with integer-compatible IDs to process.")

    full_idx_data = load_full_idx_data(args.idx)
    if not full_idx_data:
        sys.exit("Failed to load index data.")

    tasks = []
    missing = 0
    for node_id in node_ids_to_process:
        if node_id in full_idx_data:
            offset, n_records = full_idx_data[node_id]
            tasks.append((node_id, offset, n_records, args.min_af, args.min_variants,
                          args.min_allele_bq, args.variant_type, args.min_depth))
        else:
            missing += 1
    if missing:
        print(f"Warning: {missing} nodes from JSON not found in idx; skipped.")
    if not tasks:
        sys.exit("No valid tasks to run after processing JSON and idx.")

    tasks.sort(key=lambda t: t[1])  # stream-friendly

    global GLOBAL_NODE_SEQS, GLOBAL_NODE_AF
    GLOBAL_NODE_SEQS = node_sequences
    GLOBAL_NODE_AF = node_af_data

    total_tasks = len(tasks)
    print(f"\nSubmitting {total_tasks} tasks to {args.num_workers} workers...")

    total_tensors = 0
    processed = 0
    nodes_since_last_report = 0
    tensors_since_last_report = 0
    batch_start_time = time.time()

    with ProcessPoolExecutor(max_workers=args.num_workers,
                             initializer=init_worker,
                             initargs=(args.dat, args.output, need_view, is_new_format)) as executor:
        for node_id, view_data, tensor_count in executor.map(
                process_single_node_for_pileup, tasks, chunksize=max(1, args.chunksize)):
            processed += 1
            nodes_since_last_report += 1
            total_tensors += tensor_count
            tensors_since_last_report += tensor_count

            if need_view and view_data:
                node_sequence_for_view = GLOBAL_NODE_SEQS.get(node_id, "")
                display_pileup_data(view_data, str(node_id), node_sequence_for_view,
                                    args.max_view_reads,
                                    args.view if args.view != -1 else float('inf'))

            if nodes_since_last_report >= args.milestone or processed == total_tasks:
                elapsed_time = time.time() - batch_start_time
                rate = nodes_since_last_report / elapsed_time if elapsed_time > 0 else 0.0
                print(
                    f"  Processed {nodes_since_last_report} nodes (total: {processed}/{total_tasks}) in {elapsed_time:.2f}s "
                    f"({rate:.1f} nodes/sec). Batch tensors: {tensors_since_last_report}. Total tensors: {total_tensors}.")
                nodes_since_last_report = 0
                tensors_since_last_report = 0
                batch_start_time = time.time()

    print(f"\nProcessing complete. Total tensors generated: {total_tensors}.")

if __name__ == '__main__':
    main()

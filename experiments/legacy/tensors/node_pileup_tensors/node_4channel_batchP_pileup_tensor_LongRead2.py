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
# Constants
# NOTE: The fixed RECORD_STRUCT/RECORD_SIZE from the old format are no longer used
# for reading; the new format uses variable record sizes per node_length. We keep
# them defined to avoid touching other parts of the file, but reading now uses
# make_record_struct(node_len) inside process_single_node_for_pileup.
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size

# New-format: per-block header and variable-size record builder
# Block header layout: <I I H I>  -> node_id (u32), n_records (u32), flags (u16), node_length (u32)
BLOCK_HDR_PACK = struct.Struct("<I I H I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size  # 14 bytes

def make_record_struct(node_length: int) -> struct.Struct:
    """
    New-format per-node record struct:
      <h {L}s {L}s {L}s h c>
        - i16 offset
        - seq[L] bytes (null-padded)
        - bq[L] bytes (null-padded)
        - cigar[L] bytes (ASCII, null-padded)  # CHANGED: was 30, now node_length
        - i16 rq (MAPQ)
        - char strand ('+' / '-')
    """
    return struct.Struct(f"<h{node_length}s{node_length}s{node_length}shc")

BASE_TO_INDEX = {
    'A': 20, 'C': 30, 'G': 50, 'T': 70,
    'N': 10,
    '*': 90,
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
    10: 'N',
    90: '*',
    0: ' '
}

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
DEFAULT_QUALITY_PADDING = 0
DEFAULT_MAPPING_QUALITY_PADDING = -1
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1

# Globals for worker process state
worker_dat_file = None
worker_base_output_dir = None
worker_need_view = False  # gate view building

# Read-only globals (filled before pool creation; inherited by workers via fork COW)
GLOBAL_NODE_SEQS = {}
GLOBAL_NODE_AF = {}

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
# ─────────────────────────────────────────────────────────────────────────────

def af_float_to_bin(x: float) -> int:
    """Bin AF into 0..7:
       0–1e-6→0, 1e-6–1e-5→1, 1e-5–1e-4→2, 1e-4–1e-3→3,
       1e-3–1e-2→4, 1e-2–0.1→5, 0.1–0.5→6, 0.5–1.0→7"""
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

# NEW: filename canonicalization (only naming; logic unchanged)
def canonical_variant_key(v_pos, v_type, v_ref, v_alt, node_seq):
    """
    Convert current internal variant tuple into a left-anchored, 0-based filename key.

    Current internal representation in this script:
      - SNP/X: (v_pos = mismatch pos), v_ref = REF base, v_alt = ALT base
      - INS/I: v_pos = anchor pos (already left of insertion),
               v_ref = anchor base, v_alt = inserted sequence
      - DEL/D: v_pos = first deleted base (NOT the anchor),
               v_ref = deleted sequence, v_alt = "*"

    We output:
      - X : f"{v_pos}_X_{v_ref}_{v_alt}"
      - I : f"{anchor}_I_{anchorBase}_{anchorBase+inserted}"
      - D : f"{anchor}_D_{anchorBase+deleted}_{anchorBase}"
    """
    node_len = len(node_seq)
    if v_type == 'I':
        anchor_pos = v_pos
        anchor_base = (v_ref if v_ref and v_ref != "*" else
                       (node_seq[anchor_pos].upper() if 0 <= anchor_pos < node_len else "N"))
        inserted = v_alt or ""
        return f"{anchor_pos}_{v_type}_{anchor_base}_{anchor_base + inserted}"

    if v_type == 'D':
        # v_pos points to first deleted base; anchor is just before it
        anchor_pos = max(0, v_pos - 1)
        anchor_base = node_seq[anchor_pos].upper() if 0 <= anchor_pos < node_len else "N"
        deleted = v_ref or ""
        return f"{anchor_pos}_{v_type}_{anchor_base + deleted}_{anchor_base}"

    # SNP / mismatch unchanged
    return f"{v_pos}_{v_type}_{v_ref}_{v_alt}"


def calculate_window_start(variant_pos, window_size):
    center_index = window_size // 2
    return variant_pos - center_index

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def load_full_idx_data(idx_path):
    """
    New .idx entries are 26 bytes: <I Q I I H I>
      node_id (u32), offset (u64), block_size (u32), n_records (u32), flags (u16), node_length (u32)
    For backward compatibility, we also accept 22B entries (<I Q I I H>), but we still
    return only (offset, n_records) because the rest of the code expects that.
    """
    idx_data_map = {}
    print(f"Loading full index data from {idx_path}...")
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                sys.stderr.write(f"Error: Index file {idx_path} is too small.\n")
                return None
            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                sys.stderr.write(f"Error: Could not read number of nodes from {idx_path}.\n")
                return None
            num_nodes_in_idx = struct.unpack('<I', num_nodes_bytes)[0]
            print(f"  Index file reports {num_nodes_in_idx} total node entries. Reading all entries...")
            if num_nodes_in_idx == 0:
                print("Successfully loaded 0 distinct node entries.")
                return idx_data_map

            processed_entries = 0

            # Determine entry size by remaining length
            remaining = file_size - 4
            # Prefer new format (26B). Fall back to 22B if it divides exactly.
            if remaining // num_nodes_in_idx == 26:
                entry_size = 26
                unpack_fmt = '<I Q I I H I'
            elif remaining // num_nodes_in_idx == 22:
                entry_size = 22
                unpack_fmt = '<I Q I I H'
            else:
                # Try to read as new format by default
                entry_size = 26
                unpack_fmt = '<I Q I I H I'

            for i in range(num_nodes_in_idx):
                record_bytes = f.read(entry_size)
                if len(record_bytes) < entry_size:
                    sys.stderr.write(f"Error: Index file ended prematurely at record {i + 1}.\n")
                    break
                unpacked = struct.unpack(unpack_fmt, record_bytes)
                if entry_size == 26:
                    node_id_from_idx, offset, _block_size, n_records, _flags, _node_len = unpacked
                else:
                    node_id_from_idx, offset, _block_size, n_records, _flags = unpacked
                idx_data_map[node_id_from_idx] = (offset, n_records)
                processed_entries += 1
                if processed_entries % 5_000_000 == 0:
                    print(f"    Loaded {processed_entries}/{num_nodes_in_idx} index entries...")
            print(f"Successfully loaded {len(idx_data_map)} distinct node entries.")
        return idx_data_map
    except FileNotFoundError:
        sys.stderr.write(f"Error: Index file not found at {idx_path}\n")
        return None
    except Exception as e:
        sys.stderr.write(f"Error parsing full index file {idx_path}: {e}\n")
        return None

def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*':
        return []
    try:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse CIGAR string '{cigar_string}': {e}\n")
        return []

def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_quality_values, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_allele_for_indel=None):
    current_node_pos = read_offset_on_node
    current_read_pos = 0

    # NEW: track if we've seen the insertion anchor in a match block;
    # if we don't later see an 'I' at that anchor, we will return REF.
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
                        # Deletions: seeing the base at target means REF (no deletion spanning it)
                        return "REF_STATE_FOR_INDEL", quality

                    if expected_var_type == 'I':
                        # Insertions: do NOT return yet; mark that the anchor is covered,
                        # then continue scanning to see if the next op is an 'I' at this anchor.
                        saw_insertion_anchor = True
                        bq_at_anchor = quality
                        # fall through (no return) to keep scanning
                    else:
                        return allele, quality
                else:
                    # out of read range
                    if expected_var_type == 'I':
                        saw_insertion_anchor = True
                        bq_at_anchor = None
                    else:
                        return None, None
            # advance through the match
            current_node_pos += length
            current_read_pos += length

        elif op == 'I':
            # For an insertion, the anchor is the base just before current_node_pos.
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                if current_read_pos + length <= len(read_sequence):
                    qualities = read_quality_values[current_read_pos: current_read_pos + length]
                    mean_quality = sum(qualities) / len(qualities) if qualities else 0.0
                    return read_sequence[current_read_pos: current_read_pos + length].upper(), mean_quality
                return None, None
            current_read_pos += length  # ref doesn't advance on I

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
            current_node_pos += length  # ref advances on D

        elif op == 'S':
            current_read_pos += length
        elif op == 'N':
            current_node_pos += length

        # original early-stop logic kept as-is
        if current_node_pos > target_node_pos + 1 and not (
                expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):
            break

    # If we saw the insertion anchor in a match but never found an 'I' at that anchor,
    # count this read as REF for the insertion.
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
                if 0 <= win_idx < window_size:
                    if r_aln < read_seq_len:
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
        elif op in ('I', 'S'):
            read_pos += L

        if read_pos >= read_seq_len:
            break
    return bases, quals, mapqs, cigar_ops_indices

# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function
# ─────────────────────────────────────────────────────────────────────────────

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker, need_view_flag):
    global worker_dat_file, worker_base_output_dir, worker_need_view
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
        worker_need_view = bool(need_view_flag)
    except FileNotFoundError:
        sys.stderr.write(f"Error [Worker {os.getpid()}]: DAT file not found.\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()}] opening DAT file: {e}\n")
        sys.exit(1)

def process_single_node_for_pileup(task_args):
    (node_id, dat_file_offset, n_records,
     min_af_threshold, min_variants_threshold, min_allele_bq_threshold,
     variant_type_to_process) = task_args

    global worker_dat_file, worker_base_output_dir, worker_need_view
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
        # ── NEW FORMAT READ ───────────────────────────────────────────────────
        # Read the per-block header at dat_file_offset to recover node_length
        worker_dat_file.seek(dat_file_offset, os.SEEK_SET)
        block_hdr = worker_dat_file.read(BLOCK_HDR_SIZE)
        if len(block_hdr) != BLOCK_HDR_SIZE:
            raise RuntimeError(f"Cannot read block header at offset {dat_file_offset} for node {node_id}")
        nid2, nrec2, _flags, node_length_from_dat = BLOCK_HDR_PACK.unpack(block_hdr)
        if nrec2 != n_records or nid2 != node_id:
            # Not fatal; continue with what idx gave us
            pass

        # Build the per-node record struct and read all records
        rec_struct = make_record_struct(int(node_length_from_dat))
        rec_size = rec_struct.size

        # Seek to start of records for this block
        records_start = dat_file_offset + BLOCK_HDR_SIZE
        worker_dat_file.seek(records_start, os.SEEK_SET)

        bulk = worker_dat_file.read(n_records * rec_size)
        if len(bulk) < n_records * rec_size:
            n_records = len(bulk) // rec_size
            bulk = bulk[: n_records * rec_size]

        for (off_from_file, raw_seq, raw_qual, raw_cigar, mapq_val, strand_byte) in rec_struct.iter_unpack(bulk):
            if mapq_val < 10:
                continue

            try:
                # raw_seq/raw_qual are already node_length long, null-padded
                seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                qual_values = list(raw_qual.rstrip(b'\0'))
                cigar_str_original = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                # strand_byte is a single byte; keep ascii char
                strand_char = strand_byte.decode('ascii') if isinstance(strand_byte, (bytes, bytearray)) else chr(strand_byte)
            except UnicodeDecodeError:
                continue

            if not seq or len(seq) != len(qual_values):
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
                current_decoded_cigar_ops = [op for op in reversed(original_decoded_cigar_ops)] if original_decoded_cigar_ops else []
                # alignment_span_on_node = len(current_read_sequence)
                alignment_span_on_node = len(current_read_sequence)
                # Adjust span: +I, -D (reference-consuming vs non-consuming edits)
                for L, op in original_decoded_cigar_ops:
                    if op == 'I':
                        alignment_span_on_node -= L
                    elif op == 'D':
                        alignment_span_on_node += L
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
    if len(aligned_read_segments) > 100000:
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
        if variant_type_to_process == 'indel' and (v_type not in ('I', 'D')):
            continue

        alt_allele_count = ref_allele_count = other_allele_count = locus_coverage = 0
        alt_allele_base_qualities = []

        expected_ref_for_af = v_ref_from_cigar
        expected_alt_for_af = v_alt_from_cigar
        ref_allele_for_indel_context = None

        if v_type == 'D':
            expected_alt_for_af = "*"
            if 0 <= v_pos < node_len:
                expected_ref_for_af = node_sequence[v_pos]
            ref_allele_for_indel_context = v_ref_from_cigar
        elif v_type == 'I':
            expected_ref_for_af = node_sequence[v_pos] if 0 <= v_pos < node_len else "*"
            ref_allele_for_indel_context = expected_ref_for_af

        for seg in aligned_read_segments:
            allele_observed, bq = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["processed_quality_values"], seg["cigar_ops"],
                v_pos, node_sequence, v_type, ref_allele_for_indel_context)

            if allele_observed is not None:
                locus_coverage += 1
                if allele_observed == expected_alt_for_af:
                    alt_allele_count += 1
                    if bq is not None:
                        alt_allele_base_qualities.append(bq)
                elif allele_observed == expected_ref_for_af or (
                        v_type in ('I', 'D') and allele_observed == "REF_STATE_FOR_INDEL"):
                    ref_allele_count += 1
                else:
                    other_allele_count += 1

        if alt_allele_count < min_variants_threshold:
            continue

        if v_type == 'X':
            current_alt_freq_tmp = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
            if current_alt_freq_tmp < min_af_threshold:
                continue
            mean_alt_bq_tmp = sum(alt_allele_base_qualities) / len(alt_allele_base_qualities) if alt_allele_base_qualities else 0.0
            if mean_alt_bq_tmp < min_allele_bq_threshold:
                continue

        current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
        mean_alt_bq = sum(alt_allele_base_qualities) / len(alt_allele_base_qualities) if alt_allele_base_qualities else 0.0

        # variant_key_string = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
        variant_key_string = canonical_variant_key(
            v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar, node_sequence
        )

        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

        # View data (gated)
        if worker_need_view:
            pileup_data_for_view_json = []
            for read_segment_idx, seg_data in enumerate(aligned_read_segments[: TENSOR_MAX_READ_ROWS + 50]):
                row_chars_for_view = get_read_representation_in_window_for_view(
                    seg_data["cigar_ops"], seg_data["offset_on_node"], seg_data["read_sequence"],
                    window_start_pos, TENSOR_WINDOW_SIZE, node_len)
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
            if 0 <= node_pos_in_window < node_len:
                ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[node_pos_in_window].upper(), BASE_TO_INDEX['N'])

        genomead_af_row = [0] * TENSOR_WINDOW_SIZE
        if genomead_af_list:
            for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
                if 0 <= node_pos_in_window < node_len:
                    v = genomead_af_list[node_pos_in_window]
                    # Accept either digits-per-base string ('0'..'7') or float AFs
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

        # Pre-allocate numpy arrays for speed (drop Torch)
        H = 1 + TENSOR_MAX_READ_ROWS
        W = TENSOR_WINDOW_SIZE
        ch1 = np.full((H, W), PADDING_BASE_INDEX, dtype=np.int8)
        ch2 = np.full((H, W), DEFAULT_QUALITY_PADDING, dtype=np.int8)
        ch3 = np.full((H, W), MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
        ch4 = np.full((H, W), DEFAULT_MAPPING_QUALITY_PADDING, dtype=np.int8)
        ch5 = np.full((H, W), CIGAR_PADDING_INDEX, dtype=np.int8)
        ch6 = np.zeros((H, W), dtype=np.int8)

        # Row 0: reference/context rows
        ch1[0, :] = np.asarray(ref_base_indices_row, dtype=np.int8)
        ch2[0, :] = DEFAULT_QUALITY_PADDING
        ch3[0, :] = MISMATCH_CHANNEL_REF_ROW_VALUE
        ch4[0, :] = DEFAULT_MAPPING_QUALITY_PADDING
        ch5[0, :] = CIGAR_PADDING_INDEX
        ch6[0, :] = np.asarray(genomead_af_row, dtype=np.int8)

        # Reads
        reads_added = 0
        variant_window_index = v_pos - window_start_pos
        for seg_data in aligned_read_segments:
            if reads_added >= TENSOR_MAX_READ_ROWS:
                break
            mapq = max(0, min(int(seg_data["mapping_quality"]), 127))
            base_idx_row, quality_score_row, mapq_row, cigar_op_row = get_read_tensor_rows_in_window(
                seg_data["cigar_ops"], seg_data["offset_on_node"],
                seg_data["read_sequence"], seg_data["processed_quality_values"],
                mapq,
                window_start_pos, TENSOR_WINDOW_SIZE, node_len)

            if any(b != PADDING_BASE_INDEX for b in base_idx_row):
                r = 1 + reads_added
                ch1[r, :] = np.asarray(base_idx_row, dtype=np.int8)
                ch2[r, :] = np.asarray(quality_score_row, dtype=np.int8)
                # mismatch flags
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
                ch6[r, :] = ch6[0, :]  # reuse AF row
                reads_added += 1

        # Save .npy directly (NumPy only)
        tensor = np.stack([ch1, ch2, ch3, ch4, ch5, ch6], axis=0)  # (6, H, W), int8

        # Lazily create directory only if we’re actually saving
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
            sys.stderr.write(f"Error creating/saving tensor for {variant_key_string}: {e}\n")

    if variant_headers_for_summary:
        summary_path = os.path.join(worker_base_output_dir, str(node_id), "variant_summary.json")
        with open(summary_path, 'w') as f:
            json.dump({"node_id": node_id, "node_length": node_len,
                       "variants_passing_af_filter": variant_headers_for_summary}, f, indent=2)

    return node_id, (view_oriented_variant_data or {}), tensor_files_generated_for_node

# ─────────────────────────────────────────────────────────────────────────────
# Pileup Viewing Function
# ─────────────────────────────────────────────────────────────────────────────
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
# Main function
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered tensors from alignment data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat alignment file")
    parser.add_argument("idx", help=".idx index file")
    parser.add_argument("output", help="Base output directory")
    parser.add_argument("candidate_variants_json",
                        help="JSON file containing nodes and their sequences to process.")
    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Number of worker processes")
    parser.add_argument("--chunksize", type=int, default=256, help="Task chunksize for executor.map()")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                        help="Print pileups for top N variants per node (-1 for all)")
    parser.add_argument("--max_view_reads", type=int, default=20, help="Max reads to show per pileup in view mode")
    parser.add_argument("--min_af", type=float, default=0.1, help="Minimum allele frequency to process a variant")
    parser.add_argument("--min_variants", type=int, default=3, help="Alternate allele count must be >= this value")
    parser.add_argument("--min_allele_bq", type=float, default=10.0,
                        help="Minimum mean base quality of allele-supporting bases")
    parser.add_argument("--variant_type", type=str, default='all', choices=['snp', 'indel', 'all'],
                        help="Type of variants to output tensors for: 'snp', 'indel', or 'all'.")
    args = parser.parse_args()

    if not all([os.path.isfile(args.dat), os.path.isfile(args.idx), os.path.isfile(args.candidate_variants_json)]):
        sys.exit("Error: One or more input files (dat, idx, or json) were not found.")
    os.makedirs(args.output, exist_ok=True)

    # Gate view-building in workers
    need_view = (args.view is not None and args.view != 0)

    # Load candidate nodes (sequence + AF) into read-only globals
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

    # Prepare tiny tasks: do NOT include sequences/AF (workers read from globals)
    tasks = []
    missing = 0
    for node_id in node_ids_to_process:
        if node_id in full_idx_data:
            offset, n_records = full_idx_data[node_id]
            tasks.append((node_id, offset, n_records, args.min_af, args.min_variants,
                          args.min_allele_bq, args.variant_type))
        else:
            missing += 1
    if missing:
        print(f"Warning: {missing} node IDs from JSON were not found in the index file and will be skipped.")
    if not tasks:
        sys.exit("No valid tasks to run after processing JSON and index file.")

    # Sort by DAT offset → streaming I/O instead of random seeks
    tasks.sort(key=lambda t: t[1])

    # Publish read-only globals for workers (copy-on-write under fork)
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

    # Use executor.map to avoid flooding with 300k futures; tune chunksize
    with ProcessPoolExecutor(max_workers=args.num_workers,
                             initializer=init_worker,
                             initargs=(args.dat, args.output, need_view)) as executor:
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

            if nodes_since_last_report >= 100000 or processed == total_tasks:
                elapsed_time = time.time() - batch_start_time
                rate = nodes_since_last_report / elapsed_time if elapsed_time > 0 else 0.0
                print(
                    f"  Processed batch of {nodes_since_last_report} nodes (total: {processed}/{total_tasks}) in {elapsed_time:.2f}s "
                    f"({rate:.1f} nodes/sec). Tensors in batch: {tensors_since_last_report}. Total tensors: {total_tensors}.")
                nodes_since_last_report = 0
                tensors_since_last_report = 0
                batch_start_time = time.time()

    print(f"\nProcessing complete. Total tensors generated: {total_tensors}.")

if __name__ == '__main__':
    main()

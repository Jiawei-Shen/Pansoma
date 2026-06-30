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
import torch  # Added for tensor operations

# ─────────────────────────────────────────────────────────────────────────────
# Constants
RECORD_STRUCT = struct.Struct("<h150s150s20shc")  # Read offset, sequence, RAW QUALITIES, CIGAR, MAPQ, strand
RECORD_SIZE = RECORD_STRUCT.size

BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 5, ' ': 6, '-': 6}
# MODIFIED: INDEX_TO_BASE_FOR_VIEW now maps padding index 6 to '0' for display
INDEX_TO_BASE_FOR_VIEW = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N', 5: '*', 6: '0'}

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
PADDING_BASE_INDEX = BASE_TO_INDEX[' ']  # This is index 6, will be displayed as '0' due to above
DEFAULT_QUALITY_PADDING = 0
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = 0

# Globals for worker process state
worker_dat_file = None
worker_base_output_dir = None


# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
# (reverse_complement, parse_idx_file_for_single_node, load_node_sequence_from_gfa,
#  decode_cigar_to_int_ops, get_allele_from_read_at_node_pos, detect_variants_from_cigar,
#  get_read_representation_in_window_for_view, get_read_tensor_rows_in_window remain the same)
# ... [Unchanged helper functions: reverse_complement, parse_idx_file_for_single_node,
# load_node_sequence_from_gfa, decode_cigar_to_int_ops, get_allele_from_read_at_node_pos,
# detect_variants_from_cigar, get_read_representation_in_window_for_view,
# get_read_tensor_rows_in_window] ...
# (These functions are robust to window_start being outside node boundaries due to
# initialization with padding and conditional filling based on actual_node_pos)
# ─────────────────────────────────────────────────────────────────────────────

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]


def parse_idx_file_for_single_node(idx_path, target_node_id):
    node_info = None
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                print(f"❌ Error: Index file {idx_path} is too small.", file=sys.stderr)
                return None
            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                print(f"❌ Error: Could not read number of nodes from {idx_path}.", file=sys.stderr)
                return None
            num_nodes = struct.unpack('<I', num_nodes_bytes)[0]
            print(f"🔹 Index file contains {num_nodes} nodes. Searching for node {target_node_id}...")
            found = False
            for i in range(num_nodes):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                    print(f"❌ Error: Index file ended prematurely while reading record {i + 1}/{num_nodes}.",
                          file=sys.stderr)
                    break
                node_id, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                if node_id == target_node_id:
                    node_info = (offset, n_records)
                    found = True
                    print(f"✔ Found node {target_node_id} in index: Offset={offset}, N_Records={n_records}")
                    break
            if not found:
                print(f"❌ Error: Node ID {target_node_id} not found in the index file {idx_path}.", file=sys.stderr)
                return None
        return node_info
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
        sys.exit(1)


def load_node_sequence_from_gfa(gfa_path, target_node_id):
    node_sequence = None
    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file to find sequence for node {target_node_id}: {gfa_path}")
            line_counter = 0
            for line in f:
                line_counter += 1
                if line_counter % 5_000_000 == 0:
                    print(f"  Checked {line_counter:,} lines in GFA file...")
                if not line.startswith('S\t'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                try:
                    nid = int(parts[1])
                except ValueError:
                    continue
                if nid == target_node_id:
                    node_sequence = parts[2]
                    print(f"✔ Found sequence for node {target_node_id} in GFA.")
                    break
            if line_counter > 0: print(f"✔ Finished GFA scan after {line_counter:,} lines.")
            if node_sequence is None:
                print(f"❌ Error: Sequence for node ID {target_node_id} not found in GFA file {gfa_path}.",
                      file=sys.stderr)
    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
        sys.exit(1)
    return node_sequence


def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*':
        return []
    ops = []
    try:
        for length_str, op_char in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string):
            ops.append((int(length_str), op_char))
        return ops
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []


def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_allele_for_indel=None):
    current_node_pos = read_offset_on_node
    current_read_pos = 0
    for length, op in read_cigar_ops_decoded:
        if op == 'M' or op == '=' or op == 'X':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I': return "REF_STATE_FOR_INDEL"
                if expected_var_type == 'D': return "REF_STATE_FOR_INDEL"
                offset_in_block = target_node_pos - current_node_pos
                if current_read_pos + offset_in_block < len(read_sequence):
                    return read_sequence[current_read_pos + offset_in_block].upper()
                return None
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                return read_sequence[current_read_pos: current_read_pos + length].upper()
            current_read_pos += length
        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I': return "OTHER_FOR_INDEL"
                if expected_var_type == 'D':
                    deleted_seq_in_read = node_sequence[current_node_pos: current_node_pos + length]
                    if deleted_seq_in_read == expected_ref_allele_for_indel:
                        return "*"
                    else:
                        return "OTHER_FOR_INDEL"
                return "*"
            current_node_pos += length
        elif op == 'S':
            current_read_pos += length
        elif op == 'N':
            current_node_pos += length
        if current_node_pos > target_node_pos + 1 and op in ('M', '=', 'X', 'D', 'N'):
            if not (expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):
                break
    return None


def detect_variants_from_cigar(offset_on_node, cigar_ops_decoded, read_sequence, node_sequence):
    variants = []
    node_pos = offset_on_node
    read_pos = 0
    for length, op in cigar_ops_decoded:
        if op == 'M' or op == '=' or op == 'X':
            for i in range(length):
                current_node_pos = node_pos + i
                current_read_pos = read_pos + i
                if current_node_pos < len(node_sequence) and current_read_pos < len(read_sequence):
                    node_base = node_sequence[current_node_pos].upper()
                    read_base = read_sequence[current_read_pos].upper()
                    if node_base != read_base and op != '=':
                        variants.append((current_node_pos, 'X', read_base, node_base))
            node_pos += length
            read_pos += length
        elif op == 'I':
            inserted_sequence = read_sequence[read_pos: read_pos + length].upper()
            ref_anchor_pos = node_pos - 1 if node_pos > 0 else 0
            ref_base_at_anchor = node_sequence[ref_anchor_pos].upper() if 0 <= ref_anchor_pos < len(
                node_sequence) else ""
            variants.append((ref_anchor_pos, 'I', inserted_sequence, ref_base_at_anchor if ref_base_at_anchor else "*"))
            read_pos += length
        elif op == 'D':
            deleted_sequence_from_ref = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= len(
                node_sequence) else ""
            if deleted_sequence_from_ref:
                variants.append((node_pos, 'D', '*', deleted_sequence_from_ref))
            node_pos += length
        elif op == 'S':
            read_pos += length
        elif op == 'N':
            node_pos += length
    return variants


def get_read_representation_in_window_for_view(segment_cigar_ops, segment_offset_on_node, segment_read_sequence,
                                               window_start_node, window_size, node_len):
    # Initializes with ' ' which maps to PADDING_BASE_INDEX (6).
    # INDEX_TO_BASE_FOR_VIEW[6] is now '0', so this padding will appear as '0'.
    window_char_representation = [' '] * window_size
    current_node_pos_in_read = segment_offset_on_node
    current_read_pos_in_read = 0
    for cigar_len, cigar_op in segment_cigar_ops:
        if cigar_op in ('M', '=', 'X'):
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                read_aln_pos = current_read_pos_in_read + i
                if window_start_node <= node_aln_pos < window_start_node + window_size:
                    window_idx = node_aln_pos - window_start_node
                    if read_aln_pos < len(segment_read_sequence):
                        window_char_representation[window_idx] = segment_read_sequence[read_aln_pos].upper()
            current_node_pos_in_read += cigar_len
            current_read_pos_in_read += cigar_len
        elif cigar_op == 'D' or cigar_op == 'N':
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                if window_start_node <= node_aln_pos < window_start_node + window_size:
                    window_idx = node_aln_pos - window_start_node
                    window_char_representation[window_idx] = '*'  # Deletions shown as '*'
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I' or cigar_op == 'S':
            current_read_pos_in_read += cigar_len
        if current_node_pos_in_read >= window_start_node + window_size:  # Optimization
            break
    return window_char_representation


def get_read_tensor_rows_in_window(segment_cigar_ops, segment_offset_on_node,
                                   segment_read_sequence, segment_quality_str,
                                   window_start_node, current_tensor_window_size, node_len):
    # PADDING_BASE_INDEX (index 6) is used for bases. DEFAULT_QUALITY_PADDING (0) for qualities.
    base_indices_row = [PADDING_BASE_INDEX] * current_tensor_window_size
    quality_scores_row = [DEFAULT_QUALITY_PADDING] * current_tensor_window_size
    current_node_pos_in_read = segment_offset_on_node
    current_read_pos_in_read = 0
    for cigar_len, cigar_op in segment_cigar_ops:
        if current_node_pos_in_read >= window_start_node + current_tensor_window_size and cigar_op in (
                'M', 'D', 'N', '=', 'X'):  # Optimization
            break
        if current_read_pos_in_read >= len(segment_read_sequence):  # Boundary check
            break

        if cigar_op in ('M', '=', 'X'):
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                read_aln_pos = current_read_pos_in_read + i
                if read_aln_pos >= len(segment_read_sequence): break  # Boundary check for read
                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_char = segment_read_sequence[read_aln_pos].upper()
                    base_indices_row[window_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                    if read_aln_pos < len(segment_quality_str):
                        try:
                            quality_scores_row[window_idx] = ord(segment_quality_str[read_aln_pos]) - 33
                        except (TypeError, IndexError):
                            quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
                    else:
                        quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
            current_node_pos_in_read += cigar_len
            current_read_pos_in_read += cigar_len
        elif cigar_op == 'D' or cigar_op == 'N':  # Deletion or Skip in reference
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_indices_row[window_idx] = BASE_TO_INDEX['*']  # Deletion marker
                    quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I' or cigar_op == 'S':  # Insertion or Soft clip in read
            current_read_pos_in_read += cigar_len
    return base_indices_row, quality_scores_row


## ─────────────────────────────────────────────────────────────────────────────
## Worker Process Initialization and Target Function
## ─────────────────────────────────────────────────────────────────────────────

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}\n")
        sys.exit(1)


def process_single_node_for_pileup(task_args_with_af_thresh):
    node_id, dat_file_offset, n_records, node_sequence, min_af_threshold = task_args_with_af_thresh
    global worker_dat_file, worker_base_output_dir

    if worker_dat_file is None or worker_base_output_dir is None:
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()}]: Worker not initialized properly (dat_file or output_dir missing).\n")
        return node_id, None

    if not node_sequence:
        sys.stderr.write(
            f"ℹ️ [Worker {os.getpid()}]: No sequence provided for node {node_id}. Skipping tensor/summary generation.\n")
        return node_id, {}

    node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
    try:
        os.makedirs(node_specific_output_dir, exist_ok=True)
    except OSError as e:
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()}]: Could not create directory {node_specific_output_dir}: {e}\n")
        return node_id, None

    node_len = len(node_sequence)
    variant_headers_for_node = []
    view_oriented_variant_data = {}

    aligned_read_segments = []
    try:
        worker_dat_file.seek(dat_file_offset + 10)
        for record_idx in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE: break
            off_from_file, raw_seq, raw_qual, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)
            if mapq < 10: continue
            try:
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                qual_str = raw_qual.rstrip(b'\x00').decode('ascii', errors='replace')
                cigar_str_original = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue
            if len(seq) == 0 or len(seq) != len(qual_str):
                continue

            original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
            current_read_sequence = seq
            current_quality_str = qual_str
            current_decoded_cigar_ops = list(original_decoded_cigar_ops)
            current_offset_on_node = off_from_file

            if strand_char == '-':
                current_read_sequence = reverse_complement(seq)
                current_quality_str = qual_str[::-1]
                current_decoded_cigar_ops.reverse()
                alignment_span_on_node = len(current_read_sequence)  # User-confirmed fix
                if alignment_span_on_node > 0:
                    # The CIGAR defines alignment span on the reference.
                    # For reverse strand, if original offset was X, and original CIGAR aligns for L bases on node,
                    # the segment on node is [X, X+L-1].
                    # On reversed node sequence (conceptually), this would be [node_len-(X+L-1)-1, node_len-X-1].
                    # The new offset is node_len - (off_from_file + alignment_span_on_node)
                    current_offset_on_node = node_len - (off_from_file + alignment_span_on_node)
                    if current_offset_on_node < 0:
                        # This can happen if GFA node length is inconsistent with alignment reference length
                        # Or if alignment_span calculation from CIGAR is imperfect for the specific case.
                        # sys.stderr.write(f"⚠️ Warning [Worker {os.getpid()}]: Negative offset for node {node_id}, strand '-', orig_off {off_from_file}, span {alignment_span_on_node}, node_len {node_len}. Skipping.\n")
                        continue

            aligned_read_segments.append({
                "offset_on_node": current_offset_on_node,
                "read_sequence": current_read_sequence,
                "processed_quality_str": current_quality_str,
                "cigar_ops": current_decoded_cigar_ops,
                "original_cigar_str": cigar_str_original,
                "strand": strand_char
            })
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}] reading records for node {node_id}: {e}\n")
        return node_id, None

    candidate_variants = defaultdict(int)
    for segment in aligned_read_segments:
        variants_in_read = detect_variants_from_cigar(
            segment["offset_on_node"], segment["cigar_ops"],
            segment["read_sequence"], node_sequence
        )
        for v_pos, v_type, v_alt, v_ref in variants_in_read:
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    variant_processing_window_size = TENSOR_WINDOW_SIZE  # This is the 'window_size' for get_read_representation
    half_window = variant_processing_window_size // 2

    for (v_pos, v_type, v_ref_defined, v_alt_defined), _ in candidate_variants.items():
        af_alt_count, af_ref_count, af_other_count, af_locus_coverage = 0, 0, 0, 0
        ref_allele_for_indel_af_check = v_ref_defined if v_type == 'D' else \
            (node_sequence[v_pos] if v_type == 'I' and 0 <= v_pos < node_len else None)
        for segment in aligned_read_segments:
            allele = get_allele_from_read_at_node_pos(
                segment["offset_on_node"], segment["read_sequence"], segment["cigar_ops"],
                v_pos, node_sequence, v_type, ref_allele_for_indel_af_check)
            if allele is not None:
                af_locus_coverage += 1
                if v_type == 'X':
                    if allele == v_alt_defined:
                        af_alt_count += 1
                    elif allele == v_ref_defined:
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'I':
                    if allele == v_alt_defined:
                        af_alt_count += 1
                    elif allele == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'D':
                    if allele == "*":  # Alt is deletion
                        af_alt_count += 1
                    elif allele == "REF_STATE_FOR_INDEL":  # Ref state means read spans the deleted region
                        af_ref_count += 1
                    else:  # e.g. other indel
                        af_other_count += 1

        alt_freq = af_alt_count / af_locus_coverage if af_locus_coverage > 0 else 0.0
        if alt_freq < min_af_threshold: continue

        variant_key_str = f"{v_pos}_{v_type}_{v_ref_defined}_{v_alt_defined}"

        pileup_reads_data_for_view = []
        # Determine window center for view/tensor (variant position itself, or after for insertion)
        if v_type == 'I':  # For insertions, anchor base is v_pos, insertion occurs after v_pos
            window_center_on_node_view = v_pos + 1
        else:  # For SNPs, Deletions, v_pos is the first affected coordinate
            window_center_on_node_view = v_pos

        # This start can be negative or window can extend beyond node_len. Handled by downstream functions.
        window_start_on_node_view = window_center_on_node_view - half_window

        for segment in aligned_read_segments:
            # `get_read_representation_in_window_for_view` uses ' ' for padding,
            # which will be converted to PADDING_BASE_INDEX (6)
            row_chars = get_read_representation_in_window_for_view(
                segment["cigar_ops"], segment["offset_on_node"], segment["read_sequence"],
                window_start_on_node_view, variant_processing_window_size, node_len)

            if any(c != ' ' for c in row_chars):  # Only add if read has some overlap with window
                row_indices = [BASE_TO_INDEX.get(char.upper(), PADDING_BASE_INDEX) for char in row_chars]
                pileup_reads_data_for_view.append({
                    "bases": row_indices, "offset": segment["offset_on_node"],
                    "strand": segment["strand"], "cigar": segment["original_cigar_str"]})

        view_oriented_variant_data[variant_key_str] = {
            "pileup_reads_data": pileup_reads_data_for_view,
            "alt_allele_count": af_alt_count, "ref_allele_count_at_locus": af_ref_count,
            "other_allele_count_at_locus": af_other_count, "coverage_at_locus": af_locus_coverage,
            "alt_allele_frequency": round(alt_freq, 4)
        }

        tensor_ch1_bases_list = []
        tensor_ch2_qualities_list = []
        tensor_ch3_mismatches_list = []

        # Tensor window uses the same start as the view window.
        window_start_on_node_tensor = window_start_on_node_view

        # Reference row for tensor
        ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE  # PADDING_BASE_INDEX (6) for out-of-bounds
        ref_qual_scores_row = [DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE  # CH2 Ref Quals (0)
        ref_mismatch_row = [MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE  # CH3 Ref Mismatches (0)
        for i in range(TENSOR_WINDOW_SIZE):
            actual_node_pos = window_start_on_node_tensor + i
            if 0 <= actual_node_pos < node_len:  # Check if within node bounds
                ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[actual_node_pos].upper(), BASE_TO_INDEX['N'])
            # Else, it remains PADDING_BASE_INDEX
        tensor_ch1_bases_list.append(ref_base_indices_row)
        tensor_ch2_qualities_list.append(ref_qual_scores_row)
        tensor_ch3_mismatches_list.append(ref_mismatch_row)

        reads_added_to_tensor = 0
        for segment in aligned_read_segments:  # Iterate through all reads relevant to the node
            if reads_added_to_tensor >= TENSOR_MAX_READ_ROWS: break

            # get_read_tensor_rows_in_window will use PADDING_BASE_INDEX for bases and DEFAULT_QUALITY_PADDING for qualities
            # for parts of the window not covered by the read or outside the node.
            base_indices_row, quality_scores_row = get_read_tensor_rows_in_window(
                segment["cigar_ops"], segment["offset_on_node"],
                segment["read_sequence"], segment["processed_quality_str"],
                window_start_on_node_tensor, TENSOR_WINDOW_SIZE, node_len)

            # Only add if there's some actual base data, not all padding
            # (though get_read_tensor_rows_in_window doesn't explicitly filter empty, rely on downstream or CIGAR)
            # This check might be implicit if all padding leads to no mismatch contribution.

            tensor_ch1_bases_list.append(base_indices_row)
            tensor_ch2_qualities_list.append(quality_scores_row)

            mismatch_row_for_read = [MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE  # Init with 0
            for i in range(TENSOR_WINDOW_SIZE):
                read_base_idx = base_indices_row[i]
                ref_base_idx_for_comp = ref_base_indices_row[i]  # From the ref_base_indices_row defined above

                # If either read or ref position is padding, mismatch is padding value (0)
                if read_base_idx == PADDING_BASE_INDEX or ref_base_idx_for_comp == PADDING_BASE_INDEX:
                    mismatch_row_for_read[i] = MISMATCH_COMPARISON_PADDING_VALUE
                elif read_base_idx == ref_base_idx_for_comp:
                    mismatch_row_for_read[i] = 0  # Match
                else:
                    # Ensure not to compare against N if possible, though N vs N is match.
                    # '*' (deletion in read against ref base) vs ref base = mismatch
                    # ref base vs '*' (deletion in ref / read has insertion) = mismatch potentially, depends on definition
                    # Current PADDING_BASE_INDEX is 6, '*' is 5. 'N' is 4.
                    if read_base_idx != BASE_TO_INDEX['N'] and ref_base_idx_for_comp != BASE_TO_INDEX['N']:
                        mismatch_row_for_read[i] = 1  # Mismatch
                    else:  # One is N, other is not N and not padding. Treat as mismatch or specific N value.
                        mismatch_row_for_read[i] = 1  # Simplification: N is a mismatch to a concrete base.
            tensor_ch3_mismatches_list.append(mismatch_row_for_read)
            reads_added_to_tensor += 1

        # Pad tensor to TENSOR_MAX_READ_ROWS + 1 (for ref)
        for _ in range(TENSOR_MAX_READ_ROWS - reads_added_to_tensor):
            tensor_ch1_bases_list.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            tensor_ch2_qualities_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            tensor_ch3_mismatches_list.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)

        try:
            final_tensor = torch.tensor([
                tensor_ch1_bases_list,
                tensor_ch2_qualities_list,
                tensor_ch3_mismatches_list
            ], dtype=torch.int8)

            tensor_filename = f"{variant_key_str}.pth"
            tensor_filepath = os.path.join(node_specific_output_dir, tensor_filename)
            torch.save(final_tensor, tensor_filepath)

            variant_headers_for_node.append({
                "variant_key": variant_key_str,
                "tensor_file": tensor_filename,
                "alt_allele_count": af_alt_count,
                "ref_allele_count_at_locus": af_ref_count,
                "other_allele_count_at_locus": af_other_count,
                "coverage_at_locus": af_locus_coverage,
                "alt_allele_frequency": round(alt_freq, 4)
            })
        except Exception as e:
            sys.stderr.write(
                f"❌ Error [Worker {os.getpid()}]: Failed to create or save tensor for {variant_key_str} in node {node_id}: {e}\n")

    if variant_headers_for_node:
        summary_json_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        try:
            with open(summary_json_path, 'w') as sjf:
                json.dump({
                    "node_id": node_id,
                    "node_length": node_len,
                    "node_sequence": node_sequence,
                    "variants": variant_headers_for_node
                }, sjf, indent=2)
        except Exception as e:
            sys.stderr.write(f"❌ Error [Worker {os.getpid()}]: Failed to write summary JSON for node {node_id}: {e}\n")

    return node_id, view_oriented_variant_data


## ─────────────────────────────────────────────────────────────────────────────
## Pileup Viewing Function (display_pileup_data)
## ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if not node_data_for_display_view or not isinstance(node_data_for_display_view, dict):
        print(f"ℹ️ No valid pileup data to display for node {node_id_str_for_display}.", file=sys.stderr)
        return

    print(
        f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")

    if not node_data_for_display_view:
        print(f"ℹ️ No variants found or pileups generated for this node (or all filtered by AF).")
        return

    variants_displayed_count = 0
    # Sort by position then type for consistent display order
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))

    display_window_size = TENSOR_WINDOW_SIZE  # Use same window size as tensor for consistency
    half_display_window = display_window_size // 2
    padding_char_for_ref_display = '0'  # MODIFIED: Character for padding ref display

    for variant_key in sorted_variant_keys:
        if variants_displayed_count >= max_variants_to_display:
            print(
                f"\n  ... (and {len(node_data_for_display_view) - variants_displayed_count} more variants not shown due to limit)")
            break

        variant_data = node_data_for_display_view[variant_key]
        pileup_reads_display_data = variant_data.get("pileup_reads_data", [])

        v_pos = int(variant_key.split('_')[0])
        v_type = variant_key.split('_')[1]

        # Determine window center (variant position itself, or after for insertion)
        if v_type == 'I':
            window_center_on_node = v_pos + 1
        else:
            window_center_on_node = v_pos

        # This start can be negative. The display logic below handles it.
        current_window_start_on_node = window_center_on_node - half_display_window

        print(
            f"\n--- Variant: {variant_key} (Node Pos: {v_pos}, Display Window on Node: {current_window_start_on_node}-{current_window_start_on_node + display_window_size - 1}) ---")
        print(f"  Alt Count: {variant_data.get('alt_allele_count', 'N/A')}")
        print(f"  Ref Count: {variant_data.get('ref_allele_count_at_locus', 'N/A')}")
        print(f"  Other Count: {variant_data.get('other_allele_count_at_locus', 'N/A')}")
        print(f"  Coverage: {variant_data.get('coverage_at_locus', 'N/A')}")
        alt_freq_val = variant_data.get('alt_allele_frequency', 'N/A')
        if isinstance(alt_freq_val, float):
            print(f"  Alt Freq: {alt_freq_val:.4f}")
        else:
            print(f"  Alt Freq: {alt_freq_val}")

        # MODIFIED: Initialize reference display with padding_char_for_ref_display ('0')
        ref_display_parts = [padding_char_for_ref_display] * display_window_size
        marker_line_parts = [' '] * display_window_size  # Keep markers as spaces initially

        variant_display_idx_in_window = -1  # Relative to the window 0..display_window_size-1
        # Calculate where v_pos itself (or anchor for I) falls in the current window
        # v_pos is the 0-indexed coordinate on the node.
        # current_window_start_on_node is the 0-indexed node coordinate at the start of the window.
        # So, relative index is v_pos - current_window_start_on_node

        # For marker, use v_pos as the primary coordinate of interest.
        # For insertions, v_pos is the base *before* the insertion. The marker points between v_pos and v_pos+1.
        # Visually, if window shows v_pos, marker could be at v_pos or v_pos+1 depending on convention.
        # The script previously put 'I' at v_pos and '^' at v_pos+1 if space. Let's keep this.

        variant_marker_node_coord = v_pos  # The actual node coordinate the marker is related to.

        if current_window_start_on_node <= variant_marker_node_coord < current_window_start_on_node + display_window_size:
            variant_display_idx_in_window = variant_marker_node_coord - current_window_start_on_node

        for i in range(display_window_size):
            actual_node_pos_in_window = current_window_start_on_node + i
            if 0 <= actual_node_pos_in_window < len(full_node_sequence):  # Check if within actual node sequence
                ref_display_parts[i] = full_node_sequence[actual_node_pos_in_window].upper()
            # Else, it remains padding_char_for_ref_display ('0')

            if i == variant_display_idx_in_window:  # If current window slot is the variant's position
                if v_type == 'I':
                    marker_line_parts[i] = "I"  # Mark the base *before* insertion
                    # Place '^' under the conceptual gap where insertion occurs, if space.
                    if i + 1 < display_window_size: marker_line_parts[i + 1] = "^"
                    # If insertion is at the very end of window, only 'I' might show.
                else:  # SNP, DEL
                    marker_line_parts[i] = "^"  # Mark the affected base(s)

        print(f"  Node Ref: {''.join(ref_display_parts)}")
        print(f"  Marker  : {''.join(marker_line_parts)}")

        if not pileup_reads_display_data:
            print("  (No reads data in pileup for this variant's window or all filtered)")
        else:
            displayed_reads_count = 0
            for i, read_info in enumerate(pileup_reads_display_data):
                if displayed_reads_count >= max_reads_to_display_per_variant:
                    print(f"  ... (and {len(pileup_reads_display_data) - displayed_reads_count} more reads)")
                    break
                base_indices = read_info["bases"]  # These are indices
                read_offset = read_info["offset"]
                read_strand = read_info["strand"]
                read_cigar = read_info.get("cigar", "N/A")
                # Convert indices to viewable bases. PADDING_BASE_INDEX (6) will map to '0'.
                pileup_row_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in base_indices])
                print(
                    f"  Read {i + 1:3d}: {pileup_row_str}  (Offset: {read_offset}, Strand: {read_strand}, CIGAR: {read_cigar})")
                displayed_reads_count += 1
        variants_displayed_count += 1
    print("\n")


## ─────────────────────────────────────────────────────────────────────────────
## Main function
## ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered pileups for a single specified node and optionally view them.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="Base output directory for node-specific folders, tensors, and summaries.")
    parser.add_argument("--node_id", type=int, required=True, help="The specific node ID to process.")
    parser.add_argument("--gfa", help="GFA graph file path (required if node sequence cache is not used/built).")
    parser.add_argument("--load-cache", help="Load node sequence from this JSON cache file.")
    parser.add_argument("--save-cache", help="Save node sequence to this JSON cache file (used if --gfa is provided).")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                        help="Print generated pileups to console. Optionally specify N for first N variants (e.g., --view 5). If no number, shows all. If 0, shows none from view block.")
    parser.add_argument("--max_view_reads", type=int, default=20,
                        help="Max reads to display per pileup in console view.")
    parser.add_argument("--min_af", type=float, default=0.1,
                        help="Min allele frequency for a variant to be processed.")
    args = parser.parse_args()

    if not os.path.isfile(args.dat): sys.exit(f"❌ Error: DAT file not found: {args.dat}")
    if not os.path.isfile(args.idx): sys.exit(f"❌ Error: Index file not found: {args.idx}")
    if not args.load_cache and not args.gfa: sys.exit("❌ Error: Must provide --gfa or --load-cache.")
    if args.load_cache and not os.path.isfile(args.load_cache): sys.exit(
        f"❌ Error: Cache file not found: {args.load_cache}")
    if args.gfa and not os.path.isfile(args.gfa): sys.exit(f"❌ Error: GFA file not found: {args.gfa}")
    if args.min_af < 0.0 or args.min_af > 1.0: sys.exit("❌ Error: --min_af must be between 0.0 and 1.0.")

    try:
        os.makedirs(args.output, exist_ok=True)
        print(f"🔹 Base output directory: {args.output}")
    except OSError as e:
        sys.exit(f"❌ Error: Could not create base output directory {args.output}: {e}")

    target_node_id = args.node_id
    print(f"🔹 Processing single target node ID: {target_node_id}")

    start_time = time.time()
    node_dat_info = parse_idx_file_for_single_node(args.idx, target_node_id)
    if not node_dat_info: sys.exit(f"❌ Error: Failed to get index info for node {target_node_id}.")
    dat_offset, n_records = node_dat_info
    print(f"✔ Index parsing for node {target_node_id} took {time.time() - start_time:.2f}s.")

    node_sequence = None
    if args.load_cache and os.path.isfile(args.load_cache):
        start_time_cache = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                loaded_cache = json.load(cf)
            node_sequence = loaded_cache.get(str(target_node_id))
            if node_sequence:
                print(
                    f"✔ Loaded sequence for node {target_node_id} from cache in {time.time() - start_time_cache:.2f}s.")
            else:
                print(f"⚠️ Warning: Node {target_node_id} not in cache {args.load_cache}.")
                if not args.gfa: sys.exit(f"❌ Error: Node {target_node_id} not in cache and no GFA. Exiting.")
        except Exception as e:
            print(f"❌ Error loading cache {args.load_cache}: {e}", file=sys.stderr)
            if not args.gfa: sys.exit(1)
            node_sequence = None  # Ensure it's None if cache fails and GFA is required.

    if node_sequence is None and args.gfa:  # If not loaded from cache (or cache failed) and GFA is provided
        start_time_gfa = time.time()
        node_sequence = load_node_sequence_from_gfa(args.gfa, target_node_id)
        if node_sequence and args.save_cache:
            print(f"🔹 Saving sequence for node {target_node_id} to cache: {args.save_cache}...")
            try:
                existing_cache_data = {}
                if os.path.isfile(args.save_cache):
                    with open(args.save_cache, 'r') as rcf:
                        try:
                            existing_cache_data = json.load(rcf)
                        except json.JSONDecodeError:
                            print(f"⚠️ Corrupt cache {args.save_cache}, overwriting.", file=sys.stderr)
                existing_cache_data[str(target_node_id)] = node_sequence
                with open(args.save_cache, 'w') as wcf:
                    json.dump(existing_cache_data, wcf, indent=2)  # Added indent
                print(f"✔ Saved sequence for node {target_node_id} to cache.")
            except Exception as e:
                print(f"❌ Error saving to cache {args.save_cache}: {e}", file=sys.stderr)
        elif node_sequence:  # Loaded from GFA but not saving to cache
            print(f"✔ Sequence loading for node {target_node_id} from GFA took {time.time() - start_time_gfa:.2f}s.")

    if not node_sequence: sys.exit(f"❌ Error: Failed to obtain sequence for node {target_node_id}. Exiting.")

    task = (target_node_id, dat_offset, n_records, node_sequence, args.min_af)
    print(f"🔹 Prepared task for node {target_node_id} with min AF threshold: {args.min_af}.")
    print(f"🔹 Processing node {target_node_id} using 1 worker (ProcessPoolExecutor)...")
    start_proc_time = time.time()

    processed_node_id_result = None
    view_data_for_node = None

    try:
        # Using ProcessPoolExecutor for a single task here is a bit overhead but matches original structure
        # It's useful if script were to be expanded to multiple nodes in main again.
        with ProcessPoolExecutor(max_workers=1, initializer=init_worker, initargs=(args.dat, args.output)) as executor:
            future = executor.submit(process_single_node_for_pileup, task)
            # Add timeout to future.result() if concerned about hangs, though for single node it might be less critical
            processed_node_id_result, view_data_for_node = future.result()

        if view_data_for_node is None and processed_node_id_result is not None:
            print(
                f"⚠️ Warning: Processing for node {processed_node_id_result} might have encountered an issue or yielded no processable variants. Check worker logs. No view data returned.",
                file=sys.stderr)
        elif processed_node_id_result is None:  # Indicates a more severe failure in the worker
            print(
                f"⚠️ Critical Error: Worker did not return expected results for node {target_node_id}. Check for errors above.",
                file=sys.stderr)

    except Exception as pool_exc:  # Catch exceptions from submit or result
        sys.stderr.write(
            f"\n❌ An error occurred during processing node {target_node_id} via ProcessPoolExecutor: {pool_exc}\n")
        # Depending on severity, may want to sys.exit here or just note the failure.
        # For single node processing, exiting might be appropriate.
        sys.exit(1)

    total_elapsed_time = time.time() - start_proc_time
    print(f"✔ Node {target_node_id} processing (including I/O if any in worker) finished in {total_elapsed_time:.2f}s.")

    if args.view is not None:  # if --view is passed
        if view_data_for_node and isinstance(view_data_for_node, dict) and view_data_for_node:
            max_v_show = float('inf')
            if args.view == -1:  # --view without a number
                max_v_show = float('inf')
            elif args.view >= 0:  # --view N (N can be 0)
                max_v_show = args.view

            view_msg_detail = f"(max {args.max_view_reads} reads/variant)"
            if max_v_show == float('inf'):
                view_msg = f"all variants found {view_msg_detail}..."
            elif max_v_show == 0:
                view_msg = "0 variants (view suppressed by --view 0)."
            else:
                view_msg = f"first {int(max_v_show)} variants {view_msg_detail}..."

            print(f"🔹 Displaying pileups for node {target_node_id}: {view_msg}")
            if max_v_show > 0:
                display_pileup_data(view_data_for_node, str(target_node_id),
                                    node_sequence, args.max_view_reads, max_v_show)
        elif processed_node_id_result is not None:  # Node was processed, but no view_data (e.g. no variants met AF)
            print(
                f"\nℹ️ --view specified for node {target_node_id}, but no variants met AF threshold or were found to display (or an error occurred in worker resulting in no view data).")
            if node_sequence: print(f"   (Node Length: {len(node_sequence)})")
        # If processed_node_id_result is None, an error message about worker failure was already printed.

    node_summary_file = os.path.join(args.output, str(target_node_id), "variant_summary.json")
    if os.path.exists(node_summary_file):
        print(f"✅ Output generated for node {target_node_id} in {os.path.join(args.output, str(target_node_id))}")
    elif processed_node_id_result is not None and (view_data_for_node is None or not view_data_for_node):
        # This case covers:
        # 1. view_data_for_node is None (worker returned it as such, implying an issue or no variants for tensor/view)
        # 2. view_data_for_node is an empty dict (no variants met AF for view data generation)
        # process_single_node_for_pileup might still write a summary if variant_headers_for_node is empty,
        # but if variant_headers_for_node is also empty, no summary.json will be written.
        print(
            f"ℹ️ No variants met AF threshold for node {target_node_id}, or an issue prevented tensor/summary generation. No summary file at {node_summary_file}.")
    elif processed_node_id_result is None:  # Worker failed more substantially.
        print(f"❌ Critical error in worker for node {target_node_id}. No output summary expected.")
    else:  # Fallback, should ideally be covered by above.
        print(f"ℹ️ Output summary status for node {target_node_id} is unclear. Check logs and output directory.")

    print("✅ Script finished.")


if __name__ == '__main__':
    main()
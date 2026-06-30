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
# Read offset, sequence, RAW QUALITIES, CIGAR, MAPQ, strand
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size

BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 5, ' ': 6, '-': 6}  # Added '-' to map to space/padding
INDEX_TO_BASE_FOR_VIEW = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N', 5: '*', 6: ' '}

# Tensor-specific constants
TENSOR_WINDOW_SIZE = 100  # Should match 'window_size' used for variant processing
TENSOR_MAX_READ_ROWS = 200  # Max number of read rows in the tensor
PADDING_BASE_INDEX = BASE_TO_INDEX[' ']  # Numerical index for padding bases
DEFAULT_QUALITY_PADDING = 0  # Default quality for padding or non-base characters

# Global for worker process state (file handle)
worker_dat_file = None


# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions

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
    """Generates character representation for console view."""
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
                    window_char_representation[window_idx] = '*'
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I' or cigar_op == 'S':
            current_read_pos_in_read += cigar_len
        if current_node_pos_in_read >= window_start_node + window_size:
            break
    return window_char_representation


def get_read_tensor_rows_in_window(segment_cigar_ops, segment_offset_on_node,
                                   segment_read_sequence, segment_quality_str,
                                   window_start_node, current_tensor_window_size, node_len):
    """
    Generates numerical base and quality representations for a read in a window for tensor output.
    """
    base_indices_row = [PADDING_BASE_INDEX] * current_tensor_window_size
    quality_scores_row = [DEFAULT_QUALITY_PADDING] * current_tensor_window_size

    current_node_pos_in_read = segment_offset_on_node
    current_read_pos_in_read = 0

    for cigar_len, cigar_op in segment_cigar_ops:
        if current_node_pos_in_read >= window_start_node + current_tensor_window_size and cigar_op in (
        'M', 'D', 'N', '=', 'X'):
            break
        if current_read_pos_in_read >= len(segment_read_sequence):  # Consumed entire read
            break

        if cigar_op in ('M', '=', 'X'):
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                read_aln_pos = current_read_pos_in_read + i

                if read_aln_pos >= len(segment_read_sequence): break  # Safety break

                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_char = segment_read_sequence[read_aln_pos].upper()
                    base_indices_row[window_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])

                    if read_aln_pos < len(segment_quality_str):
                        try:
                            quality_scores_row[window_idx] = ord(segment_quality_str[read_aln_pos]) - 33
                        except TypeError:  # Handle if qual char is not an ord-able type somehow
                            quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
                        except IndexError:  # Should not happen if lengths match
                            quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
                    else:  # Should not happen if quality string has same length as sequence
                        quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING
            current_node_pos_in_read += cigar_len
            current_read_pos_in_read += cigar_len
        elif cigar_op == 'D' or cigar_op == 'N':  # Consumes reference
            for i in range(cigar_len):
                node_aln_pos = current_node_pos_in_read + i
                if window_start_node <= node_aln_pos < window_start_node + current_tensor_window_size:
                    window_idx = node_aln_pos - window_start_node
                    base_indices_row[window_idx] = BASE_TO_INDEX['*']
                    quality_scores_row[window_idx] = DEFAULT_QUALITY_PADDING  # No base, so no quality
            current_node_pos_in_read += cigar_len
        elif cigar_op == 'I' or cigar_op == 'S':  # Consumes query (read and qualities)
            current_read_pos_in_read += cigar_len

    return base_indices_row, quality_scores_row


# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function

def init_worker(dat_file_path_for_worker):
    global worker_dat_file
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
    except FileNotFoundError:
        print(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}", file=sys.stderr)
        sys.exit(1)


def process_single_node_for_pileup(task_args_with_af_thresh):
    node_id, dat_file_offset, n_records, node_sequence, min_af_threshold = task_args_with_af_thresh
    global worker_dat_file

    if worker_dat_file is None: return node_id, {}
    if not node_sequence: return node_id, {}

    node_len = len(node_sequence)
    final_variant_output = {}
    aligned_read_segments = []

    try:
        worker_dat_file.seek(dat_file_offset + 10)
        for record_idx in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE: break
            # Unpack including raw_qual
            off_from_file, raw_seq, raw_qual, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)
            if mapq < 10: continue
            try:
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='replace')
                qual_str = raw_qual.rstrip(b'\x00').decode('ascii', errors='replace')  # Decode qualities
                cigar_str_original = raw_cigar.rstrip(b'\x00').decode('ascii', errors='replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue

            if len(seq) != len(qual_str):  # Basic sanity check
                # print(f"⚠️ Warning: Mismatch length seq vs qual for a read on node {node_id}. Skipping read.")
                # print(f"   Seq len: {len(seq)}, Qual len: {len(qual_str)}")
                continue

            original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)

            current_read_sequence = seq
            current_quality_str = qual_str  # Store qualities
            current_decoded_cigar_ops = list(original_decoded_cigar_ops)
            current_offset_on_node = off_from_file

            if strand_char == '-':
                current_read_sequence = reverse_complement(seq)
                current_quality_str = qual_str[::-1]  # Reverse qualities
                current_decoded_cigar_ops.reverse()

                alignment_span_on_node = len(current_read_sequence)  # User's confirmed fix

                # This transformation block is from the user-confirmed working script
                if alignment_span_on_node > 0:  # Ensure span is positive
                    current_offset_on_node = node_len - alignment_span_on_node - off_from_file
                    if current_offset_on_node < 0:
                        # print(f"⚠️ Negative offset node {node_id}, strand '-', off_from_file {off_from_file}, span {alignment_span_on_node}. Skipping.")
                        continue
                        # The conditional -1 based on off_from_file !=0 was discussed but not in the last script the user provided as "fixed"
                # If that was intended, it would be applied here to current_offset_on_node.
                # For now, sticking to the last provided script logic for this part.

            aligned_read_segments.append({
                "offset_on_node": current_offset_on_node,
                "read_sequence": current_read_sequence,
                "processed_quality_str": current_quality_str,  # Store processed qualities
                "cigar_ops": current_decoded_cigar_ops,
                "original_cigar_str": cigar_str_original,
                "strand": strand_char
            })
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] reading records for node {node_id}: {e}", file=sys.stderr)
        return node_id, {}

    candidate_variants = defaultdict(int)
    # This part remains unchanged, uses segment["cigar_ops"], segment["read_sequence"], etc.
    for segment in aligned_read_segments:
        variants_in_read = detect_variants_from_cigar(
            segment["offset_on_node"], segment["cigar_ops"],
            segment["read_sequence"], node_sequence
        )
        for v_pos, v_type, v_alt, v_ref in variants_in_read:
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    # This should be TENSOR_WINDOW_SIZE for consistency
    # window_size was used for display_pileup_data marker logic, ensure it's consistent
    # For tensor generation, use TENSOR_WINDOW_SIZE explicitly
    variant_processing_window_size = TENSOR_WINDOW_SIZE
    half_window = variant_processing_window_size // 2

    for (v_pos, v_type, v_ref_defined, v_alt_defined), _ in candidate_variants.items():
        # ... (AF calculation logic remains the same) ...
        af_alt_count = 0
        af_ref_count = 0
        af_other_count = 0
        af_locus_coverage = 0
        ref_allele_for_indel_af_check = v_ref_defined if v_type == 'D' else \
            (node_sequence[v_pos] if v_type == 'I' and 0 <= v_pos < node_len else None)
        for segment in aligned_read_segments:
            allele_in_segment_at_v_pos = get_allele_from_read_at_node_pos(
                segment["offset_on_node"], segment["read_sequence"], segment["cigar_ops"],
                v_pos, node_sequence, v_type, ref_allele_for_indel_af_check
            )
            if allele_in_segment_at_v_pos is not None:
                af_locus_coverage += 1
                if v_type == 'X':
                    if allele_in_segment_at_v_pos == v_alt_defined:
                        af_alt_count += 1
                    elif allele_in_segment_at_v_pos == v_ref_defined:
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'I':
                    if allele_in_segment_at_v_pos == v_alt_defined:
                        af_alt_count += 1
                    elif allele_in_segment_at_v_pos == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1
                    else:
                        af_other_count += 1
                elif v_type == 'D':
                    if allele_in_segment_at_v_pos == "*":
                        af_alt_count += 1
                    elif allele_in_segment_at_v_pos == "REF_STATE_FOR_INDEL":
                        af_ref_count += 1
                    else:
                        af_other_count += 1
        alt_freq = af_alt_count / af_locus_coverage if af_locus_coverage > 0 else 0.0
        if alt_freq < min_af_threshold:
            continue

        # --- Prepare data for console view (pileup_reads_data) ---
        # This data is used by display_pileup_data
        pileup_reads_data_for_view = []
        if v_type == 'I':
            window_center_on_node_view = v_pos + 1
        else:
            window_center_on_node_view = v_pos
        window_start_on_node_view = window_center_on_node_view - half_window  # uses variant_processing_window_size's half_window

        for segment in aligned_read_segments:
            row_chars_for_view = get_read_representation_in_window_for_view(
                segment["cigar_ops"], segment["offset_on_node"], segment["read_sequence"],
                window_start_on_node_view, variant_processing_window_size, node_len
            )
            if any(c != ' ' for c in row_chars_for_view):
                row_indices_for_view = [BASE_TO_INDEX.get(char.upper(), BASE_TO_INDEX['N']) for char in
                                        row_chars_for_view]
                pileup_reads_data_for_view.append({
                    "bases": row_indices_for_view,  # Numerical indices for console view compatibility
                    "offset": segment["offset_on_node"],
                    "strand": segment["strand"],
                    "cigar": segment["original_cigar_str"]
                })

        # --- Tensor Generation ---
        tensor_ch1_bases = []
        tensor_ch2_qualities = []

        # Determine window start for tensor (same as for view for consistency)
        window_start_on_node_tensor = window_start_on_node_view

        # Channel 1 & 2: Reference Row
        ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        ref_qual_scores_row = [DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE
        for i in range(TENSOR_WINDOW_SIZE):
            actual_node_pos = window_start_on_node_tensor + i
            if 0 <= actual_node_pos < node_len:
                ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[actual_node_pos].upper(), BASE_TO_INDEX['N'])
            # else it remains PADDING_BASE_INDEX (and PADDING_QUALITY_VALUE)
        tensor_ch1_bases.append(ref_base_indices_row)
        tensor_ch2_qualities.append(ref_qual_scores_row)

        # Channel 1 & 2: Read Rows
        reads_added_to_tensor = 0
        for segment in aligned_read_segments:
            if reads_added_to_tensor >= TENSOR_MAX_READ_ROWS:
                break

            # Generate numerical base and quality rows for the tensor
            base_indices_row, quality_scores_row = get_read_tensor_rows_in_window(
                segment["cigar_ops"], segment["offset_on_node"],
                segment["read_sequence"], segment["processed_quality_str"],
                window_start_on_node_tensor, TENSOR_WINDOW_SIZE, node_len
            )

            # Add if the read has any coverage in the window (optional check, could always add)
            # For simplicity, we'll add it and let padding rows fill if it's all non-coverage later
            # if any(b != PADDING_BASE_INDEX for b in base_indices_row):
            tensor_ch1_bases.append(base_indices_row)
            tensor_ch2_qualities.append(quality_scores_row)
            reads_added_to_tensor += 1

        # Pad remaining read rows for tensor
        for _ in range(TENSOR_MAX_READ_ROWS - reads_added_to_tensor):
            tensor_ch1_bases.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            tensor_ch2_qualities.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)

        variant_key_str = f"{v_pos}_{v_type}_{v_ref_defined}_{v_alt_defined}"
        final_variant_output[variant_key_str] = {
            "pileup_reads_data": pileup_reads_data_for_view,  # For console view
            "tensor_channel1_bases": tensor_ch1_bases,  # New tensor output
            "tensor_channel2_qualities": tensor_ch2_qualities,  # New tensor output
            "alt_allele_count": af_alt_count,
            "ref_allele_count_at_locus": af_ref_count,
            "other_allele_count_at_locus": af_other_count,
            "coverage_at_locus": af_locus_coverage,
            "alt_allele_frequency": round(alt_freq, 4)
        }
    return node_id, final_variant_output


# ─────────────────────────────────────────────────────────────────────────────
# Pileup Viewing Function (Integrated) - Remains largely the same, uses "pileup_reads_data"
# ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if not node_data_for_display_view or not isinstance(node_data_for_display_view, dict):
        print(f"ℹ️ No valid pileup data to display for node {node_id_str_for_display}.", file=sys.stderr)
        return

    node_length = node_data_for_display_view.get("node_length")
    variants_dict = node_data_for_display_view.get("variants", {})

    print(
        f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {node_length if node_length is not None else 'N/A'}) ===")

    if not variants_dict:
        print(f"ℹ️ No variants found or pileups generated for this node (or all filtered by AF).")
        return

    variants_displayed_count = 0
    # Sort by position then type
    sorted_variant_keys = sorted(variants_dict.keys(), key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))

    # This window_size is for display marker logic, should be consistent with TENSOR_WINDOW_SIZE
    # if variant_processing_window_size was used consistently
    display_window_size = TENSOR_WINDOW_SIZE
    half_display_window = display_window_size // 2

    for variant_key in sorted_variant_keys:
        if variants_displayed_count >= max_variants_to_display:
            print(f"\n  ... (and {len(variants_dict) - variants_displayed_count} more variants not shown due to limit)")
            break

        variant_data = variants_dict[variant_key]
        # This key "pileup_reads_data" now points to the data formatted for the view
        pileup_reads_display_data = variant_data.get("pileup_reads_data", [])

        v_pos = int(variant_key.split('_')[0])
        v_type = variant_key.split('_')[1]

        if v_type == 'I':
            window_center_on_node = v_pos + 1
        else:
            window_center_on_node = v_pos
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

        ref_display_parts = []
        marker_line_parts = [' '] * display_window_size
        variant_display_idx_in_window = -1

        if v_type != 'I':
            if current_window_start_on_node <= v_pos < current_window_start_on_node + display_window_size:
                variant_display_idx_in_window = v_pos - current_window_start_on_node
        elif v_type == 'I':
            if current_window_start_on_node <= v_pos < current_window_start_on_node + display_window_size:
                variant_display_idx_in_window = v_pos - current_window_start_on_node

        for i in range(display_window_size):
            actual_node_pos_in_window = current_window_start_on_node + i
            if 0 <= actual_node_pos_in_window < len(full_node_sequence):
                ref_display_parts.append(full_node_sequence[actual_node_pos_in_window])
            else:
                ref_display_parts.append(" ")

            if i == variant_display_idx_in_window:
                if v_type == 'I':
                    marker_line_parts[i] = "I"
                    if i + 1 < display_window_size: marker_line_parts[i + 1] = "^"
                else:
                    marker_line_parts[i] = "^"

        print(f"  Node Ref: {''.join(ref_display_parts)}")
        print(f"  Marker  : {''.join(marker_line_parts)}")

        if not pileup_reads_display_data:
            print("  (No reads data in pileup for this variant's window)")
        else:
            displayed_reads_count = 0
            for i, read_info in enumerate(pileup_reads_display_data):
                if displayed_reads_count >= max_reads_to_display_per_variant:
                    print(
                        f"  ... (and {len(pileup_reads_display_data) - max_reads_to_display_per_variant} more reads for this variant's pileup window)")
                    break

                # For display, "bases" are numerical indices from pileup_reads_data_for_view
                base_indices = read_info["bases"]
                read_offset = read_info["offset"]
                read_strand = read_info["strand"]
                read_cigar = read_info.get("cigar", "N/A")

                pileup_row_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in base_indices])
                print(
                    f"  Read {i + 1:3d}: {pileup_row_str}  (Offset: {read_offset}, Strand: {read_strand}, CIGAR: {read_cigar})")
                displayed_reads_count += 1
        variants_displayed_count += 1
    print("\n")


# ─────────────────────────────────────────────────────────────────────────────
# Main function - No major changes needed here other than awareness of new JSON fields
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered pileups for a single specified node and optionally view them.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="JSON output file path for pileups")
    parser.add_argument("--node_id", type=int, required=True, help="The specific node ID to process.")
    parser.add_argument("--gfa", help="GFA graph file path (required if node sequence cache is not used/built).")
    parser.add_argument("--load-cache", help="Load node sequence from this JSON cache file.")
    parser.add_argument("--save-cache", help="Save node sequence to this JSON cache file (used if --gfa is provided).")

    parser.add_argument(
        "--view",
        nargs='?', const=-1, default=None, type=int, metavar='N',
        help="Print generated pileups to console. Optionally specify N to view the first N variants. If no N, all variants are shown."
    )
    parser.add_argument(
        "--max_view_reads", type=int, default=20,  # This is for console view, TENSOR_MAX_READ_ROWS is for JSON tensor
        help="Maximum number of reads to display per pileup matrix in console view (if --view is used)."
    )
    parser.add_argument(
        "--min_af", type=float, default=0.1,
        help="Minimum allele frequency threshold for a variant to be included in the output and pileup generation."
    )

    args = parser.parse_args()

    if not os.path.isfile(args.dat): sys.exit(f"❌ Error: DAT file not found: {args.dat}")
    if not os.path.isfile(args.idx): sys.exit(f"❌ Error: Index file not found: {args.idx}")
    if not args.load_cache and not args.gfa: sys.exit("❌ Error: Must provide --gfa or --load-cache.")
    if args.load_cache and not os.path.isfile(args.load_cache): sys.exit(
        f"❌ Error: Cache file not found: {args.load_cache}")
    if args.gfa and not os.path.isfile(args.gfa): sys.exit(f"❌ Error: GFA file not found: {args.gfa}")
    if args.gfa and args.load_cache: print("🔹 Info: Both --gfa and --load-cache provided. Cache preferred.")
    if args.min_af < 0.0 or args.min_af > 1.0: sys.exit("❌ Error: --min_af must be between 0.0 and 1.0.")

    target_node_id = args.node_id
    print(f"🔹 Processing single target node ID: {target_node_id}")

    start_time = time.time()
    node_dat_info = parse_idx_file_for_single_node(args.idx, target_node_id)
    if not node_dat_info: sys.exit(f"❌ Error: Failed to get index info for node {target_node_id}.")
    dat_offset, n_records = node_dat_info
    print(f"✔ Index parsing for node {target_node_id} took {time.time() - start_time:.2f}s.")

    node_sequence = None
    # ... (node sequence loading logic remains the same) ...
    if args.load_cache and os.path.isfile(args.load_cache):
        start_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                loaded_cache = json.load(cf)
                node_sequence = loaded_cache.get(str(target_node_id))  # Handle if node not in cache
            if node_sequence:
                print(f"✔ Loaded sequence for node {target_node_id} from cache in {time.time() - start_time:.2f}s.")
            else:
                print(f"⚠️ Warning: Node {target_node_id} not in cache {args.load_cache}.")
                if not args.gfa: sys.exit(f"❌ Error: Node {target_node_id} not in cache and no GFA. Exiting.")
        except Exception as e:
            print(f"❌ Error loading cache {args.load_cache}: {e}", file=sys.stderr)
            if not args.gfa: sys.exit(1)
            node_sequence = None  # Ensure it's None if loading failed

    if node_sequence is None and args.gfa:  # Check if still None after cache attempt
        start_time = time.time()
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
                    json.dump(existing_cache_data, wcf)
                print(f"✔ Saved sequence for node {target_node_id} to cache.")
            except Exception as e:
                print(f"❌ Error saving to cache {args.save_cache}: {e}", file=sys.stderr)
        elif node_sequence:
            print(f"✔ Sequence loading for node {target_node_id} from GFA took {time.time() - start_time:.2f}s.")

    if not node_sequence: sys.exit(f"❌ Error: Failed to obtain sequence for node {target_node_id}. Exiting.")

    node_len_for_output = len(node_sequence)
    if not (isinstance(node_len_for_output, int) and node_len_for_output >= 0):
        print(f"❌ Error: Invalid node_len_for_output ({node_len_for_output}) for node {target_node_id}. Exiting.",
              file=sys.stderr)
        sys.exit(1)

    task = (target_node_id, dat_offset, n_records, node_sequence, args.min_af)
    print(f"🔹 Prepared task for node {target_node_id} with min AF threshold: {args.min_af}.")
    output_data_for_json = {}
    print(f"🔹 Processing node {target_node_id} using 1 worker...")
    start_proc_time = time.time()

    try:
        with ProcessPoolExecutor(max_workers=1, initializer=init_worker, initargs=(args.dat,)) as executor:
            future = executor.submit(process_single_node_for_pileup, task)
            processed_node_id, variants_dict_from_worker = future.result()
            variants_dict_to_store = variants_dict_from_worker if isinstance(variants_dict_from_worker, dict) else {}
            output_data_for_json[str(processed_node_id)] = {
                "node_length": node_len_for_output,
                "variants": variants_dict_to_store
            }
            if variants_dict_from_worker is None:
                print(f"⚠️ Warning: Processing for node {processed_node_id} returned None. Stored empty variant dict.",
                      file=sys.stderr)
    except Exception as pool_exc:
        sys.exit(f"\n❌ An error occurred during processing node {target_node_id}: {pool_exc}")

    total_elapsed_time = time.time() - start_proc_time
    print(f"✔ Node {target_node_id} processing finished in {total_elapsed_time:.2f}s.")
    node_result_data_for_output_and_view = output_data_for_json.get(str(target_node_id))
    has_actual_variants = bool(node_result_data_for_output_and_view and
                               isinstance(node_result_data_for_output_and_view.get("variants"), dict) and
                               node_result_data_for_output_and_view.get("variants"))
    should_save_main_output_file = has_actual_variants or args.save_cache

    if should_save_main_output_file:
        if node_result_data_for_output_and_view:
            print(f"🔹 Writing pileup results for node {target_node_id} to JSON output: {args.output}")
            start_write_time = time.time()
            try:
                with open(args.output, 'w') as out_f:
                    json.dump(output_data_for_json, out_f, indent=2)
                print(
                    f"✔ Output written in {time.time() - start_write_time:.2f}s. ✅ Pileup JSON saved to {args.output}")
                if args.view is not None:
                    if has_actual_variants:
                        max_v_show = float('inf') if args.view == -1 else (
                            args.view if args.view >= 0 else float('inf'))
                        if args.view < -1: print("⚠️ Warning: Invalid number for --view. Showing all.", file=sys.stderr)
                        view_msg = f"all variants (max {args.max_view_reads} reads/variant)..."
                        if max_v_show != float(
                            'inf'): view_msg = f"first {int(max_v_show)} variants (max {args.max_view_reads} reads/variant)..."
                        print(f"🔹 Displaying pileups for node {target_node_id}: {view_msg}")
                        display_pileup_data(node_result_data_for_output_and_view, str(target_node_id),
                                            node_sequence, args.max_view_reads, max_v_show)
                    else:
                        print(
                            f"\nℹ️ --view specified for node {target_node_id}, but no variants met AF or found to display.")
                        if node_result_data_for_output_and_view:
                            print(f"   (Node Length: {node_result_data_for_output_and_view.get('node_length', 'N/A')})")
                        else:
                            print(f"   (No data available for node {target_node_id} to display length.)")
            except Exception as e:
                sys.exit(f"❌ Error writing/viewing output: {e}")
    else:
        print(
            f"ℹ️ No variants met AF for node {target_node_id} (or none found) and --save-cache not specified. Main output file '{args.output}' will not be created/overwritten.")
        # ... (logic for skipping output file creation/overwrite remains the same)

    print("✅ Script finished.")


if __name__ == '__main__':
    main()
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
from concurrent.futures import ProcessPoolExecutor, as_completed
import torch

# ─────────────────────────────────────────────────────────────────────────────
# Constants
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size

# Base to Index Mapping
BASE_TO_INDEX = {
    'A': 2, 'C': 3, 'G': 5, 'T': 7,  # Standard bases
    'N': 1,  # Unknown or ambiguous base
    '*': 9,  # Deletion character from CIGAR or gap
    '_PADDING_': 0  # Representing padding, mapped to index 0
}
PADDING_BASE_INDEX = 0

# CIGAR Operation to Index Mapping
CIGAR_OP_TO_INDEX = {
    'M': 1, 'I': 2, 'D': 3, 'N': 4, 'S': 5, 'H': 6, 'P': 7, '=': 8, 'X': 9,
    '_PADDING_': 0
}
CIGAR_PADDING_INDEX = 0

# Index to Base Mapping for console visualization
INDEX_TO_BASE_FOR_VIEW = {
    2: 'A', 3: 'C', 5: 'G', 7: 'T',
    1: 'N',
    9: '*',
    0: ' '  # Padding is represented by '0'
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


# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
# ─────────────────────────────────────────────────────────────────────────────

def calculate_window_start(variant_pos, window_size):
    """
    Calculates the window's start position to place a variant at the center.
    The returned start_pos can be negative, indicating left-side padding is needed.
    """
    center_index = window_size // 2
    return variant_pos - center_index


def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]


def load_full_idx_data(idx_path):
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
            if num_nodes_in_idx == 0: return idx_data_map

            processed_entries = 0
            for i in range(num_nodes_in_idx):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                    sys.stderr.write(f"Error: Index file ended prematurely at record {i + 1}.\n")
                    break
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
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


def load_multiple_node_sequences_from_gfa(gfa_path, target_node_ids_set):
    node_sequences = {}
    if not target_node_ids_set: return node_sequences
    nodes_to_find_int = {int(nid) for nid in target_node_ids_set}

    try:
        with open(gfa_path, 'r') as f:
            print(f"Reading GFA file for {len(nodes_to_find_int)} target nodes: {gfa_path}")
            for line in f:
                if line.startswith('S\t'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        try:
                            nid_int_from_gfa = int(parts[1])
                            if nid_int_from_gfa in nodes_to_find_int:
                                node_sequences[str(nid_int_from_gfa)] = parts[2]
                                if len(node_sequences) == len(nodes_to_find_int):
                                    break
                        except ValueError:
                            continue
            missing_nodes = len(nodes_to_find_int) - len(node_sequences)
            if missing_nodes > 0:
                print(f"Warning: Could not find GFA sequences for {missing_nodes} node(s).")
    except FileNotFoundError:
        sys.stderr.write(f"Error: GFA file not found at {gfa_path}\n")
        return {}
    except Exception as e:
        sys.stderr.write(f"Error reading GFA file {gfa_path}: {e}\n")
    return node_sequences


def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*': return []
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
                if expected_var_type == 'I': return "OTHER_FOR_INDEL", None
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
            deleted_sequence_from_ref = node_sequence[
                                        node_pos: node_pos + length].upper() if node_pos + length <= node_seq_len else ""
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

        if node_pos >= window_start_node + window_size and read_pos > 0: break
        if read_pos >= read_seq_len: break
    return window_chars


def get_read_tensor_rows_in_window(segment_cigar_ops, segment_offset_on_node,
                                   segment_read_sequence, segment_quality_values,
                                   window_start_node, tensor_win_size, node_len):
    bases = [PADDING_BASE_INDEX] * tensor_win_size
    quals = [DEFAULT_QUALITY_PADDING] * tensor_win_size
    cigar_ops_indices = [CIGAR_PADDING_INDEX] * tensor_win_size

    node_pos, read_pos = segment_offset_on_node, 0
    read_seq_len = len(segment_read_sequence)
    qual_len = len(segment_quality_values)

    for L, op in segment_cigar_ops:
        op_idx = CIGAR_OP_TO_INDEX.get(op, CIGAR_PADDING_INDEX)
        if node_pos >= window_start_node + tensor_win_size and read_pos > 0: break

        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = node_pos + i, read_pos + i
                win_idx = n_aln - window_start_node
                if 0 <= win_idx < tensor_win_size:
                    cigar_ops_indices[win_idx] = op_idx
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
                    bases[win_idx] = BASE_TO_INDEX['*']
                    quals[win_idx] = DEFAULT_QUALITY_PADDING
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L

        if read_pos >= read_seq_len: break
    return bases, quals, cigar_ops_indices


# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function
# ─────────────────────────────────────────────────────────────────────────────

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
    except FileNotFoundError:
        sys.stderr.write(f"Error [Worker {os.getpid()}]: DAT file not found.\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()}] opening DAT file: {e}\n")
        sys.exit(1)


def process_single_node_for_pileup(task_args):
    node_id, dat_file_offset, n_records, node_sequence, min_af_threshold, min_variants_threshold, min_allele_bq_threshold, variant_type_to_process = task_args
    global worker_dat_file, worker_base_output_dir

    tensor_files_generated_for_node = 0
    if worker_dat_file is None or worker_base_output_dir is None:
        return node_id, None, tensor_files_generated_for_node
    if not node_sequence:
        return node_id, {}, tensor_files_generated_for_node

    node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
    os.makedirs(node_specific_output_dir, exist_ok=True)

    node_len = len(node_sequence)
    aligned_read_segments = []
    try:
        worker_dat_file.seek(dat_file_offset + 10)
        for _ in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE: break

            off_from_file, raw_seq, raw_qual, raw_cigar, mapq_val, strand_byte = RECORD_STRUCT.unpack(data)
            if mapq_val < 10: continue

            try:
                seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                qual_values = list(raw_qual.rstrip(b'\0'))
                cigar_str_original = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue

            if not seq or len(seq) != len(qual_values): continue

            original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
            if not original_decoded_cigar_ops and cigar_str_original != '*': continue

            current_read_sequence = seq
            current_quality_values = qual_values
            current_decoded_cigar_ops = original_decoded_cigar_ops
            current_offset_on_node = off_from_file

            if strand_char == '-':
                current_read_sequence = reverse_complement(seq)
                current_quality_values = qual_values[::-1]
                current_decoded_cigar_ops = [op for op in
                                             reversed(original_decoded_cigar_ops)] if original_decoded_cigar_ops else []

                # Applying user-confirmed offset logic for reverse strand
                alignment_span_on_node = len(current_read_sequence)
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

    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    variant_headers_for_summary = []
    view_oriented_variant_data = {}

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
        if variant_type_to_process == 'snp' and v_type != 'X':
            continue
        if variant_type_to_process == 'indel' and v_type not in ('I', 'D'):
            continue

        alt_allele_count, ref_allele_count, other_allele_count, locus_coverage = 0, 0, 0, 0
        alt_allele_base_qualities = []

        expected_ref_for_af = v_ref_from_cigar
        expected_alt_for_af = v_alt_from_cigar
        ref_allele_for_indel_context = None

        if v_type == 'D':
            expected_alt_for_af = "*"
            if 0 <= v_pos < node_len: expected_ref_for_af = node_sequence[v_pos]
            ref_allele_for_indel_context = v_ref_from_cigar
        elif v_type == 'I':
            if 0 <= v_pos < node_len:
                expected_ref_for_af = node_sequence[v_pos]
            else:
                expected_ref_for_af = "*"
            ref_allele_for_indel_context = expected_ref_for_af

        for seg in aligned_read_segments:
            allele_observed, bq = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["processed_quality_values"], seg["cigar_ops"],
                v_pos, node_sequence, v_type, ref_allele_for_indel_context)

            if allele_observed is not None:
                locus_coverage += 1
                if allele_observed == expected_alt_for_af:
                    alt_allele_count += 1
                    if bq is not None: alt_allele_base_qualities.append(bq)
                elif allele_observed == expected_ref_for_af or (
                        v_type in ('I', 'D') and allele_observed == "REF_STATE_FOR_INDEL"):
                    ref_allele_count += 1
                else:
                    other_allele_count += 1

        # This filter applies to all variant types
        if alt_allele_count < min_variants_threshold:
            continue

        # These filters only apply to SNPs
        if v_type == 'X':
            current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
            if current_alt_freq < min_af_threshold:
                continue
            mean_alt_bq = sum(alt_allele_base_qualities) / len(
                alt_allele_base_qualities) if alt_allele_base_qualities else 0.0
            if mean_alt_bq < min_allele_bq_threshold:
                continue

        # For indels, we need to calculate these values for the summary, but we don't filter on them
        current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
        mean_alt_bq = sum(alt_allele_base_qualities) / len(
            alt_allele_base_qualities) if alt_allele_base_qualities else 0.0

        variant_key_string = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"
        window_center_pos = v_pos + 1 if v_type == 'I' else v_pos
        window_start_pos = calculate_window_start(window_center_pos, TENSOR_WINDOW_SIZE)

        pileup_data_for_view_json = []
        for read_segment_idx, seg_data in enumerate(aligned_read_segments):
            if read_segment_idx >= TENSOR_MAX_READ_ROWS + 50: break
            row_chars_for_view = get_read_representation_in_window_for_view(
                seg_data["cigar_ops"], seg_data["offset_on_node"], seg_data["read_sequence"],
                window_start_pos, TENSOR_WINDOW_SIZE, node_len)
            if any(char != ' ' for char in row_chars_for_view):
                bases_for_view = [
                    (PADDING_BASE_INDEX if char == ' ' else BASE_TO_INDEX.get(char.upper(), BASE_TO_INDEX['N'])) for
                    char in row_chars_for_view]
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

        ch1_list, ch2_list, ch3_list, ch4_list, ch5_list = [], [], [], [], []

        ref_base_indices_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i, node_pos_in_window in enumerate(range(window_start_pos, window_start_pos + TENSOR_WINDOW_SIZE)):
            if 0 <= node_pos_in_window < node_len:
                ref_base_indices_row[i] = BASE_TO_INDEX.get(node_sequence[node_pos_in_window].upper(),
                                                            BASE_TO_INDEX['N'])

        ch1_list.append(ref_base_indices_row)
        ch2_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
        ch3_list.append([MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE)
        ch4_list.append([DEFAULT_MAPPING_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
        ch5_list.append([CIGAR_PADDING_INDEX] * TENSOR_WINDOW_SIZE)

        reads_added = 0
        for seg_data in aligned_read_segments:
            if reads_added >= TENSOR_MAX_READ_ROWS: break

            base_idx_row, quality_score_row, cigar_op_row = get_read_tensor_rows_in_window(
                seg_data["cigar_ops"], seg_data["offset_on_node"],
                seg_data["read_sequence"], seg_data["processed_quality_values"],
                window_start_pos, TENSOR_WINDOW_SIZE, node_len)

            if any(b != PADDING_BASE_INDEX for b in base_idx_row):
                ch1_list.append(base_idx_row)
                ch2_list.append(quality_score_row)
                mismatch_flags_row = [
                    MISMATCH_COMPARISON_PADDING_VALUE if b == PADDING_BASE_INDEX or r == PADDING_BASE_INDEX else (
                        0 if b == r else 1) for b, r in zip(base_idx_row, ref_base_indices_row)]
                ch3_list.append(mismatch_flags_row)
                mapq = max(0, min(int(seg_data["mapping_quality"]), 127))
                ch4_list.append([mapq] * TENSOR_WINDOW_SIZE)
                ch5_list.append(cigar_op_row)
                reads_added += 1

        for _ in range(TENSOR_MAX_READ_ROWS - reads_added):
            ch1_list.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            ch2_list.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            ch3_list.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)
            ch4_list.append([DEFAULT_MAPPING_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            ch5_list.append([CIGAR_PADDING_INDEX] * TENSOR_WINDOW_SIZE)

        try:
            tensor_chw = torch.tensor([ch1_list, ch2_list, ch3_list, ch4_list, ch5_list], dtype=torch.int8)
            numpy_array_to_save = tensor_chw.numpy()
            tensor_filename_npy = f"{variant_key_string}.npy"
            tensor_filepath_npy = os.path.join(node_specific_output_dir, tensor_filename_npy)
            np.save(tensor_filepath_npy, numpy_array_to_save)

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
        summary_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        with open(summary_path, 'w') as f:
            json.dump({"node_id": node_id, "node_length": node_len,
                       "variants_passing_af_filter": variant_headers_for_summary}, f, indent=2)

    return node_id, view_oriented_variant_data, tensor_files_generated_for_node


# ─────────────────────────────────────────────────────────────────────────────
# Pileup Viewing Function
# ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if max_variants_to_display == 0: return
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
            print(f"  Read {j + 1:3d}: {bases_str} (CIGAR:{read_entry['cigar']})")

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

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--node_id", type=int, help="Specific node ID to process")
    group.add_argument("--node_id_file", help="File with node IDs to process")

    parser.add_argument("--gfa", help="GFA graph file (required if no cache)")
    parser.add_argument("--load-cache", help="Load node sequences from JSON cache")
    parser.add_argument("--save-cache", help="Save node sequences to JSON cache")

    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Number of worker processes")
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

    if not all([os.path.isfile(args.dat), os.path.isfile(args.idx)]):
        sys.exit("Error: DAT or IDX file not found.")
    if not args.gfa and not args.load_cache:
        sys.exit("Error: Must provide --gfa or --load-cache.")
    os.makedirs(args.output, exist_ok=True)

    target_node_ids = set()
    if args.node_id:
        target_node_ids.add(args.node_id)
    else:
        with open(args.node_id_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    target_node_ids.add(int(line.strip()))

    node_sequences = {}
    if args.load_cache and os.path.isfile(args.load_cache):
        with open(args.load_cache, 'r') as f:
            node_sequences = json.load(f)

    nodes_to_fetch = {nid for nid in target_node_ids if str(nid) not in node_sequences}
    if nodes_to_fetch and args.gfa:
        fetched_sequences = load_multiple_node_sequences_from_gfa(args.gfa, nodes_to_fetch)
        for nid, seq in fetched_sequences.items():
            node_sequences[str(nid)] = seq.upper()

    if args.save_cache:
        with open(args.save_cache, 'w') as f:
            json.dump(node_sequences, f, indent=2)

    full_idx_data = load_full_idx_data(args.idx)
    if not full_idx_data: sys.exit("Failed to load index data.")

    tasks = []
    for node_id in target_node_ids:
        if str(node_id) in node_sequences and node_id in full_idx_data:
            offset, n_records = full_idx_data[node_id]
            tasks.append((node_id, offset, n_records, node_sequences[str(node_id)], args.min_af, args.min_variants,
                          args.min_allele_bq, args.variant_type))

    if not tasks:
        sys.exit("No valid tasks to run after checking for sequences and index entries.")

    print(f"\nSubmitting {len(tasks)} tasks to {args.num_workers} workers...")
    total_tensors = 0
    nodes_processed_since_last_report = 0
    tensors_since_last_report = 0
    batch_start_time = time.time()

    with ProcessPoolExecutor(max_workers=args.num_workers, initializer=init_worker,
                             initargs=(args.dat, args.output)) as executor:
        future_to_node = {executor.submit(process_single_node_for_pileup, task): task[0] for task in tasks}

        for i, future in enumerate(as_completed(future_to_node)):
            node_id = future_to_node[future]
            nodes_processed_since_last_report += 1
            try:
                _, view_data, tensor_count = future.result()
                total_tensors += tensor_count
                tensors_since_last_report += tensor_count

                if args.view is not None and view_data:
                    display_pileup_data(view_data, str(node_id), node_sequences[str(node_id)], args.max_view_reads,
                                        args.view if args.view != -1 else float('inf'))
            except Exception as e:
                print(f"Error processing node {node_id}: {e}", file=sys.stderr)

            if nodes_processed_since_last_report >= 1000 or (i + 1) == len(tasks):
                elapsed_time = time.time() - batch_start_time
                rate = nodes_processed_since_last_report / elapsed_time if elapsed_time > 0 else 0
                print(
                    f"  Processed batch of {nodes_processed_since_last_report} nodes (total: {i + 1}/{len(tasks)}) in {elapsed_time:.2f}s "
                    f"({rate:.1f} nodes/sec). Tensors in batch: {tensors_since_last_report}. Total tensors: {total_tensors}.")
                nodes_processed_since_last_report = 0
                tensors_since_last_report = 0
                batch_start_time = time.time()

    print(f"\nProcessing complete. Total tensors generated: {total_tensors}.")


if __name__ == '__main__':
    main()
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

# New Base to Index Mapping as per your request
BASE_TO_INDEX = {
    'A': 2, 'C': 3, 'G': 5, 'T': 7,  # Standard bases
    'N': 1,  # Unknown or ambiguous base
    '*': 9,  # Deletion character from CIGAR or gap
    '_PADDING_': 0  # Representing padding, mapped to index 0
}
PADDING_BASE_INDEX = 0  # Explicitly set padding index to 0

# Updated Index to Base Mapping for console visualization (--view)
INDEX_TO_BASE_FOR_VIEW = {
    2: 'A', 3: 'C', 5: 'G', 7: 'T',
    1: 'N',
    9: '*',
    0: ' '  # Padding index 0 will be visualized as a space
}

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200  # Max reads per tensor (excluding reference row)
DEFAULT_QUALITY_PADDING = 0
DEFAULT_MAPPING_QUALITY_PADDING = -1  # Padding for mapping quality channel (ref row / empty read rows)
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1  # Padding for mismatch channel

# Globals for worker process state
worker_dat_file = None
worker_base_output_dir = None


# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
# ─────────────────────────────────────────────────────────────────────────────

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
                sys.stderr.write(f"Error: Index file {idx_path} is too small (size: {file_size} bytes).\n")
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
                    sys.stderr.write(
                        f"Error: Index file ended prematurely while reading record {i + 1}/{num_nodes_in_idx}. Loaded {processed_entries} entries.\n")
                    break
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                idx_data_map[node_id_from_idx] = (offset, n_records)
                processed_entries += 1
                if processed_entries > 0 and processed_entries % 5_000_000 == 0:
                    print(f"    Loaded {processed_entries}/{num_nodes_in_idx} index entries...")
            if processed_entries != num_nodes_in_idx and len(
                    record_bytes) == 22:
                sys.stderr.write(
                    f"Warning: Index header reported {num_nodes_in_idx} entries, but {processed_entries} were processed from file content.\n")
            print(f"Successfully loaded {len(idx_data_map)} distinct node entries from index file {idx_path}.")
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
            print(f"Reading GFA file to find sequences for {len(nodes_to_find_int)} target nodes: {gfa_path}")
            line_counter = 0
            found_count_gfa = 0
            for line in f:
                line_counter += 1
                if line_counter % 10_000_000 == 0:
                    print(
                        f"  Checked {line_counter:,} lines in GFA file... {len(nodes_to_find_int) - found_count_gfa} nodes remaining to find.")
                if not line.startswith('S\t'): continue
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                try:
                    nid_int_from_gfa = int(parts[1])
                except ValueError:
                    continue

                if nid_int_from_gfa in nodes_to_find_int:
                    node_sequences[str(nid_int_from_gfa)] = parts[2]
                    found_count_gfa += 1
                    if found_count_gfa == len(target_node_ids_set):
                        print(
                            f"Found all {len(target_node_ids_set)} requested node sequences in GFA after checking {line_counter:,} lines.")
                        break

            if found_count_gfa < len(target_node_ids_set):
                print(
                    f"Finished GFA scan ({line_counter:,} lines). Found {found_count_gfa}/{len(target_node_ids_set)} sequences.")
                missing_nodes = [nid for nid in target_node_ids_set if str(nid) not in node_sequences]
                if missing_nodes:
                    print(
                        f"Warning: Could not find GFA sequences for {len(missing_nodes)} node ID(s). Examples: {missing_nodes[:5]}")
            elif found_count_gfa == len(target_node_ids_set):
                print(f"Successfully found all {found_count_gfa} sequences in GFA.")

    except FileNotFoundError:
        sys.stderr.write(f"Error: GFA file not found at {gfa_path}\n")
        return {}
    except Exception as e:
        sys.stderr.write(f"Error reading GFA file {gfa_path}: {e}\n")
        return node_sequences
    return node_sequences


def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*': return []
    ops = []
    try:
        for length_str, op_char in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string):
            ops.append((int(length_str), op_char))
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse CIGAR string '{cigar_string}': {e}\n")
        return []
    return ops


def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_allele_for_indel=None):
    current_node_pos = read_offset_on_node
    current_read_pos = 0
    for length, op in read_cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type in ('I', 'D'): return "REF_STATE_FOR_INDEL"
                offset_in_block = target_node_pos - current_node_pos
                if current_read_pos + offset_in_block < len(read_sequence):
                    return read_sequence[current_read_pos + offset_in_block].upper()
                return None
            current_node_pos += length
            current_read_pos += length
        elif op == 'I':
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                if current_read_pos + length <= len(read_sequence):
                    return read_sequence[current_read_pos: current_read_pos + length].upper()
                return None
            current_read_pos += length
        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I': return "OTHER_FOR_INDEL"
                if expected_var_type == 'D':
                    if 0 <= current_node_pos < len(node_sequence) and current_node_pos + length <= len(node_sequence):
                        deleted_seq_in_ref_context = node_sequence[current_node_pos: current_node_pos + length]
                        if deleted_seq_in_ref_context == expected_ref_allele_for_indel:
                            return "*"
                        else:
                            return "OTHER_FOR_INDEL"
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
                if window_start_node <= n_aln < window_start_node + window_size:
                    win_idx = n_aln - window_start_node
                    if r_aln < read_seq_len:
                        window_chars[win_idx] = segment_read_sequence[r_aln].upper()
            node_pos += L
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                n_aln = node_pos + i
                if window_start_node <= n_aln < window_start_node + window_size:
                    window_chars[n_aln - window_start_node] = '*'
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L

        if node_pos >= window_start_node + window_size and op in ('M', '=', 'X', 'D', 'N'):
            break
        if read_pos >= read_seq_len and op in ('M', '=', 'X', 'I', 'S'):
            break
    return window_chars


def get_read_tensor_rows_in_window(segment_cigar_ops, segment_offset_on_node,
                                   segment_read_sequence, segment_quality_str,
                                   window_start_node, tensor_win_size, node_len):
    bases = [PADDING_BASE_INDEX] * tensor_win_size
    quals = [DEFAULT_QUALITY_PADDING] * tensor_win_size
    node_pos, read_pos = segment_offset_on_node, 0
    read_seq_len = len(segment_read_sequence)
    qual_str_len = len(segment_quality_str)

    for L, op in segment_cigar_ops:
        if node_pos >= window_start_node + tensor_win_size and op in ('M', 'D', 'N', '=', 'X'): break

        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = node_pos + i, read_pos + i
                if r_aln >= read_seq_len: break

                if window_start_node <= n_aln < window_start_node + tensor_win_size:
                    win_idx = n_aln - window_start_node
                    base_char = segment_read_sequence[r_aln].upper()
                    bases[win_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                    if r_aln < qual_str_len:
                        try:
                            quals[win_idx] = ord(segment_quality_str[r_aln]) - 33
                        except (TypeError, ValueError):
                            quals[win_idx] = DEFAULT_QUALITY_PADDING
                    else:
                        quals[win_idx] = DEFAULT_QUALITY_PADDING
            node_pos += L
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                n_aln = node_pos + i
                if window_start_node <= n_aln < window_start_node + tensor_win_size:
                    win_idx = n_aln - window_start_node
                    bases[win_idx] = BASE_TO_INDEX['*']
                    quals[win_idx] = DEFAULT_QUALITY_PADDING
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L

        if read_pos >= read_seq_len and op in ('M', '=', 'X', 'I', 'S'): break
    return bases, quals


# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function
# ─────────────────────────────────────────────────────────────────────────────

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker):
    global worker_dat_file, worker_base_output_dir
    try:
        worker_dat_file = open(dat_file_path_for_worker, 'rb')
        worker_base_output_dir = base_output_dir_for_worker
    except FileNotFoundError:
        sys.stderr.write(f"Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}\n")
        sys.exit(1)


def process_single_node_for_pileup(task_args_with_af_thresh):
    node_id, dat_file_offset, n_records, node_sequence, min_af_threshold = task_args_with_af_thresh
    global worker_dat_file, worker_base_output_dir

    tensor_files_generated_for_node = 0
    if worker_dat_file is None or worker_base_output_dir is None:
        sys.stderr.write(f"Error [Worker {os.getpid()} for Node {node_id}]: Worker not initialized properly.\n")
        return node_id, None, tensor_files_generated_for_node
    if not node_sequence:
        sys.stderr.write(f"Info [Worker {os.getpid()} for Node {node_id}]: No sequence provided. Skipping.\n")
        return node_id, {}, tensor_files_generated_for_node

    node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
    try:
        os.makedirs(node_specific_output_dir, exist_ok=True)
    except OSError as e:
        sys.stderr.write(
            f"Error [Worker {os.getpid()} for Node {node_id}]: Could not create directory {node_specific_output_dir}: {e}\n")
        return node_id, None, tensor_files_generated_for_node

    node_len = len(node_sequence)
    view_oriented_variant_data = {}
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
                qual_str = raw_qual.rstrip(b'\0').decode('ascii', 'replace')
                cigar_str_original = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand_char = strand_byte.decode('ascii')
            except UnicodeDecodeError:
                continue

            if not seq or len(seq) != len(qual_str): continue

            original_decoded_cigar_ops = decode_cigar_to_int_ops(cigar_str_original)
            if not original_decoded_cigar_ops and cigar_str_original != '*': continue

            current_read_sequence = seq
            current_quality_str = qual_str
            current_decoded_cigar_ops = list(original_decoded_cigar_ops)
            current_offset_on_node = off_from_file

            if strand_char == '-':
                current_read_sequence = reverse_complement(seq)
                current_quality_str = qual_str[::-1]
                current_decoded_cigar_ops = [op for op in
                                             reversed(original_decoded_cigar_ops)] if original_decoded_cigar_ops else []
                alignment_span_on_node = len(current_read_sequence)
                current_offset_on_node = node_len - alignment_span_on_node - off_from_file
                if current_offset_on_node < 0:
                    continue

            aligned_read_segments.append({
                "offset_on_node": current_offset_on_node,
                "read_sequence": current_read_sequence,
                "processed_quality_str": current_quality_str,
                "cigar_ops": current_decoded_cigar_ops,
                "original_cigar_str": cigar_str_original,
                "strand": strand_char,
                "mapping_quality": mapq_val
            })
    except Exception as e:
        sys.stderr.write(f"Error [Worker {os.getpid()} for Node {node_id}] during DAT record processing: {e}\n")
        return node_id, None, tensor_files_generated_for_node

    if not aligned_read_segments:
        return node_id, {}, tensor_files_generated_for_node

    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    variant_headers_for_summary = []

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
        alt_allele_count, ref_allele_count, other_allele_count, locus_coverage = 0, 0, 0, 0

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
            allele_observed = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["cigar_ops"],
                v_pos, node_sequence, v_type, ref_allele_for_indel_context)

            if allele_observed is not None:
                locus_coverage += 1
                if allele_observed == expected_alt_for_af:
                    alt_allele_count += 1
                elif allele_observed == expected_ref_for_af or \
                        (v_type in ('I', 'D') and allele_observed == "REF_STATE_FOR_INDEL"):
                    ref_allele_count += 1
                else:
                    other_allele_count += 1

        current_alt_freq = alt_allele_count / locus_coverage if locus_coverage > 0 else 0.0
        if current_alt_freq < min_af_threshold: continue

        variant_key_string = f"{v_pos}_{v_type}_{v_ref_from_cigar}_{v_alt_from_cigar}"

        # ───────────────────────────────────────────────────────────────────────
        # REVISION: Adaptive, Centered Window Calculation
        # The window is now perfectly centered on the variant. Its size is
        # determined by the shorter distance to either end of the node's sequence,
        # ensuring the window never goes out of bounds and requires no padding.
        # ───────────────────────────────────────────────────────────────────────
        max_radius = TENSOR_WINDOW_SIZE // 2

        # For insertions, the variant position is the base *before* the insertion.
        center_pos = v_pos + 1 if v_type == 'I' else v_pos

        dist_to_start = center_pos
        dist_to_end = node_len - 1 - center_pos

        # The actual radius is limited by the shortest side, or the max window size.
        actual_radius = min(dist_to_start, dist_to_end, max_radius)

        adaptive_tensor_width = 2 * actual_radius + 1
        window_start_pos = center_pos - actual_radius

        # Prepare data for visualization using the new adaptive window
        pileup_data_for_view_json = []
        for read_segment_idx, seg_data in enumerate(aligned_read_segments):
            if read_segment_idx >= TENSOR_MAX_READ_ROWS + 50: break

            row_chars_for_view = get_read_representation_in_window_for_view(
                seg_data["cigar_ops"], seg_data["offset_on_node"], seg_data["read_sequence"],
                window_start_pos, adaptive_tensor_width, node_len)

            if any(char_in_row != ' ' for char_in_row in row_chars_for_view):
                pileup_data_for_view_json.append({
                    "bases": [BASE_TO_INDEX.get(char.upper(), BASE_TO_INDEX['N']) for char in row_chars_for_view],
                    "offset": seg_data["offset_on_node"],
                    "strand": seg_data["strand"],
                    "cigar": seg_data["original_cigar_str"]
                })

        # Store window parameters along with view data
        view_oriented_variant_data[variant_key_string] = {
            "window_start": window_start_pos,
            "window_width": adaptive_tensor_width,
            "pileup_reads_data": pileup_data_for_view_json[:TENSOR_MAX_READ_ROWS],
            "alt_allele_count": alt_allele_count, "ref_allele_count_at_locus": ref_allele_count,
            "other_allele_count_at_locus": other_allele_count, "coverage_at_locus": locus_coverage,
            "alt_allele_frequency": round(current_alt_freq, 4)
        }

        # Tensor Data Preparation
        ch1, ch2, ch3, ch4 = [], [], [], []

        ref_row = [PADDING_BASE_INDEX] * adaptive_tensor_width
        for i in range(adaptive_tensor_width):
            abs_pos = window_start_pos + i
            ref_row[i] = BASE_TO_INDEX.get(node_sequence[abs_pos].upper(), BASE_TO_INDEX['N'])

        ch1.append(ref_row)
        ch2.append([DEFAULT_QUALITY_PADDING] * adaptive_tensor_width)
        ch3.append([MISMATCH_CHANNEL_REF_ROW_VALUE] * adaptive_tensor_width)
        ch4.append([DEFAULT_MAPPING_QUALITY_PADDING] * adaptive_tensor_width)

        reads_added = 0
        for seg in aligned_read_segments:
            if reads_added >= TENSOR_MAX_READ_ROWS: break
            base_row, qual_row = get_read_tensor_rows_in_window(
                seg["cigar_ops"], seg["offset_on_node"], seg["read_sequence"],
                seg["processed_quality_str"], window_start_pos, adaptive_tensor_width, node_len)

            if any(b != PADDING_BASE_INDEX for b in base_row):
                ch1.append(base_row)
                ch2.append(qual_row)

                mismatch_row = [MISMATCH_COMPARISON_PADDING_VALUE] * adaptive_tensor_width
                for i in range(adaptive_tensor_width):
                    if base_row[i] != PADDING_BASE_INDEX and ref_row[i] != PADDING_BASE_INDEX:
                        mismatch_row[i] = 0 if base_row[i] == ref_row[i] else 1
                ch3.append(mismatch_row)

                mapq = max(0, min(int(seg["mapping_quality"]), 127))
                ch4.append([mapq] * adaptive_tensor_width)
                reads_added += 1

        for _ in range(TENSOR_MAX_READ_ROWS - reads_added):
            ch1.append([PADDING_BASE_INDEX] * adaptive_tensor_width)
            ch2.append([DEFAULT_QUALITY_PADDING] * adaptive_tensor_width)
            ch3.append([MISMATCH_COMPARISON_PADDING_VALUE] * adaptive_tensor_width)
            ch4.append([DEFAULT_MAPPING_QUALITY_PADDING] * adaptive_tensor_width)

        try:
            tensor = torch.tensor([ch1, ch2, ch3, ch4], dtype=torch.int8).permute(1, 2, 0)
            np.save(os.path.join(node_specific_output_dir, f"{variant_key_string}.npy"), tensor.numpy())

            variant_headers_for_summary.append({
                "variant_key": variant_key_string, "tensor_file": f"{variant_key_string}.npy",
                "alt_allele_count": alt_allele_count, "ref_allele_count_at_locus": ref_allele_count,
                "other_allele_count_at_locus": other_allele_count, "coverage_at_locus": locus_coverage,
                "alt_allele_frequency": round(current_alt_freq, 4)
            })
            tensor_files_generated_for_node += 1
        except Exception as e:
            sys.stderr.write(f"Error creating tensor for {variant_key_string}: {e}\n")

    if variant_headers_for_summary:
        summary_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        try:
            with open(summary_path, 'w') as f:
                json.dump({
                    "node_id": node_id, "node_length": node_len,
                    "node_sequence_preview": node_sequence[:100] + ("..." if node_len > 100 else ""),
                    "variants_passing_af_filter": variant_headers_for_summary
                }, f, indent=2)
        except Exception as e:
            sys.stderr.write(f"Error writing summary for Node {node_id}: {e}\n")

    return node_id, view_oriented_variant_data, tensor_files_generated_for_node


# ─────────────────────────────────────────────────────────────────────────────
# Pileup Viewing Function
# ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if max_variants_to_display == 0:
        return
    if not node_data_for_display_view:
        return

    print(f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")

    variants_displayed_count = 0
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))

    for variant_key in sorted_variant_keys:
        if variants_displayed_count >= max_variants_to_display:
            print(f"\n  ... (and {len(node_data_for_display_view) - variants_displayed_count} more variants not shown)")
            break

        variant_data = node_data_for_display_view[variant_key]
        v_pos_str, v_type = variant_key.split('_')[:2]
        v_pos = int(v_pos_str)

        # REVISION: Read window parameters directly from the data
        window_start_pos = variant_data.get('window_start', 0)
        display_window_size = variant_data.get('window_width', 0)
        if display_window_size == 0: continue

        center_pos = v_pos + 1 if v_type == 'I' else v_pos

        print(f"\n--- Variant: {variant_key} (Node Pos: {v_pos}, Type: {v_type}) ---")
        print(
            f"  Display Window (0-based node coords): {window_start_pos}-{window_start_pos + display_window_size - 1}")
        alt_freq = variant_data.get('alt_allele_frequency', 0.0)
        print(f"  Alt Count: {variant_data.get('alt_allele_count', 'N/A')} | "
              f"Ref Count: {variant_data.get('ref_allele_count_at_locus', 'N/A')} | "
              f"Coverage: {variant_data.get('coverage_at_locus', 'N/A')} | "
              f"Alt Freq: {alt_freq:.4f}")

        ref_display_chars = [' '] * display_window_size
        marker_line_chars = [' '] * display_window_size

        # The variant position relative to the start of our new window
        variant_pos_in_window = center_pos - window_start_pos

        for i in range(display_window_size):
            absolute_node_pos = window_start_pos + i
            if 0 <= absolute_node_pos < len(full_node_sequence):
                ref_display_chars[i] = full_node_sequence[absolute_node_pos]

        if 0 <= variant_pos_in_window < display_window_size:
            marker_line_chars[variant_pos_in_window] = "I" if v_type == 'I' else "^"
            if v_type == 'I':
                if variant_pos_in_window + 1 < display_window_size:
                    marker_line_chars[variant_pos_in_window + 1] = "^"

        print(f"  Node Ref: {''.join(ref_display_chars)}")
        print(f"  Marker  : {''.join(marker_line_chars)}")

        pileup_reads_for_variant = variant_data.get("pileup_reads_data", [])
        if not pileup_reads_for_variant:
            print("  (No reads in window for display)")
        else:
            for i, read_entry in enumerate(pileup_reads_for_variant):
                if i >= max_reads_to_display_per_variant:
                    print(f"  ... ({len(pileup_reads_for_variant) - i} more reads not shown)")
                    break
                bases_str_for_read = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in read_entry["bases"]])
                print(
                    f"  Read {i + 1:3d}: {bases_str_for_read}  (Off:{read_entry['offset']},Str:{read_entry['strand']},CIG:{read_entry.get('cigar', 'N/A')})")
        variants_displayed_count += 1
    print()


# ─────────────────────────────────────────────────────────────────────────────
# Main function
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered 4-channel .npy tensors from .dat alignment files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="Base output directory")
    # ... (rest of main function is unchanged)
    input_node_group = parser.add_mutually_exclusive_group(required=True)
    input_node_group.add_argument("--node_id", type=int, help="The specific node ID to process.")
    input_node_group.add_argument("--node_id_file",
                                  help="Path to a text file containing node IDs (one per line) to process.")

    parser.add_argument("--gfa", help="GFA graph file path (required if node sequences are not cached).")
    parser.add_argument("--load-cache", help="Load node sequences from this JSON cache file.")
    parser.add_argument("--save-cache",
                        help="Save/update node sequences to this JSON cache file (used if --gfa is provided).")

    parser.add_argument("--num_workers", type=int, default=None,
                        help="Number of worker processes. Defaults to a heuristic based on os.cpu_count().")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N_VARIANTS',
                        help="Print generated pileups to console for displayed variants. "
                             "Provide no value or -1 to view all variants passing AF for processed nodes. "
                             "Provide an integer N > 0 to view the first N variants per node. "
                             "Provide 0 to disable all stdout pileup views. "
                             "Note: This only affects console view, JSON summaries are still generated.")
    parser.add_argument("--max_view_reads", type=int, default=20,
                        help="Maximum number of reads to display per pileup in console view.")
    parser.add_argument("--min_af", type=float, default=0.1,
                        help="Minimum allele frequency for a variant to be processed for tensor generation and JSON summary.")
    args = parser.parse_args()

    if not os.path.isfile(args.dat): sys.exit(f"Error: DAT file not found: {args.dat}")
    if not os.path.isfile(args.idx): sys.exit(f"Error: Index file not found: {args.idx}")
    if not args.load_cache and not args.gfa:
        sys.exit("Error: Must provide --gfa or --load-cache to obtain node sequences.")
    if args.load_cache and not os.path.isfile(args.load_cache) and os.path.exists(
            args.load_cache):
        sys.exit(f"Error: Cache path '{args.load_cache}' exists but is not a file.")
    if args.gfa and not os.path.isfile(args.gfa):
        sys.exit(f"Error: GFA file not found: {args.gfa}")
    if not (0.0 <= args.min_af <= 1.0):
        sys.exit("Error: --min_af must be between 0.0 and 1.0.")

    effective_num_workers = args.num_workers if args.num_workers and args.num_workers > 0 else (os.cpu_count() or 1)
    print(f"Using {effective_num_workers} worker process(es) for parallel processing.")
    try:
        os.makedirs(args.output, exist_ok=True)
        print(f"Base output directory: {args.output}")
    except OSError as e:
        sys.exit(f"Error: Could not create base output directory {args.output}: {e}")

    target_node_ids_int_set = set()
    if args.node_id_file:
        try:
            with open(args.node_id_file, 'r') as f_nodes:
                for line_num, line in enumerate(f_nodes, 1):
                    line = line.strip()
                    if line and not line.startswith("#"):
                        try:
                            target_node_ids_int_set.add(int(line))
                        except ValueError:
                            sys.stderr.write(
                                f"Warning: Invalid non-integer node ID '{line}' in {args.node_id_file} at line {line_num}. Skipping.\n")
            if not target_node_ids_int_set:
                sys.exit(f"Error: No valid node IDs found or specified in {args.node_id_file}.")
            print(f"Will process {len(target_node_ids_int_set)} unique node ID(s) from file: {args.node_id_file}")
        except FileNotFoundError:
            sys.exit(f"Error: Node ID file not found: {args.node_id_file}")
    elif args.node_id is not None:
        target_node_ids_int_set.add(args.node_id)
        print(f"Will process single target node ID: {args.node_id}")

    if not target_node_ids_int_set:
        sys.exit("Info: No target node IDs specified. Exiting.")

    node_sequences_map_str_keys = {}
    if args.load_cache and os.path.isfile(args.load_cache):
        cache_load_start_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                node_sequences_map_str_keys = json.load(cf)
            print(
                f"Loaded {len(node_sequences_map_str_keys)} sequences from cache '{args.load_cache}' in {time.time() - cache_load_start_time:.2f}s.")
        except Exception as e:
            sys.stderr.write(
                f"Warning: Error loading cache {args.load_cache}: {e}. Proceeding without cached sequences if GFA is provided.\n")
            node_sequences_map_str_keys = {}

    nodes_needing_gfa_fetch = {nid_int for nid_int in target_node_ids_int_set if
                               str(nid_int) not in node_sequences_map_str_keys}
    if nodes_needing_gfa_fetch and args.gfa:
        print(f"{len(nodes_needing_gfa_fetch)} node(s) require sequence fetching from GFA: {args.gfa}")
        gfa_load_start_time = time.time()
        fetched_sequences_map = load_multiple_node_sequences_from_gfa(args.gfa, nodes_needing_gfa_fetch)
        node_sequences_map_str_keys.update(fetched_sequences_map)
        print(
            f"Fetched {len(fetched_sequences_map)} new sequences from GFA in {time.time() - gfa_load_start_time:.2f}s. Total sequences in map: {len(node_sequences_map_str_keys)}.")
    elif nodes_needing_gfa_fetch:
        sys.stderr.write(
            f"Warning: {len(nodes_needing_gfa_fetch)} node(s) need sequences from GFA, but --gfa argument was not provided. These nodes will be skipped.\n")

    idx_data_load_start_time = time.time()
    full_idx_data_map = load_full_idx_data(args.idx)
    if full_idx_data_map is None:
        sys.exit(f"Error: Failed to load index data from {args.idx}. Cannot proceed.")
    print(f"Index data with {len(full_idx_data_map)} entries loaded in {time.time() - idx_data_load_start_time:.2f}s.")

    tasks_for_submission = []
    skipped_nodes_count_pre_submit = 0
    print(f"Preparing tasks for {len(target_node_ids_int_set)} target nodes...")
    task_prep_s_time = time.time()

    for i, node_id_val_int in enumerate(target_node_ids_int_set):
        if (i + 1) % 50000 == 0 and i > 0: print(
            f"  Prepared tasks for {i + 1}/{len(target_node_ids_int_set)} nodes...")

        node_sequence_val = node_sequences_map_str_keys.get(str(node_id_val_int))
        node_dat_info_tuple = full_idx_data_map.get(node_id_val_int)

        if not node_sequence_val or not node_dat_info_tuple:
            skipped_nodes_count_pre_submit += 1
            continue
        tasks_for_submission.append(
            (node_id_val_int, node_dat_info_tuple[0], node_dat_info_tuple[1], node_sequence_val, args.min_af))

    print(f"Task preparation completed in {time.time() - task_prep_s_time:.2f}s.")
    if skipped_nodes_count_pre_submit > 0:
        print(
            f"Warning: Skipped {skipped_nodes_count_pre_submit} nodes before submission due to missing sequence or index data.")
    if not tasks_for_submission:
        sys.exit("Info: No valid tasks to process after filtering. Exiting.")

    overall_processing_start_time = time.time()
    total_tensor_files_generated_all_nodes = 0
    nodes_completed_by_worker = 0
    nodes_with_actual_output = 0

    batch_start_time = time.time()
    nodes_in_batch_count = 0
    tensors_in_batch_count = 0

    print(f"\nSubmitting {len(tasks_for_submission)} node tasks to {effective_num_workers} worker(s)...")
    with ProcessPoolExecutor(max_workers=effective_num_workers, initializer=init_worker,
                             initargs=(args.dat, args.output)) as executor:
        future_to_node_map = {executor.submit(process_single_node_for_pileup, task): task[0] for task in
                              tasks_for_submission}

        for future_idx, completed_future in enumerate(as_completed(future_to_node_map)):
            processed_task_count_total = future_idx + 1
            original_node_id_for_future = future_to_node_map[completed_future]

            try:
                returned_node_id, view_data_dict, tensor_files_count_for_node = completed_future.result()
                nodes_completed_by_worker += 1

                if returned_node_id is None:
                    sys.stderr.write(
                        f"Error: Worker failed for node {original_node_id_for_future} (returned None for node_id).\n")
                else:
                    total_tensor_files_generated_all_nodes += tensor_files_count_for_node
                    tensors_in_batch_count += tensor_files_count_for_node

                    summary_file_path = os.path.join(args.output, str(returned_node_id), "variant_summary.json")
                    if tensor_files_count_for_node > 0 or os.path.exists(summary_file_path):
                        nodes_with_actual_output += 1

                    if args.view is not None and args.view != 0 and view_data_dict:
                        node_seq_for_view = node_sequences_map_str_keys.get(str(returned_node_id))
                        if node_seq_for_view:
                            max_variants_for_this_node = float('inf') if args.view == -1 else args.view
                            display_pileup_data(view_data_dict, str(returned_node_id), node_seq_for_view,
                                                args.max_view_reads, max_variants_for_this_node)
            except Exception as exc:
                nodes_completed_by_worker += 1
                sys.stderr.write(
                    f"Error processing node {original_node_id_for_future} (exception from worker future): {exc}\n")

            nodes_in_batch_count += 1
            if nodes_in_batch_count >= 1000 or processed_task_count_total == len(tasks_for_submission):
                if nodes_in_batch_count > 0:
                    current_batch_duration = time.time() - batch_start_time
                    processing_rate = nodes_in_batch_count / current_batch_duration if current_batch_duration > 0 else 0
                    print(
                        f"  Processed batch of {nodes_in_batch_count} nodes (total completed: {processed_task_count_total}/{len(tasks_for_submission)}) "
                        f"in {current_batch_duration:.2f}s ({processing_rate:.2f} nodes/sec). "
                        f"Generated {tensors_in_batch_count} .npy files in this batch.")
                    batch_start_time = time.time()
                    nodes_in_batch_count = 0
                    tensors_in_batch_count = 0

    print("\n══════════ PROCESSING COMPLETE ══════════")
    if args.save_cache and node_sequences_map_str_keys:
        print(f"\nSaving {len(node_sequences_map_str_keys)} sequences to cache: {args.save_cache}...")
        try:
            with open(args.save_cache, 'w') as wcf:
                json.dump(node_sequences_map_str_keys, wcf, indent=2)
            print(f"Sequences saved to cache.")
        except Exception as e:
            sys.stderr.write(f"Error saving node sequences to cache {args.save_cache}: {e}\n")
    elif args.save_cache:
        print(f"Info: --save-cache specified, but no sequences in memory to save (map is empty).")

    print(f"\nFinal Summary:")
    print(f"  Total unique node IDs targeted: {len(target_node_ids_int_set)}")
    if skipped_nodes_count_pre_submit > 0: print(
        f"  Nodes skipped before submission (missing sequence/index): {skipped_nodes_count_pre_submit}")
    print(f"  Tasks submitted to workers: {len(tasks_for_submission)}")
    print(
        f"  Node tasks completed by workers (includes tasks that may have errored in worker): {nodes_completed_by_worker}/{len(tasks_for_submission)}")
    print(f"  Nodes with output files (summary/tensors) generated: {nodes_with_actual_output}")
    print(f"  Total .npy tensor files generated across all nodes: {total_tensor_files_generated_all_nodes}")
    print(f"Parallel processing phase finished in {time.time() - overall_processing_start_time:.2f} seconds.")


if __name__ == '__main__':
    main()
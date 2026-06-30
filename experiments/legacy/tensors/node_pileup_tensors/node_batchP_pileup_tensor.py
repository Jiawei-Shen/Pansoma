#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np  # Ensure NumPy is imported
from collections import defaultdict
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
import torch

# ─────────────────────────────────────────────────────────────────────────────
# Constants
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 5, ' ': 6, '-': 6}
INDEX_TO_BASE_FOR_VIEW = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N', 5: '*', 6: ' '}
TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
PADDING_BASE_INDEX = BASE_TO_INDEX[' ']
DEFAULT_QUALITY_PADDING = 0
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1

# Globals for worker process state (set by initializer)
worker_dat_file = None
worker_base_output_dir = None


# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions (Assumed to be the same as in your provided script)
# ─────────────────────────────────────────────────────────────────────────────

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]


def load_full_idx_data(idx_path):
    idx_data_map = {}
    print(f"🔹 Loading full index data from {idx_path}...")
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                sys.stderr.write(f"❌ Error: Index file {idx_path} is too small (size: {file_size} bytes).\n")
                return None
            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                sys.stderr.write(f"❌ Error: Could not read number of nodes from {idx_path}.\n")
                return None
            num_nodes_in_idx = struct.unpack('<I', num_nodes_bytes)[0]
            print(f"  Index file reports {num_nodes_in_idx} total node entries. Reading all entries...")
            if num_nodes_in_idx == 0: return idx_data_map
            expected_min_size = 4 + (num_nodes_in_idx * 22)
            if file_size < expected_min_size:
                sys.stderr.write(
                    f"⚠️ Warning: Index file size ({file_size} bytes) is smaller than expected ({expected_min_size} bytes) for {num_nodes_in_idx} records. File may be truncated.\n")

            processed_entries = 0
            for i in range(num_nodes_in_idx):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                    sys.stderr.write(
                        f"❌ Error: Index file ended prematurely while reading record {i + 1}/{num_nodes_in_idx}. Loaded {processed_entries} entries.\n")
                    break
                node_id_from_idx, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                idx_data_map[node_id_from_idx] = (offset, n_records)
                processed_entries += 1
                if processed_entries > 0 and processed_entries % 2000000 == 0:
                    print(f"    Loaded {processed_entries}/{num_nodes_in_idx} index entries...")
            if processed_entries != num_nodes_in_idx and len(record_bytes) == 22:
                sys.stderr.write(
                    f"⚠️ Warning: Read {processed_entries} entries, but index header indicated {num_nodes_in_idx}.\n")
            print(f"✔ Successfully loaded {len(idx_data_map)} distinct node entries from index file {idx_path}.")
        return idx_data_map
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: Index file not found at {idx_path}\n")
        return None
    except Exception as e:
        sys.stderr.write(f"❌ Error parsing full index file {idx_path}: {e}\n")
        return None


def load_multiple_node_sequences_from_gfa(gfa_path, target_node_ids_set):
    node_sequences = {}
    if not target_node_ids_set: return node_sequences
    nodes_to_find = target_node_ids_set.copy()
    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file to find sequences for {len(nodes_to_find)} nodes: {gfa_path}")
            line_counter = 0
            for line in f:
                line_counter += 1
                if line_counter % 10_000_000 == 0:
                    print(
                        f"  Checked {line_counter:,} lines in GFA file... {len(nodes_to_find)} nodes remaining to find.")
                if not line.startswith('S\t'): continue
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                try:
                    nid_int_from_gfa = int(parts[1])
                except ValueError:
                    continue
                if nid_int_from_gfa in nodes_to_find:
                    node_sequences[str(nid_int_from_gfa)] = parts[2]
                    nodes_to_find.remove(nid_int_from_gfa)
                    if not nodes_to_find:
                        print(
                            f"✔ Found all {len(target_node_ids_set)} requested node sequences in GFA after checking {line_counter:,} lines.")
                        break
            found_count = len(node_sequences)
            requested_count = len(target_node_ids_set)
            if not nodes_to_find:
                if found_count != requested_count:
                    sys.stderr.write(
                        f"⚠️ GFA Load: Mismatch - found_count {found_count}, requested {requested_count}, but all marked found.\n")
            else:
                print(f"✔ Finished GFA scan ({line_counter:,} lines). Found {found_count}/{requested_count} sequences.")
                print(
                    f"⚠️ Warning: Could not find GFA sequences for {len(nodes_to_find)} node ID(s). Examples: {list(nodes_to_find)[:5]}")
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: GFA file not found at {gfa_path}\n")
        return {}
    except Exception as e:
        sys.stderr.write(f"❌ Error reading GFA file {gfa_path}: {e}\n")
        return node_sequences
    return node_sequences


def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*': return []
    ops = []
    try:
        for length_str, op_char in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string):
            ops.append((int(length_str), op_char))
    except Exception as e:
        sys.stderr.write(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}\n")
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
                    deleted_seq_in_read_context = node_sequence[current_node_pos: current_node_pos + length]
                    if deleted_seq_in_read_context == expected_ref_allele_for_indel:
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
    node_pos, read_pos = offset_on_node, 0
    node_seq_len, read_seq_len = len(node_sequence), len(read_sequence)
    for length, op in cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            for i in range(length):
                cur_node_p, cur_read_p = node_pos + i, read_pos + i
                if cur_node_p < node_seq_len and cur_read_p < read_seq_len:
                    node_base, read_base = node_sequence[cur_node_p].upper(), read_sequence[cur_read_p].upper()
                    if node_base != read_base and op != '=':  # Mismatches (X or M where bases differ)
                        variants.append((cur_node_p, 'X', read_base, node_base))
                else:
                    break
            node_pos += length;
            read_pos += length
        elif op == 'I':
            ins_seq = read_sequence[read_pos: read_pos + length].upper()
            anchor_pos = node_pos - 1 if node_pos > 0 else 0
            anchor_base = node_sequence[anchor_pos].upper() if 0 <= anchor_pos < node_seq_len else "*"
            variants.append((anchor_pos, 'I', ins_seq, anchor_base))
            read_pos += length
        elif op == 'D':
            del_seq = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= node_seq_len else ""
            if del_seq:
                variants.append((node_pos, 'D', "*", del_seq))
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
    for L, op in segment_cigar_ops:
        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = node_pos + i, read_pos + i
                if window_start_node <= n_aln < window_start_node + window_size:
                    win_idx = n_aln - window_start_node
                    if r_aln < len(segment_read_sequence): window_chars[win_idx] = segment_read_sequence[r_aln].upper()
            node_pos += L;
            read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                n_aln = node_pos + i
                if window_start_node <= n_aln < window_start_node + window_size:
                    window_chars[n_aln - window_start_node] = '*'
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L
        if node_pos >= window_start_node + window_size and op in ('M', '=', 'X', 'D', 'N'): break
    return window_chars


def get_read_tensor_rows_in_window(segment_cigar_ops, segment_offset_on_node,
                                   segment_read_sequence, segment_quality_str,
                                   window_start_node, tensor_win_size, node_len):
    bases = [PADDING_BASE_INDEX] * tensor_win_size
    quals = [DEFAULT_QUALITY_PADDING] * tensor_win_size
    node_pos, read_pos = segment_offset_on_node, 0
    for L, op in segment_cigar_ops:
        if node_pos >= window_start_node + tensor_win_size and op in ('M', 'D', 'N', '=', 'X'): break
        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = node_pos + i, read_pos + i
                if r_aln >= len(segment_read_sequence): break
                if window_start_node <= n_aln < window_start_node + tensor_win_size:
                    win_idx = n_aln - window_start_node
                    base_char = segment_read_sequence[r_aln].upper()
                    bases[win_idx] = BASE_TO_INDEX.get(base_char, BASE_TO_INDEX['N'])
                    if r_aln < len(segment_quality_str):
                        try:
                            quals[win_idx] = ord(segment_quality_str[r_aln]) - 33
                        except:
                            quals[win_idx] = DEFAULT_QUALITY_PADDING  # Keep default on error
                    else:
                        quals[win_idx] = DEFAULT_QUALITY_PADDING
            node_pos += L;
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
        if read_pos >= len(segment_read_sequence) and op in ('M', '=', 'X', 'I', 'S'): break
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
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path_for_worker}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path_for_worker}: {e}\n")
        sys.exit(1)


def process_single_node_for_pileup(task_args_with_af_thresh):
    node_id, dat_file_offset, n_records, node_sequence, min_af_threshold = task_args_with_af_thresh
    global worker_dat_file, worker_base_output_dir

    npy_files_generated_for_node = 0  # Renamed from pth_files_generated_for_node
    if worker_dat_file is None or worker_base_output_dir is None:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Worker not initialized.\n")
        return node_id, None, npy_files_generated_for_node
    if not node_sequence:
        sys.stderr.write(f"ℹ️ [Worker {os.getpid()} for Node {node_id}]: No sequence. Skipping.\n")
        return node_id, {}, npy_files_generated_for_node

    node_specific_output_dir = os.path.join(worker_base_output_dir, str(node_id))
    try:
        os.makedirs(node_specific_output_dir, exist_ok=True)
    except OSError as e:
        sys.stderr.write(
            f"❌ Error [Worker {os.getpid()} for Node {node_id}]: Dir creation {node_specific_output_dir}: {e}\n")
        return node_id, None, npy_files_generated_for_node

    node_len = len(node_sequence)
    view_oriented_variant_data = {}
    aligned_read_segments = []
    try:
        worker_dat_file.seek(dat_file_offset + 10)
        for _ in range(n_records):
            data = worker_dat_file.read(RECORD_SIZE)
            if len(data) < RECORD_SIZE: break
            off, raw_seq, raw_qual, raw_cigar, mapq, strand_b = RECORD_STRUCT.unpack(data)
            if mapq < 10: continue
            try:
                seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                qual = raw_qual.rstrip(b'\0').decode('ascii', 'replace')
                cigar_orig = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand = strand_b.decode('ascii')
            except UnicodeDecodeError:
                continue
            if not seq or len(seq) != len(qual): continue
            cigar_ops_orig = decode_cigar_to_int_ops(cigar_orig)
            if not cigar_ops_orig and cigar_orig != '*': continue

            cur_seq, cur_qual, cur_cigar_ops, cur_offset = seq, qual, list(cigar_ops_orig), off
            if strand == '-':
                cur_seq, cur_qual = reverse_complement(seq), qual[::-1]
                cur_cigar_ops = [op for op in reversed(cigar_ops_orig)] if cigar_ops_orig else []
                read_len_span = len(cur_seq)
                if read_len_span > 0:
                    cur_offset = node_len - read_len_span - off
                    if cur_offset < 0: continue
            aligned_read_segments.append({
                "offset_on_node": cur_offset, "read_sequence": cur_seq,
                "processed_quality_str": cur_qual, "cigar_ops": cur_cigar_ops,
                "original_cigar_str": cigar_orig, "strand": strand})
    except Exception as e:
        sys.stderr.write(f"❌ Error [Worker {os.getpid()} for Node {node_id}] reading DAT: {e}\n")
        return node_id, None, npy_files_generated_for_node

    if not aligned_read_segments: return node_id, {}, npy_files_generated_for_node

    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    variant_headers = []
    half_win = TENSOR_WINDOW_SIZE // 2
    for (v_pos, v_type, v_ref_cig, v_alt_cig), _ in candidate_variants.items():
        alt_c, ref_c, other_c, locus_cov = 0, 0, 0, 0
        indel_ref_check = None
        exp_ref_af, exp_alt_af = v_ref_cig, v_alt_cig  # Default for SNV
        if v_type == 'D':
            exp_alt_af = "*"
            indel_ref_check = v_ref_cig  # The deleted sequence from ref
            if v_pos < node_len: exp_ref_af = node_sequence[v_pos]  # Anchor base
        elif v_type == 'I':
            indel_ref_check = node_sequence[v_pos] if 0 <= v_pos < node_len else None  # Anchor base
            # exp_ref_af is v_ref_cig (anchor base)
            # exp_alt_af is v_alt_cig (inserted sequence)

        for seg in aligned_read_segments:
            allele = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["cigar_ops"],
                v_pos, node_sequence, v_type, indel_ref_check)
            if allele is not None:
                locus_cov += 1
                if allele == exp_alt_af:
                    alt_c += 1
                elif allele == exp_ref_af or (v_type in 'ID' and allele == "REF_STATE_FOR_INDEL"):
                    ref_c += 1
                else:
                    other_c += 1

        alt_freq = alt_c / locus_cov if locus_cov > 0 else 0.0
        if alt_freq < min_af_threshold: continue

        key_str = f"{v_pos}_{v_type}_{v_ref_cig}_{v_alt_cig}"
        win_center = v_pos + 1 if v_type == 'I' else v_pos
        win_start = max(0, win_center - half_win)

        view_reads_data = []
        for s_idx, seg in enumerate(aligned_read_segments):
            if s_idx >= TENSOR_MAX_READ_ROWS * 2: break  # Limit reads for view generation for performance
            row_chars = get_read_representation_in_window_for_view(
                seg["cigar_ops"], seg["offset_on_node"], seg["read_sequence"],
                win_start, TENSOR_WINDOW_SIZE, node_len)
            if any(c != ' ' for c in row_chars):
                view_reads_data.append({
                    "bases": [BASE_TO_INDEX.get(c.upper(), BASE_TO_INDEX['N']) for c in row_chars],
                    "offset": seg["offset_on_node"], "strand": seg["strand"], "cigar": seg["original_cigar_str"]})
        view_oriented_variant_data[key_str] = {
            "pileup_reads_data": view_reads_data[:TENSOR_MAX_READ_ROWS],
            "alt_allele_count": alt_c, "ref_allele_count_at_locus": ref_c,
            "other_allele_count_at_locus": other_c, "coverage_at_locus": locus_cov,
            "alt_allele_frequency": round(alt_freq, 4)}

        ch1, ch2, ch3 = [], [], []
        ref_tensor_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i in range(TENSOR_WINDOW_SIZE):
            abs_p = win_start + i
            if 0 <= abs_p < node_len: ref_tensor_row[i] = BASE_TO_INDEX.get(node_sequence[abs_p].upper(),
                                                                            BASE_TO_INDEX['N'])
        ch1.append(ref_tensor_row)
        ch2.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
        ch3.append([MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE)

        reads_added = 0
        for seg in aligned_read_segments:
            if reads_added >= TENSOR_MAX_READ_ROWS: break
            base_r, qual_r = get_read_tensor_rows_in_window(
                seg["cigar_ops"], seg["offset_on_node"], seg["read_sequence"],
                seg["processed_quality_str"], win_start, TENSOR_WINDOW_SIZE, node_len)
            if any(b != PADDING_BASE_INDEX for b in base_r):
                ch1.append(base_r);
                ch2.append(qual_r)
                mismatch_r = [MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE
                for i in range(TENSOR_WINDOW_SIZE):
                    if base_r[i] == PADDING_BASE_INDEX or ref_tensor_row[i] == PADDING_BASE_INDEX: continue
                    mismatch_r[i] = 0 if base_r[i] == ref_tensor_row[i] else 1
                ch3.append(mismatch_r)
                reads_added += 1
        for _ in range(TENSOR_MAX_READ_ROWS - reads_added):
            ch1.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            ch2.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            ch3.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)
        try:
            # Create PyTorch tensor
            final_tensor = torch.tensor([ch1, ch2, ch3], dtype=torch.int8)

            # Convert to NumPy array
            numpy_array_to_save = final_tensor.numpy()

            # Save as .npy file
            tensor_filename = f"{key_str}.npy"  # Changed extension
            tensor_filepath = os.path.join(node_specific_output_dir, tensor_filename)
            np.save(tensor_filepath, numpy_array_to_save)  # Use numpy.save

            variant_headers.append({
                "variant_key": key_str, "tensor_file": tensor_filename,  # Reflects .npy
                "alt_allele_count": alt_c, "ref_allele_count_at_locus": ref_c,
                "other_allele_count_at_locus": other_c, "coverage_at_locus": locus_cov,
                "alt_allele_frequency": round(alt_freq, 4)})
            npy_files_generated_for_node += 1  # Increment for .npy file
        except Exception as e:
            sys.stderr.write(f"❌ Error [Worker Node {node_id}]: Tensor save as .npy for {key_str}: {e}\n")

    if variant_headers:
        summary_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        try:
            with open(summary_path, 'w') as sjf:
                json.dump({"node_id": node_id, "node_length": node_len,
                           "node_sequence_preview": node_sequence[:100] + ("..." if node_len > 100 else ""),
                           "variants_passing_af_filter": variant_headers}, sjf, indent=2)
        except Exception as e:
            sys.stderr.write(f"❌ Error [Worker Node {node_id}]: Summary JSON write: {e}\n")

    return node_id, view_oriented_variant_data, npy_files_generated_for_node


# ─────────────────────────────────────────────────────────────────────────────
# Pileup Viewing Function
# ─────────────────────────────────────────────────────────────────────────────
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if not node_data_for_display_view:
        print(f"ℹ️ No pileup data to display for node {node_id_str_for_display}.", file=sys.stderr)
        return
    print(f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")
    if not node_data_for_display_view:
        print(f"ℹ️ No variants met AF threshold or found for node {node_id_str_for_display}.")
        return

    variants_displayed_count = 0
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))
    display_window_size, half_display_window = TENSOR_WINDOW_SIZE, TENSOR_WINDOW_SIZE // 2

    for variant_key in sorted_variant_keys:
        if variants_displayed_count >= max_variants_to_display:
            print(
                f"\n  ... (and {len(node_data_for_display_view) - variants_displayed_count} more variants for node {node_id_str_for_display} not shown)")
            break
        variant_data = node_data_for_display_view[variant_key]
        v_pos, v_type = int(variant_key.split('_')[0]), variant_key.split('_')[1]
        window_center = v_pos + 1 if v_type == 'I' else v_pos
        window_start = max(0, window_center - half_display_window)

        print(f"\n--- Variant: {variant_key} (Node Pos: {v_pos}, Type: {v_type}) ---")
        print(f"  Display Window: {window_start}-{window_start + display_window_size - 1}")
        for k, v_name in [("alt_allele_count", "Alt"), ("ref_allele_count_at_locus", "Ref"),
                          ("other_allele_count_at_locus", "Other"), ("coverage_at_locus", "Cov")]:
            print(f"  {v_name} Count: {variant_data.get(k, 'N/A')}", end=" | ")
        alt_freq = variant_data.get('alt_allele_frequency', 'N/A')
        print(f"Alt Freq: {alt_freq:.4f}" if isinstance(alt_freq, float) else f"Alt Freq: {alt_freq}")

        ref_disp, marker_disp = [' '] * display_window_size, [' '] * display_window_size
        var_idx_win = (v_pos - window_start) if window_start <= v_pos < window_start + display_window_size else -1
        for i in range(display_window_size):
            abs_p = window_start + i
            if 0 <= abs_p < len(full_node_sequence): ref_disp[i] = full_node_sequence[abs_p]
            if i == var_idx_win:
                marker_disp[i] = "I" if v_type == 'I' else "^"
                if v_type == 'I' and i + 1 < display_window_size:
                    marker_disp[i + 1] = "^"
                elif v_type == 'I' and i == display_window_size - 1:
                    marker_disp[i] = ">"
        print(f"  Node Ref: {''.join(ref_disp)}")
        print(f"  Marker  : {''.join(marker_disp)}")

        pileup_reads = variant_data.get("pileup_reads_data", [])
        if not pileup_reads:
            print("  (No reads in window for display)")
        else:
            for i, read_info in enumerate(pileup_reads):
                if i >= max_reads_to_display_per_variant:
                    print(f"  ... ({len(pileup_reads) - i} more reads not shown)")
                    break
                bases_str = "".join([INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in read_info["bases"]])
                print(
                    f"  Read {i + 1:3d}: {bases_str}  (Off:{read_info['offset']},Str:{read_info['strand']},CIG:{read_info.get('cigar', 'N/A')})")
        variants_displayed_count += 1
    print()


# ─────────────────────────────────────────────────────────────────────────────
# Main function
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered .npy tensors and JSON summaries for specified node(s) using parallel processing.",
        # Modified description
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="Base output directory")
    input_node_group = parser.add_mutually_exclusive_group(required=True)
    input_node_group.add_argument("--node_id", type=int, help="Specific node ID to process.")
    input_node_group.add_argument("--node_id_file", help="File with node IDs (one per line).")
    parser.add_argument("--gfa", help="GFA graph file path.")
    parser.add_argument("--load-cache", help="Load node sequences from JSON cache.")
    parser.add_argument("--save-cache", help="Save/update node sequences to JSON cache.")
    parser.add_argument("--num_workers", type=int, default=None,
                        help="Number of worker processes. Defaults to CPU cores.")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N_VARIANTS',
                        help="Print pileups. 0 to disable all per-node stdout and pileup views. No value or -1 for all variants. N (>0) for first N variants/node.")
    parser.add_argument("--max_view_reads", type=int, default=20, help="Max reads per pileup in view.")
    parser.add_argument("--min_af", type=float, default=0.1, help="Min allele frequency for variant processing.")
    args = parser.parse_args()

    for f_path in [args.dat, args.idx]:
        if not os.path.isfile(f_path): sys.exit(f"❌ Error: File not found: {f_path}")
    if not args.load_cache and not args.gfa: sys.exit("❌ Must provide --gfa or --load-cache.")
    if args.load_cache and not os.path.isfile(args.load_cache) and os.path.exists(args.load_cache):
        sys.exit(f"❌ Cache path '{args.load_cache}' is not a file.")
    if args.gfa and not os.path.isfile(args.gfa): sys.exit(f"❌ GFA file not found: {args.gfa}")
    if not (0.0 <= args.min_af <= 1.0): sys.exit("❌ --min_af must be between 0.0 and 1.0.")

    num_workers = args.num_workers if args.num_workers and args.num_workers > 0 else (os.cpu_count() or 1)
    # This print remains as it describes the overall setup
    print(f"🔹 Using {num_workers} worker process(es) for parallel processing.")
    os.makedirs(args.output, exist_ok=True)
    print(f"🔹 Base output directory: {args.output}")

    target_node_ids_set = set()
    if args.node_id_file:
        try:
            with open(args.node_id_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        try:
                            target_node_ids_set.add(int(line))
                        except ValueError:
                            sys.stderr.write(f"⚠️ Invalid node ID '{line}'. Skipping.\n")
            if not target_node_ids_set: sys.exit(f"❌ No valid IDs in {args.node_id_file}.")
            print(f"🔹 Will process {len(target_node_ids_set)} unique ID(s) from {args.node_id_file}")
        except FileNotFoundError:
            sys.exit(f"❌ Node ID file not found: {args.node_id_file}")
    elif args.node_id is not None:
        target_node_ids_set.add(args.node_id)
        print(f"🔹 Will process single ID: {args.node_id}")
    if not target_node_ids_set: sys.exit("ℹ️ No target IDs. Exiting.")

    node_sequences_map = {}
    if args.load_cache and os.path.isfile(args.load_cache):
        s_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                node_sequences_map = json.load(cf)
            print(f"✔ Loaded {len(node_sequences_map)} sequences from cache in {time.time() - s_time:.2f}s.")
        except Exception as e:
            sys.stderr.write(f"⚠️ Error loading cache: {e}.\n")

    nodes_needing_gfa_sequence = {nid for nid in target_node_ids_set if str(nid) not in node_sequences_map}
    if nodes_needing_gfa_sequence and args.gfa:
        print(f"🔹 {len(nodes_needing_gfa_sequence)} node(s) require GFA sequence fetching.")
        s_time = time.time()
        fetched = load_multiple_node_sequences_from_gfa(args.gfa, nodes_needing_gfa_sequence)
        node_sequences_map.update(fetched)
        print(
            f"✔ Fetched {len(fetched)} new sequences from GFA in {time.time() - s_time:.2f}s. Total map size: {len(node_sequences_map)}.")
    elif nodes_needing_gfa_sequence:
        sys.stderr.write(f"⚠️ {len(nodes_needing_gfa_sequence)} nodes need GFA sequences, but --gfa not provided.\n")

    idx_load_start_time = time.time()
    full_idx_map = load_full_idx_data(args.idx)
    if full_idx_map is None: sys.exit(f"❌ Failed to load index data from {args.idx}.")
    print(f"✔ Index data ({len(full_idx_map)} entries) loaded in {time.time() - idx_load_start_time:.2f}s.")

    tasks_to_submit, skipped_nodes_pre_submit = [], set()
    print(f"🔹 Preparing tasks for {len(target_node_ids_set)} target nodes...")
    task_prep_start_time = time.time()
    for i, node_id_int in enumerate(target_node_ids_set):
        if (i + 1) % 50000 == 0: print(f"  Prepared tasks for {i + 1}/{len(target_node_ids_set)} nodes...")
        node_seq = node_sequences_map.get(str(node_id_int))
        node_dat_info = full_idx_map.get(node_id_int)
        if not node_seq or not node_dat_info:
            skipped_nodes_pre_submit.add(node_id_int)
            continue
        tasks_to_submit.append((node_id_int, node_dat_info[0], node_dat_info[1], node_seq, args.min_af))
    print(f"✔ Task preparation completed in {time.time() - task_prep_start_time:.2f}s.")
    if skipped_nodes_pre_submit:
        print(
            f"⚠️ Skipped {len(skipped_nodes_pre_submit)} nodes (missing seq/index). Examples: {list(skipped_nodes_pre_submit)[:3]}")
    if not tasks_to_submit: sys.exit("ℹ️ No valid tasks to process.")

    overall_parallel_start_time = time.time()  # Renamed for clarity
    total_npy_files_generated = 0  # Renamed from total_pth
    processed_nodes_count = 0
    successful_nodes_with_output = 0
    results_for_viewing = {}

    # For batch reporting in parallel phase
    batch_start_time_parallel = time.time()
    nodes_in_current_batch_parallel = 0
    npy_in_current_batch_parallel = 0

    verbose_node_output_parallel = not (args.view == 0)  # For per-node messages from main thread if any were needed

    print(f"\n🔹 Submitting {len(tasks_to_submit)} node tasks to {num_workers} worker(s)...")
    with ProcessPoolExecutor(max_workers=num_workers, initializer=init_worker,
                             initargs=(args.dat, args.output)) as executor:
        future_to_node_id = {executor.submit(process_single_node_for_pileup, task): task[0] for task in tasks_to_submit}

        for future_idx, future in enumerate(as_completed(future_to_node_id)):
            current_completed_total = future_idx + 1
            orig_node_id = future_to_node_id[future]

            try:
                res_node_id, view_data, npy_count_for_node = future.result()  # Expecting npy_count now
                processed_nodes_count += 1
                if res_node_id is None:
                    sys.stderr.write(f"❌ Worker failed for node {orig_node_id} (returned None ID).\n")
                else:
                    total_npy_files_generated += npy_count_for_node
                    npy_in_current_batch_parallel += npy_count_for_node
                    summary_file = os.path.join(args.output, str(res_node_id), "variant_summary.json")
                    if npy_count_for_node > 0 or os.path.exists(summary_file):
                        successful_nodes_with_output += 1
                    if args.view is not None and args.view != 0 and view_data:  # Only collect if view is active and not 0
                        seq_for_view = node_sequences_map.get(str(res_node_id))
                        if seq_for_view: results_for_viewing[res_node_id] = (view_data, seq_for_view)
            except Exception as exc:
                processed_nodes_count += 1  # Still counts as a task that completed (with error)
                sys.stderr.write(f"❌ Error processing node {orig_node_id} (future exception): {exc}\n")

            nodes_in_current_batch_parallel += 1
            if nodes_in_current_batch_parallel == 1000 or current_completed_total == len(tasks_to_submit):
                if nodes_in_current_batch_parallel > 0:
                    batch_duration = time.time() - batch_start_time_parallel
                    rate = nodes_in_current_batch_parallel / batch_duration if batch_duration > 0 else 0
                    print(
                        f"  Processed batch of {nodes_in_current_batch_parallel} nodes ({current_completed_total}/{len(tasks_to_submit)} total) "
                        f"in {batch_duration:.2f}s ({rate:.2f} nodes/sec). "
                        f"Generated {npy_in_current_batch_parallel} .npy files in this batch.")
                    batch_start_time_parallel = time.time()
                    nodes_in_current_batch_parallel = 0
                    npy_in_current_batch_parallel = 0

    if args.view is not None and args.view != 0 and results_for_viewing:  # Check args.view != 0 here
        print("\n══════════ VIEWING PILEUPS ══════════")
        for node_id_view in sorted(results_for_viewing.keys()):
            view_data, node_seq = results_for_viewing[node_id_view]
            max_v = float('inf') if args.view == -1 else args.view
            if args.view < -1: max_v = float('inf')  # ensure correct default for other negatives
            if max_v > 0 or max_v == float('inf'):  # Only call display if we intend to show something
                display_pileup_data(view_data, str(node_id_view), node_seq, args.max_view_reads, max_v)
    elif args.view is not None and args.view != 0:  # View was requested (not 0) but no results
        print(f"ℹ️ --view specified, but no pileup data gathered for display (or all nodes failed).")
    elif args.view == 0:
        print(f"ℹ️ --view 0 specified: Pileup display disabled for all nodes.")

    print("\n══════════ PROCESSING COMPLETE ══════════")
    if args.save_cache and node_sequences_map:
        print(f"\n🔹 Saving {len(node_sequences_map)} sequences to cache: {args.save_cache}...")
        try:
            with open(args.save_cache, 'w') as wcf:
                json.dump(node_sequences_map, wcf, indent=2)
            print(f"✔ Sequences saved to cache.")
        except Exception as e:
            sys.stderr.write(f"❌ Error saving to cache {args.save_cache}: {e}\n")
    elif args.save_cache:
        print(f"ℹ️ --save-cache: No sequences in map to save.")

    print(f"\nSummary:")
    print(f"  Targeted: {len(target_node_ids_set)} unique node IDs.")
    if skipped_nodes_pre_submit: print(f"  Skipped pre-submission: {len(skipped_nodes_pre_submit)} node(s).")
    print(f"  Submitted to workers: {len(tasks_to_submit)} node(s).")
    print(f"  Worker tasks completed (incl. errors): {processed_nodes_count}/{len(tasks_to_submit)}.")
    print(f"  Output files generated for: {successful_nodes_with_output} node(s).")
    print(f"  Total .npy tensor files: {total_npy_files_generated}.")  # Updated to .npy
    print(f"🏁 Parallel processing phase finished in {time.time() - overall_parallel_start_time:.2f} seconds.")


if __name__ == '__main__':
    main()
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
# ProcessPoolExecutor removed as this version is serial
import torch  # Still used for tensor creation before converting to NumPy
import os

# ─────────────────────────────────────────────────────────────────────────────
# Dynamic per-node record struct (support NEW + OLD .dat formats)

def make_record_struct_old(node_length: int) -> struct.Struct:
    """
    OLD format per-read record (variable by node_length):
      <h {L}s {L}s {L}s h c>
        - i16 offset
        - seq[L]
        - bq[L]
        - cigar[L]
        - i16 mapq
        - char strand
    """
    return struct.Struct(f"<h{node_length}s{node_length}s{node_length}shc")

def make_record_struct_latest(max_read_len: int, max_cigar_len: int) -> struct.Struct:
    """
    LATEST format per-read record (variable by per-block maxima):
      <h {R}s {R}s {C}s h c>
        - i16 offset
        - seq[R]
        - bq[R]
        - cigar[C]
        - i16 mapq
        - char strand
    """
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")

# .dat block headers
BLOCK_HDR_PACK_LATEST = struct.Struct("<I I H I I")   # 18 bytes: node_id, n_records, flags, max_read_len, max_cigar_len
BLOCK_HDR_PACK_PADDED = struct.Struct("<I I H 2x I")  # 16 bytes: node_id, n_records, flags, pad, node_length
BLOCK_HDR_PACK_14B    = struct.Struct("<I I H I")     # 14 bytes: node_id, n_records, flags, node_length

# ─────────────────────────────────────────────────────────────────────────────
# Constants (unchanged semantics for downstream code)
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '*': 5, ' ': 6, '-': 6}
INDEX_TO_BASE_FOR_VIEW = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N', 5: '*', 6: ' '}  # As provided by user

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
PADDING_BASE_INDEX = BASE_TO_INDEX[' ']
DEFAULT_QUALITY_PADDING = 0
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def _read_block_header(dat_file, block_start):
    """
    Read .dat block header at block_start; returns:
      (node_id, n_records, flags, lengths, header_size_bytes, fmt)
    where:
      - fmt == 'latest' and lengths == (R, C)  for 18B header
      - fmt == 'old'    and lengths == L       for 16B/14B headers
    """
    # Try latest 18B (R,C)
    dat_file.seek(block_start, os.SEEK_SET)
    hdr = dat_file.read(BLOCK_HDR_PACK_LATEST.size)
    if len(hdr) == BLOCK_HDR_PACK_LATEST.size:
        nid, nrec, flags, R, C = BLOCK_HDR_PACK_LATEST.unpack(hdr)
        if (1 <= R <= 1_000_000) and (1 <= C <= 1_000_000) and (0 <= nrec < 10_000_000):
            return nid, nrec, flags, (R, C), BLOCK_HDR_PACK_LATEST.size, 'latest'

    # Try old 16B padded (L)
    dat_file.seek(block_start, os.SEEK_SET)
    hdr = dat_file.read(BLOCK_HDR_PACK_PADDED.size)
    if len(hdr) == BLOCK_HDR_PACK_PADDED.size:
        nid, nrec, flags, L = BLOCK_HDR_PACK_PADDED.unpack(hdr)
        if (1 <= L <= 1_000_000) and (0 <= nrec < 10_000_000):
            return nid, nrec, flags, L, BLOCK_HDR_PACK_PADDED.size, 'old'

    # Fallback: old 14B (L)
    dat_file.seek(block_start, os.SEEK_SET)
    hdr = dat_file.read(BLOCK_HDR_PACK_14B.size)
    if len(hdr) != BLOCK_HDR_PACK_14B.size:
        raise RuntimeError(f"Cannot read block header at offset {block_start}")
    nid, nrec, flags, L = BLOCK_HDR_PACK_14B.unpack(hdr)
    if not (1 <= L <= 1_000_000 and 0 <= nrec < 10_000_000):
        raise RuntimeError(f"Suspicious block header at {block_start}: node_len={L}, nrec={nrec}")
    return nid, nrec, flags, L, BLOCK_HDR_PACK_14B.size, 'old'

def load_full_idx_data(idx_path):
    """
    Supports:
      • Latest fixed-size entry: 30B  <I Q I I H I I>
      • New fixed-size entry:    26B  <I Q I I H I>
      • Older fixed-size entry:  22B  <I Q I I H>
      • Legacy variable-size entries with metadata_len (fallback)
    Returns: dict { node_id: (offset, n_records) }
    """
    idx_data_map = {}
    print(f"🔹 Loading full index data from {idx_path}...")
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                sys.stderr.write(f"❌ Error: Index file {idx_path} is too small.\n")
                return None
            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                sys.stderr.write(f"❌ Error: Could not read number of nodes from {idx_path}.\n")
                return None
            num_nodes_in_idx = struct.unpack('<I', num_nodes_bytes)[0]
            if num_nodes_in_idx == 0:
                print("✔ Index has 0 entries.")
                return idx_data_map

            remaining = file_size - 4
            # Choose strategy by per-entry size if possible
            if remaining % num_nodes_in_idx == 0:
                per = remaining // num_nodes_in_idx
                if per == 30:
                    strategy = "fixed30"
                elif per == 26:
                    strategy = "fixed26"
                elif per == 22:
                    strategy = "fixed22"
                else:
                    strategy = "legacy"
            else:
                strategy = "legacy"

            print(f"  Detected index entry layout: {strategy}")

            for i in range(num_nodes_in_idx):
                if strategy == "fixed30":
                    rec = f.read(30)
                    if len(rec) != 30:
                        sys.stderr.write(f"❌ Error: Truncated 30B index entry at {i}.\n"); break
                    node_id, offset, _bsize, nrec, _flags, _R, _C = struct.unpack("<I Q I I H I I", rec)
                elif strategy == "fixed26":
                    rec = f.read(26)
                    if len(rec) != 26:
                        sys.stderr.write(f"❌ Error: Truncated 26B index entry at {i}.\n"); break
                    node_id, offset, _bsize, nrec, _flags, _nlen = struct.unpack("<I Q I I H I", rec)
                elif strategy == "fixed22":
                    rec = f.read(22)
                    if len(rec) != 22:
                        sys.stderr.write(f"❌ Error: Truncated 22B index entry at {i}.\n"); break
                    node_id, offset, _bsize, nrec, _flags = struct.unpack("<I Q I I H", rec)
                else:
                    # Legacy with metadata_len
                    rec = f.read(22)
                    if len(rec) != 22:
                        sys.stderr.write(f"❌ Error: Truncated legacy index header at {i}.\n"); break
                    node_id, offset, _bsize, nrec, meta_len = struct.unpack("<I Q I I H", rec)
                    if meta_len:
                        skipped = f.read(meta_len)
                        if len(skipped) != meta_len:
                            sys.stderr.write(f"❌ Error: Truncated legacy metadata at {i}.\n"); break
                idx_data_map[node_id] = (offset, nrec)
                if (i + 1) % 2_000_000 == 0:
                    print(f"    Loaded {i + 1}/{num_nodes_in_idx} index entries...")

            print(f"✔ Successfully loaded {len(idx_data_map)} entries.")
        return idx_data_map
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: Index file not found: {idx_path}\n")
        return None
    except Exception as e:
        sys.stderr.write(f"❌ Error parsing index {idx_path}: {e}\n")
        return None

def load_multiple_node_sequences_from_gfa(gfa_path, target_node_ids_set):
    node_sequences = {}
    if not target_node_ids_set: return node_sequences
    nodes_to_find = target_node_ids_set.copy()
    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file for {len(nodes_to_find)} target node(s): {gfa_path}")
            line_counter = 0
            for line in f:
                line_counter += 1
                if line_counter % 10_000_000 == 0:
                    print(f"  Checked {line_counter:,} GFA lines... {len(nodes_to_find)} nodes remaining.")
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
                            f"✔ Found all {len(target_node_ids_set)} requested node sequences in GFA ({line_counter:,} lines checked).")
                        break
            if nodes_to_find:
                print(
                    f"✔ GFA scan complete ({line_counter:,} lines). Found {len(node_sequences)}/{len(target_node_ids_set)} sequences.")
                if nodes_to_find:
                    print(f"⚠️ Missing example IDs: {list(nodes_to_find)[:5]}")
    except FileNotFoundError:
        sys.stderr.write(f"❌ Error: GFA file not found: {gfa_path}\n")
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
    except Exception:
        return []
    return ops

def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None, expected_ref_allele_for_indel=None):
    current_node_pos = read_offset_on_node
    current_read_pos = 0

    # Track whether we've covered the insertion anchor in a match;
    # if no 'I' follows at this anchor, we’ll return REF at the end.
    saw_insertion_anchor = False

    for length, op in read_cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            if current_node_pos <= target_node_pos < current_node_pos + length:
                offset_in_block = target_node_pos - current_node_pos
                read_idx = current_read_pos + offset_in_block

                if read_idx < len(read_sequence):
                    allele = read_sequence[read_idx].upper()

                    if expected_var_type == 'D':
                        # Deletion: covering base ⇒ REF (no deletion spanning it)
                        return "REF_STATE_FOR_INDEL"

                    if expected_var_type == 'I':
                        # Insertion: mark anchor seen and KEEP SCANNING for a following 'I'
                        saw_insertion_anchor = True
                        # do not return yet
                    else:
                        # SNP
                        return allele
                else:
                    # Out of read range; still treat as anchor seen for insertions
                    if expected_var_type == 'I':
                        saw_insertion_anchor = True
                    else:
                        return None
            current_node_pos += length
            current_read_pos += length

        elif op == 'I':
            # Insertion occurs BETWEEN bases; its anchor is (current_node_pos - 1)
            if expected_var_type == 'I' and (current_node_pos - 1) == target_node_pos:
                return read_sequence[current_read_pos: current_read_pos + length].upper()
            current_read_pos += length  # ref doesn't advance

        elif op == 'D':
            if current_node_pos <= target_node_pos < current_node_pos + length:
                if expected_var_type == 'I':
                    return "OTHER_FOR_INDEL"
                if expected_var_type == 'D':
                    ref_deleted_segment = node_sequence[current_node_pos: current_node_pos + length]
                    if ref_deleted_segment == expected_ref_allele_for_indel:
                        return "*"
                    else:
                        return "OTHER_FOR_INDEL"
                return "*"
            current_node_pos += length

        elif op == 'S':
            current_read_pos += length
        elif op == 'N':
            current_node_pos += length

        # Preserve original early-stop condition
        if current_node_pos > target_node_pos + 1 and op in ('M', '=', 'X', 'D', 'N'):
            if not (expected_var_type == 'I' and (current_node_pos - 1) <= target_node_pos):
                break

    # If we saw the insertion anchor but never encountered an 'I' at that anchor → REF
    if expected_var_type == 'I' and saw_insertion_anchor:
        return "REF_STATE_FOR_INDEL"

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
                    if node_base != read_base and read_base != 'N' and node_base != 'N':
                        variants.append((cur_node_p, 'X', read_base, node_base))
                else:
                    break
            node_pos += length; read_pos += length
        elif op == 'I':
            ins_seq = read_sequence[read_pos: read_pos + length].upper()
            anchor_pos = node_pos - 1 if node_pos > 0 else 0
            anchor_base = node_sequence[anchor_pos].upper() if 0 <= anchor_pos < node_seq_len else "*"
            if ins_seq and 'N' not in ins_seq:
                variants.append((anchor_pos, 'I', ins_seq, anchor_base))
            read_pos += length
        elif op == 'D':
            del_seq = node_sequence[node_pos: node_pos + length].upper() if node_pos + length <= node_seq_len else ""
            if del_seq and 'N' not in del_seq:
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
            node_pos += L; read_pos += L
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
                                   segment_read_sequence, segment_quality_values,
                                   window_start_node, tensor_win_size, node_len):
    """Accepts quality values as a list/bytes of ints (raw), or str (ASCII PHRED)."""
    bases = [PADDING_BASE_INDEX] * tensor_win_size
    quals = [DEFAULT_QUALITY_PADDING] * tensor_win_size
    node_pos, read_pos = segment_offset_on_node, 0

    def _qual_at(idx):
        if isinstance(segment_quality_values, (bytes, bytearray, list, tuple)):
            if 0 <= idx < len(segment_quality_values):
                try:
                    return int(segment_quality_values[idx])
                except Exception:
                    return DEFAULT_QUALITY_PADDING
            return DEFAULT_QUALITY_PADDING
        # fallback: string of ASCII +33
        if 0 <= idx < len(segment_quality_values):
            try:
                return max(0, ord(segment_quality_values[idx]) - 33)
            except Exception:
                return DEFAULT_QUALITY_PADDING
        return DEFAULT_QUALITY_PADDING

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
                    quals[win_idx] = _qual_at(r_aln)
            node_pos += L; read_pos += L
        elif op in ('D', 'N'):
            for i in range(L):
                n_aln = node_pos + i
                if window_start_node <= n_aln < window_start_node + tensor_win_size:
                    win_idx = n_aln - window_start_node
                    bases[win_idx] = BASE_TO_INDEX['*']  # quals remain default (padding)
            node_pos += L
        elif op in ('I', 'S'):
            read_pos += L
        if read_pos >= len(segment_read_sequence) and op in ('M', '=', 'X', 'I', 'S'): break
    return bases, quals

# ─────────────────────────────────────────────────────────────────────────────
# Core Processing Logic (Serial)
def process_node_serially(dat_file_path, base_output_dir,
                          node_id, dat_file_offset, n_records_from_idx, node_sequence,
                          min_af_threshold, variant_type):
    npy_files_generated = 0
    view_oriented_variant_data = {}
    node_specific_output_dir = os.path.join(base_output_dir, str(node_id))
    try:
        os.makedirs(node_specific_output_dir, exist_ok=True)
    except OSError as e:
        sys.stderr.write(f"❌ Error creating dir {node_specific_output_dir} for Node {node_id}: {e}\n")
        return node_id, None, npy_files_generated

    if not node_sequence:
        sys.stderr.write(f"ℹ️ Node {node_id}: No sequence. Skipping.\n")
        return node_id, {}, npy_files_generated

    node_len = len(node_sequence)
    aligned_read_segments = []
    try:
        with open(dat_file_path, 'rb') as dat_f:
            # Read block header at offset to recover lengths (and authoritative n_records)
            nid2, n_records, _flags, lengths, hdr_size, fmt = _read_block_header(dat_f, dat_file_offset)

            # Build dynamic record struct
            if fmt == 'latest':
                R, C = lengths
                rec_struct = make_record_struct_latest(int(R), int(C))
            else:
                L = int(lengths)
                rec_struct = make_record_struct_old(L)

            rec_size = rec_struct.size

            # Seek to records start
            records_start = dat_file_offset + hdr_size
            dat_f.seek(records_start, os.SEEK_SET)

            for _ in range(n_records):
                data = dat_f.read(rec_size)
                if len(data) < rec_size: break
                off, raw_seq, raw_qual, raw_cigar, mapq, strand_b = rec_struct.unpack(data)
                if mapq < 10:  # same filter
                    continue
                try:
                    seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                    qual_vals = list(raw_qual.rstrip(b'\0'))  # raw bytes → list[int]
                    cigar_orig = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                    strand = strand_b.decode('ascii') if isinstance(strand_b, (bytes, bytearray)) else chr(strand_b)
                except UnicodeDecodeError:
                    continue
                if not cigar_orig or len(seq) != len(qual_vals):
                    continue

                cigar_ops_orig = decode_cigar_to_int_ops(cigar_orig)
                if not cigar_ops_orig and cigar_orig != '*': continue

                cur_seq, cur_quals, cur_cigar_ops, cur_offset = seq, qual_vals, list(cigar_ops_orig), off

                if strand == '-':
                    cur_seq = reverse_complement(seq)
                    cur_quals = qual_vals[::-1]
                    cur_cigar_ops = [op for op in reversed(
                        cigar_ops_orig)] if cigar_ops_orig else []
                    # alignment_span_on_node = len(cur_seq)
                    alignment_span_on_node = len(cur_seq)
                    # Adjust span: +I, -D (reference-consuming vs non-consuming edits)
                    for Lx, opx in cigar_ops_orig:
                        if opx == 'I':
                            alignment_span_on_node -= Lx
                        elif opx == 'D':
                            alignment_span_on_node += Lx
                    cur_offset = node_len - alignment_span_on_node - off
                    if cur_offset < 0: continue

                aligned_read_segments.append({
                    "offset_on_node": cur_offset,
                    "read_sequence": cur_seq,
                    "processed_quality_values": cur_quals,
                    "cigar_ops": cur_cigar_ops,
                    "original_cigar_str": cigar_orig,
                    "strand": strand
                })
    except FileNotFoundError:
        sys.stderr.write(f"❌ DAT file {dat_file_path} not found for Node {node_id}.\n")
        return node_id, None, npy_files_generated
    except Exception as e:
        sys.stderr.write(f"❌ Error reading DAT for Node {node_id}: {e}\n")
        return node_id, None, npy_files_generated

    # Build candidate variants
    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            # Apply variant-type selection early
            if variant_type == 'snp' and v_type != 'X':
                continue
            if variant_type == 'indel' and v_type not in ('I', 'D'):
                continue
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    variant_headers = []
    half_win = TENSOR_WINDOW_SIZE // 2

    for (v_pos, v_type, v_ref, v_alt), _ in candidate_variants.items():
        alt_c, ref_c, other_c, locus_cov = 0, 0, 0, 0
        indel_ref_check = None
        if v_type == 'D':
            indel_ref_check = v_ref
        elif v_type == 'I' and 0 <= v_pos < node_len:
            indel_ref_check = node_sequence[v_pos]

        for seg in aligned_read_segments:
            allele = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["cigar_ops"],
                v_pos, node_sequence, v_type, indel_ref_check if v_type in 'ID' else v_ref)
            if allele is not None:
                locus_cov += 1
                if v_type == 'X':
                    if allele == v_alt:
                        alt_c += 1
                    elif allele == v_ref:
                        ref_c += 1
                    else:
                        other_c += 1
                elif v_type == 'I':
                    if allele == v_alt:
                        alt_c += 1
                    elif allele == "REF_STATE_FOR_INDEL":
                        ref_c += 1
                    else:
                        other_c += 1
                elif v_type == 'D':
                    if allele == "*":
                        alt_c += 1
                    elif allele == "REF_STATE_FOR_INDEL":
                        ref_c += 1
                    else:
                        other_c += 1

        alt_freq = alt_c / locus_cov if locus_cov > 0 else 0.0
        if alt_freq < min_af_threshold:
            continue

        # Filename key
        if v_type == 'D':
            deleted_seq = v_alt if (v_alt and v_alt != '*') else v_ref  # pick the non-* token
            anchor_base = node_sequence[v_pos] if 0 <= v_pos < node_len else 'N'
            key_str = f"{v_pos}_{v_type}_{anchor_base}{deleted_seq}_{anchor_base}"
        else:
            key_str = f"{v_pos}_{v_type}_{v_ref}_{v_alt}"

        win_center = v_pos + 1 if v_type == 'I' else v_pos
        win_start = max(0, win_center - half_win)

        # View data (pileup) rows
        view_reads_data = []
        for seg in aligned_read_segments[:TENSOR_MAX_READ_ROWS]:
            row_chars = get_read_representation_in_window_for_view(
                seg["cigar_ops"], seg["offset_on_node"], seg["read_sequence"],
                win_start, TENSOR_WINDOW_SIZE, node_len)
            if any(c != ' ' for c in row_chars):
                view_reads_data.append({
                    "bases": [BASE_TO_INDEX.get(c.upper(), BASE_TO_INDEX['N']) for c in row_chars],
                    "offset": seg["offset_on_node"], "strand": seg["strand"], "cigar": seg["original_cigar_str"]})

        view_oriented_variant_data[key_str] = {
            "pileup_reads_data": view_reads_data, "alt_allele_count": alt_c,
            "ref_allele_count_at_locus": ref_c, "other_allele_count_at_locus": other_c,
            "coverage_at_locus": locus_cov, "alt_allele_frequency": round(alt_freq, 4)}

        # Build tensors
        ch1_bases, ch2_quals, ch3_mismatches = [], [], []
        ref_bases_tensor = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i in range(TENSOR_WINDOW_SIZE):
            abs_pos = win_start + i
            if 0 <= abs_pos < node_len:
                ref_bases_tensor[i] = BASE_TO_INDEX.get(node_sequence[abs_pos].upper(), BASE_TO_INDEX['N'])
        ch1_bases.append(ref_bases_tensor)
        ch2_quals.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
        ch3_mismatches.append([MISMATCH_CHANNEL_REF_ROW_VALUE] * TENSOR_WINDOW_SIZE)

        tensor_reads_added = 0
        for seg in aligned_read_segments:
            if tensor_reads_added >= TENSOR_MAX_READ_ROWS: break
            base_row, qual_row = get_read_tensor_rows_in_window(
                seg["cigar_ops"], seg["offset_on_node"], seg["read_sequence"],
                seg["processed_quality_values"], win_start, TENSOR_WINDOW_SIZE, node_len)
            if any(b != PADDING_BASE_INDEX for b in base_row):
                ch1_bases.append(base_row)
                ch2_quals.append(qual_row)
                mismatch_row = [MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE
                for i in range(TENSOR_WINDOW_SIZE):
                    r_idx, ref_idx = base_row[i], ref_bases_tensor[i]
                    if r_idx == PADDING_BASE_INDEX or ref_idx == PADDING_BASE_INDEX: continue
                    mismatch_row[i] = 0 if r_idx == ref_idx else 1
                ch3_mismatches.append(mismatch_row)
                tensor_reads_added += 1

        # pad to fixed height
        for _ in range(TENSOR_MAX_READ_ROWS - tensor_reads_added):
            ch1_bases.append([PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE)
            ch2_quals.append([DEFAULT_QUALITY_PADDING] * TENSOR_WINDOW_SIZE)
            ch3_mismatches.append([MISMATCH_COMPARISON_PADDING_VALUE] * TENSOR_WINDOW_SIZE)

        try:
            tensor = torch.tensor([ch1_bases, ch2_quals, ch3_mismatches], dtype=torch.int8)
            np_array = tensor.numpy()
            fname = f"{key_str}.npy"
            np.save(os.path.join(node_specific_output_dir, fname), np_array)
            variant_headers.append({
                "variant_key": key_str, "tensor_file": fname, "alt_allele_count": alt_c,
                "ref_allele_count_at_locus": ref_c, "other_allele_count_at_locus": other_c,
                "coverage_at_locus": locus_cov, "alt_allele_frequency": round(alt_freq, 4)})
        except Exception as e:
            sys.stderr.write(f"❌ Node {node_id}: Tensor error for {key_str}: {e}\n")

    npy_files_generated = len(variant_headers)
    if variant_headers:
        summary_path = os.path.join(node_specific_output_dir, "variant_summary.json")
        try:
            with open(summary_path, 'w') as f:
                json.dump({"node_id": node_id, "node_length": node_len,
                           "node_sequence_preview": node_sequence[:100] + ("..." if node_len > 100 else ""),
                           "variants": variant_headers}, f, indent=2)
        except Exception as e:
            sys.stderr.write(f"❌ Node {node_id}: Summary JSON error: {e}\n")
    return node_id, view_oriented_variant_data, npy_files_generated

# ─────────────────────────────────────────────────────────────────────────────
# Pileup Viewing Function (unchanged)
def display_pileup_data(node_data_for_display_view, node_id_str_for_display, full_node_sequence,
                        max_reads_to_display_per_variant, max_variants_to_display=float('inf')):
    if not node_data_for_display_view:
        print(f"ℹ️ No valid pileup data to display for node {node_id_str_for_display}.", file=sys.stderr)
        return
    print(f"\n=== Displaying Pileups for Node ID: {node_id_str_for_display} (Length: {len(full_node_sequence)}) ===")

    variants_displayed_count = 0
    sorted_variant_keys = sorted(node_data_for_display_view.keys(),
                                 key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))
    display_window_size, half_display_window = TENSOR_WINDOW_SIZE, TENSOR_WINDOW_SIZE // 2

    for variant_key in sorted_variant_keys:
        if variants_displayed_count >= max_variants_to_display:
            print(
                f"\n  ... ({len(node_data_for_display_view) - variants_displayed_count} more variants for node {node_id_str_for_display} not shown)")
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
# Main
def main():
    parser = argparse.ArgumentParser(
        description="Serial variant tensor generator (supports new and old .dat/.idx formats; select SNP/INDEL).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat file path")
    parser.add_argument("idx", help=".idx file path")
    parser.add_argument("output", help="Base output directory")
    input_node_group = parser.add_mutually_exclusive_group(required=True)
    input_node_group.add_argument("--node_id", type=int, help="Specific node ID to process.")
    input_node_group.add_argument("--node_id_file", help="File with node IDs (one per line).")
    parser.add_argument("--gfa", help="GFA graph file path.")
    parser.add_argument("--load-cache", help="Load sequences from JSON cache.")
    parser.add_argument("--save-cache", help="Save/update sequences to JSON cache.")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N_VARIANTS',
                        help="Print pileups. 0 to disable. No value or -1 for all. N for first N variants/node.")
    parser.add_argument("--max_view_reads", type=int, default=20, help="Max reads per pileup in view.")
    parser.add_argument("--min_af", type=float, default=0.1, help="Min allele frequency for processing.")
    parser.add_argument("--variant_type", type=str, default="all", choices=["snp", "indel", "all"],
                        help="Select 'snp' for mismatches only, 'indel' for I/D only, or 'all'.")
    args = parser.parse_args()

    for f_path in [args.dat, args.idx]:
        if not os.path.isfile(f_path): sys.exit(f"❌ Error: File not found: {f_path}")
    if not args.load_cache and not args.gfa: sys.exit("❌ Must provide --gfa or --load-cache.")
    if args.load_cache and not os.path.isfile(args.load_cache) and os.path.exists(args.load_cache):
        sys.exit(f"❌ Cache path '{args.load_cache}' is not a file.")
    if args.gfa and not os.path.isfile(args.gfa): sys.exit(f"❌ GFA file not found: {args.gfa}")
    if not (0.0 <= args.min_af <= 1.0): sys.exit("❌ --min_af must be between 0.0 and 1.0.")
    os.makedirs(args.output, exist_ok=True)
    print(f"🔹 Base output directory: {args.output}")

    # Variant type mapping
    vtype = args.variant_type
    # Normalize to internal tags
    if vtype == 'snp':
        variant_type = 'snp'
    elif vtype == 'indel':
        variant_type = 'indel'
    else:
        variant_type = 'all'

    target_node_ids = set()
    if args.node_id_file:
        try:
            with open(args.node_id_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        try:
                            target_node_ids.add(int(line))
                        except ValueError:
                            sys.stderr.write(f"⚠️ Invalid node ID '{line}'. Skipping.\n")
            if not target_node_ids: sys.exit(f"❌ No valid IDs in {args.node_id_file}.")
            print(f"🔹 Will process {len(target_node_ids)} unique ID(s) from {args.node_id_file}")
        except FileNotFoundError:
            sys.exit(f"❌ Node ID file not found: {args.node_id_file}")
    elif args.node_id is not None:
        target_node_ids.add(args.node_id)
        print(f"🔹 Will process single ID: {args.node_id}")
    if not target_node_ids: sys.exit("ℹ️ No target IDs. Exiting.")

    overall_start_time = time.time()
    idx_data = load_full_idx_data(args.idx)
    if idx_data is None: sys.exit("❌ Critical error loading IDX data.")

    seq_map = {}
    if args.load_cache and os.path.isfile(args.load_cache):
        s_tm = time.time()
        try:
            with open(args.load_cache, 'r') as f:
                seq_map = json.load(f)
            print(f"✔ Loaded {len(seq_map)} sequences from cache in {time.time() - s_tm:.2f}s.")
        except Exception as e:
            sys.stderr.write(f"⚠️ Error loading cache: {e}. Proceeding without.\n")

    needed_gfa_seqs = {nid for nid in target_node_ids if str(nid) not in seq_map}
    if needed_gfa_seqs and args.gfa:
        s_tm = time.time()
        fetched = load_multiple_node_sequences_from_gfa(args.gfa, needed_gfa_seqs)
        if fetched: seq_map.update(fetched)
        print(f"✔ GFA: Fetched {len(fetched)} new sequences in {time.time() - s_tm:.2f}s. Total in map: {len(seq_map)}.")
    elif needed_gfa_seqs:
        sys.stderr.write(f"⚠️ {len(needed_gfa_seqs)} nodes need GFA sequences, but --gfa not given.\n")

    processed_count, successful_output_count, total_npy_count = 0, 0, 0
    sorted_target_ids = sorted(list(target_node_ids))

    # Serial loop
    batch_start_time_serial = time.time()
    nodes_in_current_batch_serial = 0
    npy_in_current_batch_serial = 0

    for i, node_id_int in enumerate(sorted_target_ids):
        node_id_str = str(node_id_int)
        print(f"\n══════ PROCESSING NODE: {node_id_int} ({i + 1}/{len(sorted_target_ids)}) ══════")
        node_loop_start_time = time.time()

        dat_info = idx_data.get(node_id_int)
        seq = seq_map.get(node_id_str)
        if not dat_info:
            sys.stderr.write(f"❌ Node {node_id_int}: Not in IDX data. Skipping.\n"); continue
        if not seq:
            sys.stderr.write(f"❌ Node {node_id_int}: Sequence unavailable. Skipping.\n"); continue

        offset, n_recs_from_idx = dat_info
        print(f"  Node {node_id_int}: Len={len(seq)}bp, Records(idx)={n_recs_from_idx}, MinAF={args.min_af}, Select={args.variant_type}")

        try:
            _, view_data, npy_this_node = process_node_serially(
                args.dat, args.output, node_id_int, offset, n_recs_from_idx, seq, args.min_af, variant_type)

            processed_count += 1
            total_npy_count += npy_this_node
            npy_in_current_batch_serial += npy_this_node

            summary_path = os.path.join(args.output, node_id_str, "variant_summary.json")
            if npy_this_node > 0 or os.path.exists(summary_path):
                print(f"  ✅ Output for node {node_id_int}. Summary: {summary_path if os.path.exists(summary_path) else 'Not found'}")
                successful_output_count += 1
            elif view_data is not None and not view_data:
                print(f"  ℹ️ No variants met AF for node {node_id_int}.")
            elif view_data is None:
                print(f"  ⚠️ Processing error for node {node_id_int}.")

            if args.view is not None and args.view != 0:
                if view_data:
                    max_v = float('inf') if args.view == -1 else args.view
                    if args.view < -1: max_v = float('inf')
                    if max_v > 0 or max_v == float('inf'):
                        v_msg = f"first {int(max_v)}" if max_v != float('inf') else "all"
                        print(f"  🔹 Displaying {v_msg} pileups for node {node_id_int} (max {args.max_view_reads} reads/var)...")
                        display_pileup_data(view_data, node_id_str, seq, args.max_view_reads, max_v)
                elif view_data is not None:
                    print(f"  ℹ️ --view: No variants met AF for node {node_id_int} to display.")
            elif args.view == 0:
                print(f"  ℹ️ --view 0 specified: Pileup display disabled for node {node_id_int}.")

            print(f"  ✔ Node {node_id_int} processing took {time.time() - node_loop_start_time:.2f}s.")

        except Exception as e:
            sys.stderr.write(f"❌ CRITICAL ERROR for node {node_id_int} in main loop: {e}\n")

        nodes_in_current_batch_serial += 1
        if nodes_in_current_batch_serial == 1000 or (i + 1) == len(sorted_target_ids):
            if nodes_in_current_batch_serial > 0:
                batch_duration_serial = time.time() - batch_start_time_serial
                nodes_per_sec_serial = nodes_in_current_batch_serial / batch_duration_serial if batch_duration_serial > 0 else 0
                print(f"\n--- Batch Report (Serial Processing) ---")
                print(f"  Processed batch of {nodes_in_current_batch_serial} nodes.")
                print(f"  Total processed in run so far: {i + 1}/{len(sorted_target_ids)}")
                print(f"  Time for this batch: {batch_duration_serial:.2f}s ({nodes_per_sec_serial:.2f} nodes/sec)")
                print(f"  .npy files in this batch: {npy_in_current_batch_serial}")
                print(f"  Cumulative .npy files generated: {total_npy_count}")
                print(f"--------------------------------------\n")

                batch_start_time_serial = time.time()
                nodes_in_current_batch_serial = 0
                npy_in_current_batch_serial = 0

    print("\n═══════════ OVERALL PROCESSING COMPLETE ═══════════")
    if args.save_cache and seq_map:
        print(f"\n🔹 Saving {len(seq_map)} sequences to cache: {args.save_cache}...")
        try:
            with open(args.save_cache, 'w') as f:
                json.dump(seq_map, f, indent=2)
            print(f"✔ Sequences saved.")
        except Exception as e:
            sys.stderr.write(f"❌ Error saving cache: {e}\n")
    elif args.save_cache:
        print(f"ℹ️ --save-cache: No sequences to save.")

    total_time = time.time() - overall_start_time
    print(f"\nSummary:")
    print(f"  Attempted: {len(sorted_target_ids)} node(s).")
    print(f"  Successfully processed (data ops completed): {processed_count} node(s).")
    skipped = len(sorted_target_ids) - processed_count
    if skipped > 0: print(f"  Skipped (missing data/errors): {skipped} node(s).")
    print(f"  Output files (summary/npy) for: {successful_output_count} node(s).")
    print(f"  Total .npy files generated: {total_npy_count}.")
    print(f"🏁 Total script time: {total_time:.2f} seconds.")

if __name__ == '__main__':
    main()

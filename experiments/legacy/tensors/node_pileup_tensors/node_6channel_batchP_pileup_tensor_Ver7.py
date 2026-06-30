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
from concurrent.futures import ProcessPoolExecutor, as_completed  # ← added as_completed

# ─────────────────────────────────────────────────────────────────────────────
# Latest format constants

DAT_GLOBAL_MAGIC = b"MYFMT\x01"
GLOBAL_VER_PACK = struct.Struct("<BBI16s")      # major, minor, block_count, reserved[16]

# Latest per-block header: nid, nrec, flags, max_read_len (R), max_cigar_len (C)
BLOCK_HDR_PACK = struct.Struct("<I I H I I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size  # 18

def make_record_struct(R: int, C: int) -> struct.Struct:
    # Per-record layout (latest): <h R s  R s  C s  h c>
    return struct.Struct(f"<h{R}s{R}s{C}shc")

# Latest .idx entry: 30B  <I Q I I H I I>
IDX_ENTRY_SIZE = 30
IDX_ENTRY_PACK = struct.Struct("<I Q I I H I I")  # nid, offset, block_size, n_records, flags, R, C

# NEW: drop INDELs longer than this (bp)
MAX_INDEL_LEN = 50

# ─────────────────────────────────────────────────────────────────────────────
# Encoding tables & tensor constants
BASE_TO_INDEX = {'A': 20, 'C': 30, 'G': 50, 'T': 70,
                 'N': 5, '*': 1, 'I': 90, '_PADDING_': 0}
PADDING_BASE_INDEX = 0

CIGAR_OP_TO_INDEX = {'M': 10, 'N': 20, 'S': 30, 'I': 40, 'D': 50, 'H': 60,
                     'P': 70, '=': 80, 'X': 90, '_PADDING_': 0}
CIGAR_PADDING_INDEX = 0

INDEX_TO_BASE_FOR_VIEW = {20: 'A', 30: 'C', 50: 'G', 70: 'T',
                          5: 'N', 1: '*', 90: 'I', 0: ' '}

TENSOR_WINDOW_SIZE = 100
TENSOR_MAX_READ_ROWS = 200
DEFAULT_QUALITY_PADDING = 0
DEFAULT_MAPPING_QUALITY_PADDING = -1
MISMATCH_CHANNEL_REF_ROW_VALUE = 0
MISMATCH_COMPARISON_PADDING_VALUE = -1

# Worker globals
worker_dat_file = None
worker_base_output_dir = None
worker_need_view = False

# Read-only globals (broadcast to workers; set in parent)
GLOBAL_NODE_SEQS = {}
GLOBAL_NODE_AF = {}

# ─────────────────────────────────────────────────────────────────────────────
# Logging helper

def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers

def read_dat_global_header(dat_path):
    with open(dat_path, "rb") as f:
        magic = f.read(len(DAT_GLOBAL_MAGIC))
        if magic != DAT_GLOBAL_MAGIC:
            raise RuntimeError(f"Invalid .dat magic in {dat_path!r}: {magic!r}")
        major, minor, block_count, _ = GLOBAL_VER_PACK.unpack(f.read(GLOBAL_VER_PACK.size))
        return major, minor, block_count

def load_full_idx_data_latest(idx_path):
    """
    Strictly parse latest .idx format (30 bytes/entry):
      <I Q I I H I I> → nid, offset, block_size, n_records, flags, R, C
    Returns: dict[nid] = (offset, n_records)
    """
    m = {}
    with open(idx_path, "rb") as f:
        size = os.fstat(f.fileno()).st_size
        if size < 4:
            raise RuntimeError("Idx too small")
        (count,) = struct.unpack("<I", f.read(4))
        expected = 4 + count * IDX_ENTRY_SIZE
        if size != expected:
            if (size - 4) // IDX_ENTRY_SIZE != count:
                raise RuntimeError(f"Idx size mismatch: file={size}, count={count}, entry={IDX_ENTRY_SIZE}")
        for i in range(count):
            rec = f.read(IDX_ENTRY_SIZE)
            if len(rec) != IDX_ENTRY_SIZE:
                raise RuntimeError(f"Truncated idx at entry {i+1}")
            nid, offset, _blk_sz, nrec, _flags, _R, _C = IDX_ENTRY_PACK.unpack(rec)
            m[nid] = (offset, nrec)
    return m

def decode_cigar_to_int_ops(cigar_string):
    if not cigar_string or cigar_string == '*':
        return []
    try:
        return [(int(L), op) for L, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception:
        return []

def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]

def af_float_to_bin(x: float) -> int:
    try:
        x = float(x)
    except Exception:
        return 0
    if x <= 0.0: return 0
    if x < 1e-6: return 0
    if x < 1e-5: return 1
    if x < 1e-4: return 2
    if x < 1e-3: return 3
    if x < 1e-2: return 4
    if x < 0.1:  return 5
    if x < 0.5:  return 6
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
            return f"0_D_{deleted}_{v_alt}"
        else:
            anchor_pos = v_pos - 1
            anchor_base = node_seq[anchor_pos].upper() if 0 <= anchor_pos < node_len else "N"
            return f"{anchor_pos}_D_{anchor_base + deleted}_{anchor_base}"

    return f"{v_pos}_{v_type}_{v_ref}_{v_alt}"

def calculate_window_start(variant_pos, window_size):
    return variant_pos - (window_size // 2)

def detect_variants_from_cigar(offset_on_node, cigar_ops_decoded, read_sequence, node_sequence):
    """
    Variant discovery with extra filters:
      - Ignore SNVs if either ref/read base is 'N'
      - Drop indels longer than MAX_INDEL_LEN
      - Ignore indels if inserted/deleted sequence contains 'N'
      - Ignore insertions if the anchor base is 'N'
    """
    variants = []
    node_pos, read_pos = offset_on_node, 0
    node_seq_len, read_seq_len = len(node_sequence), len(read_sequence)

    for L, op in cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            for i in range(L):
                npos, rpos = node_pos + i, read_pos + i
                if npos < node_seq_len and rpos < read_seq_len:
                    nb = node_sequence[npos].upper()
                    rb = read_sequence[rpos].upper()

                    # NEW: ignore SNVs involving N on either side
                    if nb == 'N' or rb == 'N':
                        continue

                    if nb != rb and op != '=':
                        variants.append((npos, 'X', rb, nb))
                else:
                    break
            node_pos += L
            read_pos += L

        elif op == 'I':
            # NEW: drop long insertions
            if L > MAX_INDEL_LEN:
                read_pos += L
                continue

            ins = read_sequence[read_pos: read_pos + L].upper()

            # NEW: ignore insertions containing N (or empty)
            if not ins or 'N' in ins:
                read_pos += L
                continue

            ref_anchor = node_pos - 1 if node_pos > 0 else 0
            anchor_base = node_sequence[ref_anchor].upper() if 0 <= ref_anchor < node_seq_len else "*"

            # NEW: ignore insertions if anchor base is N
            if anchor_base == 'N':
                read_pos += L
                continue

            variants.append((ref_anchor, 'I', ins, anchor_base))
            read_pos += L

        elif op == 'D':
            # NEW: drop long deletions
            if L > MAX_INDEL_LEN:
                node_pos += L
                continue

            del_seq = node_sequence[node_pos: node_pos + L].upper() if node_pos + L <= node_seq_len else ""

            # NEW: ignore deletions containing N (or empty)
            if del_seq and ('N' not in del_seq):
                variants.append((node_pos, 'D', "*", del_seq))

            node_pos += L

        elif op == 'S':
            read_pos += L
        elif op == 'N':
            node_pos += L

    return variants

def get_read_representation_in_window_for_view(cops, offset_on_node, read_seq,
                                               win_start, win_size, node_len):
    out = [' '] * win_size
    npos, rpos = offset_on_node, 0
    rlen = len(read_seq)
    for L, op in cops:
        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = npos + i, rpos + i
                idx = n_aln - win_start
                if 0 <= idx < win_size and r_aln < rlen:
                    out[idx] = read_seq[r_aln].upper()
            npos += L; rpos += L
        elif op in ('D', 'N'):
            for i in range(L):
                idx = (npos + i) - win_start
                if 0 <= idx < win_size:
                    out[idx] = '*'
            npos += L
        elif op in ('I', 'S'):
            rpos += L
        if npos >= win_start + win_size and rpos > 0:
            break
        if rpos >= rlen:
            break
    return out

def get_read_tensor_rows_in_window(cops, offset_on_node, read_seq, qual_vals,
                                   mapq, win_start, win_size, node_len):
    bases = [PADDING_BASE_INDEX] * win_size
    quals = [DEFAULT_QUALITY_PADDING] * win_size
    mapqs = [DEFAULT_MAPPING_QUALITY_PADDING] * win_size
    cig = [CIGAR_PADDING_INDEX] * win_size

    npos, rpos = offset_on_node, 0
    rlen, qlen = len(read_seq), len(qual_vals)

    for L, op in cops:
        opi = CIGAR_OP_TO_INDEX.get(op, CIGAR_PADDING_INDEX)
        if npos >= win_start + win_size and rpos > 0:
            break
        if op in ('M', '=', 'X'):
            for i in range(L):
                n_aln, r_aln = npos + i, rpos + i
                idx = n_aln - win_start
                if 0 <= idx < win_size:
                    cig[idx] = opi
                    mapqs[idx] = mapq
                    if r_aln < rlen:
                        b = read_seq[r_aln].upper()
                        bases[idx] = BASE_TO_INDEX.get(b, BASE_TO_INDEX['N'])
                        if r_aln < qlen:
                            quals[idx] = qual_vals[r_aln]
            npos += L; rpos += L
        elif op in ('D', 'N'):
            for i in range(L):
                idx = (npos + i) - win_start
                if 0 <= idx < win_size:
                    cig[idx] = opi
                    mapqs[idx] = mapq
                    bases[idx] = BASE_TO_INDEX['*']
                    quals[idx] = DEFAULT_QUALITY_PADDING
            npos += L
        elif op in ('I', 'S'):
            rpos += L
        if rpos >= rlen:
            break
    return bases, quals, mapqs, cig

def get_allele_from_read_at_node_pos(read_offset_on_node, read_sequence,
                                     read_quality_values, read_cigar_ops_decoded,
                                     target_node_pos, node_sequence,
                                     expected_var_type=None,
                                     expected_ref_allele_for_indel=None):
    npos = read_offset_on_node
    rpos = 0
    saw_ins_anchor = False
    bq_anchor = None

    # NEW: keep track of the last aligned base-quality before any deletion
    last_aligned_bq = None
    DEFAULT_DEL_ANCHOR_BQ = 50

    for L, op in read_cigar_ops_decoded:
        if op in ('M', '=', 'X'):
            # If target is inside this match block, handle as before
            if npos <= target_node_pos < npos + L:
                off = target_node_pos - npos
                ridx = rpos + off
                if ridx < len(read_sequence):
                    allele = read_sequence[ridx].upper()
                    bq = read_quality_values[ridx] if ridx < len(read_quality_values) else 0
                    if expected_var_type == 'D':
                        return "REF_STATE_FOR_INDEL", bq
                    if expected_var_type == 'I':
                        saw_ins_anchor = True
                        bq_anchor = bq
                    else:
                        return allele, bq
                else:
                    if expected_var_type == 'I':
                        saw_ins_anchor = True
                        bq_anchor = None
                    else:
                        return None, None

            # Update last_aligned_bq with the last base in this aligned stretch
            if L > 0:
                last_idx = rpos + L - 1
                if 0 <= last_idx < len(read_quality_values):
                    last_aligned_bq = read_quality_values[last_idx]
                else:
                    last_aligned_bq = None

            npos += L
            rpos += L

        elif op == 'I':
            if expected_var_type == 'I' and (npos - 1) == target_node_pos:
                if rpos + L <= len(read_sequence):
                    qs = read_quality_values[rpos: rpos + L]
                    mq = sum(qs) / len(qs) if qs else 0.0
                    return read_sequence[rpos: rpos + L].upper(), mq
                return None, None
            rpos += L

        elif op == 'D':
            if npos <= target_node_pos < npos + L:
                # NEW: use anchor base BQ (last aligned base before this deletion), else default 50
                del_anchor_bq = last_aligned_bq if last_aligned_bq is not None else DEFAULT_DEL_ANCHOR_BQ

                if expected_var_type == 'I':
                    return "OTHER_FOR_INDEL", del_anchor_bq

                if expected_var_type == 'D':
                    if 0 <= npos < len(node_sequence) and npos + L <= len(node_sequence):
                        del_ref = node_sequence[npos: npos + L]
                        return ("*" if del_ref == expected_ref_allele_for_indel else "OTHER_FOR_INDEL"), del_anchor_bq
                    return "OTHER_FOR_INDEL", del_anchor_bq

                return "*", del_anchor_bq

            npos += L

        elif op == 'S':
            rpos += L

        elif op == 'N':
            npos += L

        if npos > target_node_pos + 1 and not (expected_var_type == 'I' and (npos - 1) <= target_node_pos):
            break

    if expected_var_type == 'I' and saw_ins_anchor:
        return "REF_STATE_FOR_INDEL", bq_anchor
    return None, None

# ─────────────────────────────────────────────────────────────────────────────
# Worker init & core

def init_worker(dat_file_path_for_worker, base_output_dir_for_worker, need_view_flag):
    global worker_dat_file, worker_base_output_dir, worker_need_view
    worker_dat_file = open(dat_file_path_for_worker, 'rb')
    worker_base_output_dir = base_output_dir_for_worker
    worker_need_view = bool(need_view_flag)

def _non_M_length(cops, variant_type_to_process=None):
    """
    Weighted non-M length used for sorting reads:

    - If variant_type_to_process == 'snp':
        X (mismatch) ops get weight 2.
    - If variant_type_to_process == 'indel':
        I and D ops get weight 2.
    - Otherwise ('all' or None):
        all non-M ops weight 1 (original behavior).
    """
    total = 0
    for L, op in (cops or []):
        if op == 'M':
            continue

        weight = 1
        if variant_type_to_process == 'snp' and op == 'X':
            weight = 2
        elif variant_type_to_process == 'indel' and op in ('I', 'D'):
            weight = 2

        total += L * weight
    return total


def process_single_node_for_pileup(task_args):
    (node_id, dat_file_offset, n_records,
     min_af_threshold, min_variants_threshold, min_allele_bq_threshold,
     variant_type_to_process) = task_args

    global worker_dat_file, worker_need_view
    global GLOBAL_NODE_SEQS, GLOBAL_NODE_AF

    node_sequence = GLOBAL_NODE_SEQS.get(node_id, "")
    if not node_sequence:
        return node_id, {}, [], []  # no tensors, no meta

    genomead_af_list = GLOBAL_NODE_AF.get(node_id, [])
    node_len = len(node_sequence)
    aligned_read_segments = []

    try:
        worker_dat_file.seek(dat_file_offset, os.SEEK_SET)
        hdr = worker_dat_file.read(BLOCK_HDR_SIZE)
        if len(hdr) != BLOCK_HDR_SIZE:
            raise RuntimeError(f"Cannot read block header at {dat_file_offset} for node {node_id}")
        nid2, nrec2, _flags, R, C = BLOCK_HDR_PACK.unpack(hdr)
        if nid2 != node_id or nrec2 != n_records:
            pass
        R = max(1, int(R))
        C = max(1, int(C))
        rec_struct = make_record_struct(R, C)
        rec_size = rec_struct.size

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
                seq = raw_seq.rstrip(b'\0').decode('ascii', 'replace')
                qual_values = list(raw_qual.rstrip(b'\0'))
                cigar_str = raw_cigar.rstrip(b'\0').decode('ascii', 'replace')
                strand_char = strand_byte.decode('ascii') if isinstance(strand_byte, (bytes, bytearray)) else chr(strand_byte)
            except UnicodeDecodeError:
                continue
            if not cigar_str or len(seq) != len(qual_values):
                continue
            cop = decode_cigar_to_int_ops(cigar_str)
            if not cop and cigar_str != '*':
                continue

            rseq = seq
            rqual = qual_values
            rcop = cop
            roff = off_from_file

            if strand_char == '-':
                rseq = reverse_complement(seq)
                rqual = qual_values[::-1]
                rcop = list(reversed(cop)) if cop else []
                span = len(rseq)
                for Lx, opx in cop:
                    if opx == 'I':
                        span -= Lx
                    elif opx == 'D':
                        span += Lx
                roff = node_len - span - off_from_file
                if roff < 0:
                    continue

            aligned_read_segments.append({
                "offset_on_node": roff,
                "read_sequence": rseq,
                "processed_quality_values": rqual,
                "cigar_ops": rcop,
                "original_cigar_str": cigar_str,
                "strand": strand_char,
                "mapping_quality": mapq_val
            })
    except Exception as e:
        sys.stderr.write(f"Error [Node {node_id}]: {e}\n")
        return node_id, None, [], []

    aligned_read_segments.sort(
        key=lambda s: _non_M_length(s["cigar_ops"], variant_type_to_process),
        reverse=True
    )
    if len(aligned_read_segments) > TENSOR_MAX_READ_ROWS:
        aligned_read_segments = aligned_read_segments[: (TENSOR_MAX_READ_ROWS * 2)]

    if not aligned_read_segments:
        return node_id, {}, [], []

    # Variant discovery
    candidate_variants = defaultdict(int)
    for seg in aligned_read_segments:
        for v_pos, v_type, v_alt, v_ref in detect_variants_from_cigar(
                seg["offset_on_node"], seg["cigar_ops"], seg["read_sequence"], node_sequence):
            candidate_variants[(v_pos, v_type, v_ref, v_alt)] += 1

    view_oriented_variant_data = {} if worker_need_view else None
    tensors_for_node = []
    meta_for_node = []

    for (v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar), _ in candidate_variants.items():
        if (variant_type_to_process == 'snp' and v_type != 'X') or \
           (variant_type_to_process == 'indel' and v_type not in ('I', 'D')):
            continue

        alt = ref = other = cov = 0
        alt_bqs = []

        exp_ref, exp_alt = v_ref_from_cigar, v_alt_from_cigar
        ref_for_indel_ctx = None
        if v_type == 'D':
            exp_alt = "*"
            if 0 <= v_pos < node_len:
                exp_ref = node_sequence[v_pos]
            ref_for_indel_ctx = v_ref_from_cigar
        elif v_type == 'I':
            exp_ref = node_sequence[v_pos] if 0 <= v_pos < node_len else "*"
            ref_for_indel_ctx = exp_ref

        for seg in aligned_read_segments:
            allele, bq = get_allele_from_read_at_node_pos(
                seg["offset_on_node"], seg["read_sequence"], seg["processed_quality_values"],
                seg["cigar_ops"], v_pos, node_sequence, v_type, ref_for_indel_ctx)
            if allele is None:
                continue
            cov += 1
            if allele == exp_alt:
                alt += 1
                if bq is not None:
                    alt_bqs.append(bq)
            elif allele == exp_ref or (v_type in ('I', 'D') and allele == "REF_STATE_FOR_INDEL"):
                ref += 1
            else:
                other += 1

        if alt < min_variants_threshold:
            continue
        # if v_type == 'X':
        af_tmp = alt / cov if cov > 0 else 0.0
        if af_tmp < min_af_threshold:
            continue
        bq_tmp = sum(alt_bqs) / len(alt_bqs) if alt_bqs else 0.0
        if bq_tmp < min_allele_bq_threshold:
            continue

        af = alt / cov if cov > 0 else 0.0
        mean_bq = sum(alt_bqs) / len(alt_bqs) if alt_bqs else 0.0

        key = canonical_variant_key(v_pos, v_type, v_ref_from_cigar, v_alt_from_cigar, node_sequence)

        center = v_pos + 1 if v_type == 'I' else v_pos
        win_start = calculate_window_start(center, TENSOR_WINDOW_SIZE)

        # Optional view data (debug)
        if worker_need_view:
            pile = []
            for seg in aligned_read_segments[: TENSOR_MAX_READ_ROWS + 50]:
                row = get_read_representation_in_window_for_view(
                    seg["cigar_ops"], seg["offset_on_node"], seg["read_sequence"],
                    win_start, TENSOR_WINDOW_SIZE, node_len)
                if any(ch != ' ' for ch in row):
                    bases_for_view = [
                        (PADDING_BASE_INDEX if ch == ' ' else BASE_TO_INDEX.get(ch.upper(), BASE_TO_INDEX['N']))
                        for ch in row
                    ]
                    pile.append({
                        "bases": bases_for_view,
                        "offset": seg["offset_on_node"],
                        "strand": seg["strand"],
                        "cigar": seg["original_cigar_str"]
                    })
            view_oriented_variant_data[key] = {
                "pileup_reads_data": pile[:TENSOR_MAX_READ_ROWS],
                "alt_allele_count": alt,
                "ref_allele_count_at_locus": ref,
                "other_allele_count_at_locus": other,
                "coverage_at_locus": cov,
                "alt_allele_frequency": round(af, 4),
                "mean_alt_allele_base_quality": round(mean_bq, 2)
            }

        # Build tensors
        ref_row = [PADDING_BASE_INDEX] * TENSOR_WINDOW_SIZE
        for i, p in enumerate(range(win_start, win_start + TENSOR_WINDOW_SIZE)):
            if 0 <= p < node_len:
                ref_row[i] = BASE_TO_INDEX.get(node_sequence[p].upper(), BASE_TO_INDEX['N'])

        af_row = [0] * TENSOR_WINDOW_SIZE
        g_list = GLOBAL_NODE_AF.get(node_id, [])
        if g_list:
            for i, p in enumerate(range(win_start, win_start + TENSOR_WINDOW_SIZE)):
                if 0 <= p < node_len:
                    v = g_list[p]
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
                    af_row[i] = int(b)

        H, W = 1 + TENSOR_MAX_READ_ROWS, TENSOR_WINDOW_SIZE
        ch1 = np.full((H, W), PADDING_BASE_INDEX, dtype=np.int8)
        ch2 = np.full((H, W), DEFAULT_QUALITY_PADDING, dtype=np.int8)
        ch3 = np.full((H, W), MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
        ch4 = np.full((H, W), DEFAULT_MAPPING_QUALITY_PADDING, dtype=np.int8)
        ch5 = np.full((H, W), CIGAR_PADDING_INDEX, dtype=np.int8)
        ch6 = np.zeros((H, W), dtype=np.int8)

        ch1[0, :] = np.asarray(ref_row, dtype=np.int8)
        ch6[0, :] = np.asarray(af_row, dtype=np.int8)

        reads_added = 0
        variant_win_idx = v_pos - win_start

        for seg in aligned_read_segments:
            if reads_added >= TENSOR_MAX_READ_ROWS:
                break
            mapq = max(0, min(int(seg["mapping_quality"]), 127))
            b_row, q_row, mq_row, cig_row = get_read_tensor_rows_in_window(
                seg["cigar_ops"], seg["offset_on_node"], seg["read_sequence"],
                seg["processed_quality_values"], mapq, win_start,
                TENSOR_WINDOW_SIZE, node_len
            )

            # Encode insertion ALT as 'I' at the anchor with mean BQ
            if v_type == 'I' and 0 <= variant_win_idx < TENSOR_WINDOW_SIZE:
                allele_ins, mean_bq_ins = get_allele_from_read_at_node_pos(
                    seg["offset_on_node"],
                    seg["read_sequence"],
                    seg["processed_quality_values"],
                    seg["cigar_ops"],
                    v_pos,
                    node_sequence,
                    expected_var_type='I',
                    expected_ref_allele_for_indel=(node_sequence[v_pos] if 0 <= v_pos < node_len else "*"),
                )
                if isinstance(allele_ins, str) and allele_ins == v_alt_from_cigar:
                    b_row[variant_win_idx] = BASE_TO_INDEX['I']
                    if mean_bq_ins is not None:
                        q_row[variant_win_idx] = int(max(0, min(127, round(mean_bq_ins))))

            if any(b != PADDING_BASE_INDEX for b in b_row):
                r = 1 + reads_added
                ch1[r, :] = np.asarray(b_row, dtype=np.int8)
                ch2[r, :] = np.asarray(q_row, dtype=np.int8)

                ref_vec = ch1[0, :].astype(np.int16)
                read_vec = ch1[r, :].astype(np.int16)
                flags = np.full(W, MISMATCH_COMPARISON_PADDING_VALUE, dtype=np.int8)
                mask = (read_vec != PADDING_BASE_INDEX) & (ref_vec != PADDING_BASE_INDEX)
                flags[mask] = (read_vec[mask] != ref_vec[mask]).astype(np.int8)
                if 0 <= variant_win_idx < W and mask[variant_win_idx] and flags[variant_win_idx] == 1:
                    flags[variant_win_idx] = 5
                ch3[r, :] = flags

                ch4[r, :] = np.asarray(mq_row, dtype=np.int8)
                ch5[r, :] = np.asarray(cig_row, dtype=np.int8)
                ch6[r, :] = ch6[0, :]

                reads_added += 1

        tensor = np.stack([ch1, ch2, ch3, ch4, ch5, ch6], axis=0)

        meta = {
            "node_id": node_id,
            "variant_key": key,
            "v_pos": int(v_pos),
            "v_type": v_type,
            "v_ref": v_ref_from_cigar,
            "v_alt": v_alt_from_cigar,
            "alt_allele_count": int(alt),
            "ref_allele_count_at_locus": int(ref),
            "other_allele_count_at_locus": int(other),
            "coverage_at_locus": int(cov),
            "alt_allele_frequency": float(round(af, 4)),
            "mean_alt_allele_base_quality": float(round(mean_bq, 2)),
        }

        tensors_for_node.append(tensor)
        meta_for_node.append(meta)

    return node_id, (view_oriented_variant_data or {}), tensors_for_node, meta_for_node

# ─────────────────────────────────────────────────────────────────────────────
# View (unchanged, optional)

def display_pileup_data(node_data_for_display_view, node_id_str_for_display,
                        full_node_sequence, max_reads_to_display_per_variant,
                        max_variants_to_display=float('inf')):
    if max_variants_to_display == 0 or not node_data_for_display_view:
        print(f"Info: No pileup data for node {node_id_str_for_display}.")
        return
    print(f"\n=== Displaying Pileups for Node {node_id_str_for_display} (Len: {len(full_node_sequence)}) ===")
    keys = sorted(node_data_for_display_view.keys(),
                  key=lambda x: (int(x.split('_')[0]), x.split('_')[1]))
    for i, k in enumerate(keys):
        if i >= max_variants_to_display:
            print(f"\n  ... ({len(keys) - i} more variants not shown due to --view limit)")
            break
        v_pos, v_type = int(k.split('_')[0]), k.split('_')[1]
        win_start = calculate_window_start(v_pos + 1 if v_type == 'I' else v_pos, TENSOR_WINDOW_SIZE)
        print(f"\n--- Variant: {k} ---")
        ref = ''.join(full_node_sequence[j] if 0 <= j < len(full_node_sequence) else '0'
                      for j in range(win_start, win_start + TENSOR_WINDOW_SIZE))
        print(f"  Node Ref: {ref}")
        marker_idx = v_pos - win_start
        mark = [' '] * TENSOR_WINDOW_SIZE
        if 0 <= marker_idx < TENSOR_WINDOW_SIZE:
            mark[marker_idx] = '^'
        print(f"  Marker  : {''.join(mark)}")
        pile = node_data_for_display_view[k].get("pileup_reads_data", [])
        for j, row in enumerate(pile[:max_reads_to_display_per_variant]):
            bases_str = ''.join(INDEX_TO_BASE_FOR_VIEW.get(idx, '?') for idx in row["bases"])
            print(f"  Read {j+1:3d}: {bases_str} (CIGAR:{row['cigar']}, STRAND:{row.get('strand', '?')})")
        vd = node_data_for_display_view[k]
        print(f"  Alt Count: {vd.get('alt_allele_count','N/A')}, "
              f"Ref Count: {vd.get('ref_allele_count_at_locus','N/A')}, "
              f"Coverage: {vd.get('coverage_at_locus','N/A')}")
        print(f"  Alt Freq: {vd.get('alt_allele_frequency',0.0):.4f}, "
              f"Mean Alt BQ: {vd.get('mean_alt_allele_base_quality',0.0):.2f}")

# ─────────────────────────────────────────────────────────────────────────────
# Main (Stage 1: tensors only)

def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered tensors from latest .dat/.idx format (R,C per block) into sharded NPYs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("dat", help=".dat alignment file (latest format)")
    parser.add_argument("idx", help=".idx index file (30B entries)")
    parser.add_argument("output", help="Base output directory (will contain shard_XXXXX_data.npy and variant_summary.ndjson)")
    parser.add_argument("candidate_variants_json", help="JSON with nodes and sequences to process.")
    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Worker processes")
    parser.add_argument("--chunksize", type=int, default=64, help="Task chunksize")
    parser.add_argument("--view", nargs='?', const=-1, default=None, type=int, metavar='N',
                        help="Print pileups for top N variants per node (-1 for all)")
    parser.add_argument("--max_view_reads", type=int, default=20, help="Max reads per pileup in view mode")
    parser.add_argument("--min_af", type=float, default=0.05, help="Minimum allele frequency")
    parser.add_argument("--min_variants", type=int, default=3, help="Minimum ALT-supporting reads")
    parser.add_argument("--min_allele_bq", type=float, default=10.0, help="Minimum mean BQ for ALT")
    parser.add_argument("--variant_type", choices=['snp', 'indel', 'all'], default='all',
                        help="Which variants to output tensors for")
    parser.add_argument("--shard_size", type=int, default=4096,
                        help="Number of variant tensors per full shard")
    parser.add_argument("--log_every_nodes", type=int, default=100000,
                        help="Log progress every N nodes")
    # NEW: customizable shard prefix
    parser.add_argument(
        "--shard_prefix",
        type=str,
        default="shard",
        help="Prefix for shard filenames (prefix_00000_data.npy)"
    )

    args = parser.parse_args()

    if not all(os.path.isfile(p) for p in (args.dat, args.idx, args.candidate_variants_json)):
        sys.exit("Error: dat/idx/json not found.")

    # Validate latest .dat header
    try:
        major, minor, blocks = read_dat_global_header(args.dat)
        log(f".dat header: version {major}.{minor}, blocks={blocks} (expect latest format)")
    except Exception as e:
        sys.exit(f"Error reading .dat header: {e}")

    # Load idx (latest only)
    try:
        idx_map = load_full_idx_data_latest(args.idx)
    except Exception as e:
        sys.exit(f"Error reading .idx (latest): {e}")

    # Input nodes JSON
    node_sequences, node_af_data, node_ids = {}, {}, set()
    try:
        with open(args.candidate_variants_json, 'r') as f:
            data = json.load(f)
            for node in data.get('nodes', []):
                nid_str = node.get('node_id')
                seq = node.get('sequence')
                af = node.get('genomead_af', [])
                if nid_str and seq:
                    try:
                        nid = int(nid_str)
                        node_sequences[nid] = seq.upper()
                        node_af_data[nid] = af
                        node_ids.add(nid)
                    except ValueError:
                        pass
    except Exception as e:
        sys.exit(f"Error parsing JSON: {e}")

    log(f"Loaded {len(node_sequences)} candidate nodes from JSON.")

    tasks = []
    missing = 0
    for nid in node_ids:
        if nid in idx_map:
            off, nrec = idx_map[nid]
            tasks.append((nid, off, nrec,
                          args.min_af, args.min_variants,
                          args.min_allele_bq, args.variant_type))
        else:
            missing += 1
    if missing:
        log(f"Warning: {missing} nodes from JSON not found in idx; skipped.")
    if not tasks:
        sys.exit("No tasks to run after JSON/idx filtering.")

    tasks.sort(key=lambda t: t[1])  # by offset
    os.makedirs(args.output, exist_ok=True)

    global GLOBAL_NODE_SEQS, GLOBAL_NODE_AF
    GLOBAL_NODE_SEQS, GLOBAL_NODE_AF = node_sequences, node_af_data

    need_view = (args.view is not None and args.view != 0)
    total_nodes = len(tasks)
    log(f"Submitting {total_nodes} nodes to {args.num_workers} workers...")

    shard_size = max(1, int(args.shard_size))
    shard_idx = 0
    buffer_tensors = []
    buffer_meta = []
    total_variants = 0

    summary_path = os.path.join(args.output, "variant_summary.ndjson")
    summary_f = open(summary_path, "w")

    t0 = time.time()
    processed_nodes = 0
    batch_nodes = 0
    batch_variants = 0

    def flush_full_shard():
        """
        Take exactly `shard_size` tensors/meta from the front of the buffers,
        write them as one shard, and keep any leftovers in the buffers.
        """
        nonlocal shard_idx, buffer_tensors, buffer_meta, total_variants

        if len(buffer_tensors) < shard_size:
            return

        # take exactly shard_size items
        chunk_tensors = buffer_tensors[:shard_size]
        chunk_meta = buffer_meta[:shard_size]

        # keep leftovers in buffer
        buffer_tensors = buffer_tensors[shard_size:]
        buffer_meta = buffer_meta[shard_size:]

        xs = np.stack(chunk_tensors, axis=0)  # (shard_size, 6, 201, 100)
        N = xs.shape[0]
        data_path = os.path.join(args.output, f"{args.shard_prefix}_{shard_idx:05d}_data.npy")
        log(f"Saving shard {shard_idx} data: {N} tensors -> {data_path}")
        np.save(data_path, xs)

        # meta for this shard
        for i, m in enumerate(chunk_meta):
            m_out = dict(m)
            m_out["shard_index"] = shard_idx
            m_out["index_within_shard"] = i
            summary_f.write(json.dumps(m_out) + "\n")

        total_variants += N
        shard_idx += 1

    def flush_remainder():
        """
        Write any remaining tensors as a final shard (may be smaller than shard_size).
        """
        nonlocal shard_idx, buffer_tensors, buffer_meta, total_variants

        if not buffer_tensors:
            return

        xs = np.stack(buffer_tensors, axis=0)
        N = xs.shape[0]
        data_path = os.path.join(args.output, f"{args.shard_prefix}_{shard_idx:05d}_data.npy")
        log(f"Saving FINAL shard {shard_idx} data: {N} tensors -> {data_path}")
        np.save(data_path, xs)

        for i, m in enumerate(buffer_meta):
            m_out = dict(m)
            m_out["shard_index"] = shard_idx
            m_out["index_within_shard"] = i
            summary_f.write(json.dumps(m_out) + "\n")

        total_variants += N
        shard_idx += 1

        buffer_tensors.clear()
        buffer_meta.clear()

    # ───── revised ex.map(...) → bounded submit + as_completed ────
    max_tasks_in_flight = max(1, args.num_workers * 4)
    task_iter = iter(tasks)

    with ProcessPoolExecutor(max_workers=args.num_workers,
                             initializer=init_worker,
                             initargs=(args.dat, args.output, need_view)) as ex:
        futures = {}

        # Seed initial tasks
        try:
            for _ in range(max_tasks_in_flight):
                t = next(task_iter)
                f = ex.submit(process_single_node_for_pileup, t)
                futures[f] = True
        except StopIteration:
            pass

        try:
            while futures:
                for f in as_completed(list(futures.keys())):
                    del futures[f]
                    nid, view_data, node_tensors, node_meta = f.result()

                    processed_nodes += 1
                    batch_nodes += 1

                    if node_tensors:
                        buffer_tensors.extend(node_tensors)
                        buffer_meta.extend(node_meta)
                        n_new = len(node_tensors)
                        batch_variants += n_new

                    # Optional live view print
                    if need_view and view_data:
                        seq_for_view = GLOBAL_NODE_SEQS.get(nid, "")
                        display_pileup_data(view_data, str(nid), seq_for_view,
                                            args.max_view_reads,
                                            args.view if args.view != -1 else float('inf'))

                    # Flush any full shards from the buffer
                    while len(buffer_tensors) >= shard_size:
                        flush_full_shard()

                    # Progress logging (unchanged logic)
                    if (processed_nodes % args.log_every_nodes == 0) or (processed_nodes == total_nodes):
                        dt = time.time() - t0
                        rate_nodes = batch_nodes / dt if dt > 0 else 0.0
                        rate_vars = batch_variants / dt if dt > 0 else 0.0
                        log(f"Batch: +{batch_nodes} nodes (+{batch_variants} variants) "
                            f"in {dt:.2f}s → {rate_nodes:.1f} nodes/s, {rate_vars:.1f} vars/s. "
                            f"[total {processed_nodes}/{total_nodes} nodes, "
                            f"{total_variants} variants, {shard_idx} shards written so far]")
                        batch_nodes = 0
                        batch_variants = 0
                        t0 = time.time()

                    # Refill the queue up to max_tasks_in_flight
                    try:
                        while len(futures) < max_tasks_in_flight:
                            t = next(task_iter)
                            nf = ex.submit(process_single_node_for_pileup, t)
                            futures[nf] = True
                    except StopIteration:
                        pass

                    if not futures:
                        break
        except KeyboardInterrupt:
            log("KeyboardInterrupt detected, shutting down executor...")
            raise

    # Final shard flush for leftovers (may be < shard_size)
    flush_remainder()
    summary_f.close()

    log("Done.")
    log(f"Total nodes processed   : {processed_nodes}")
    log(f"Total variant tensors   : {total_variants}")
    log(f"Total shards written    : {shard_idx}")
    log(f"Variant summary (NDJSON): {summary_path}")

if __name__ == "__main__":
    main()

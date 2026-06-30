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
import gzip
from typing import Any, Dict, Iterable, Iterator, List, Optional, Set, Tuple

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
MISMATCH_COMPARISON_PADDING_VALUE = -1

# Worker globals
worker_dat_file = None
worker_base_output_dir = None
worker_need_view = False

# Read-only globals (broadcast to workers; set in parent)
GLOBAL_NODE_SEQS: Dict[int, str] = {}
GLOBAL_NODE_AF: Dict[int, Any] = {}

# ─────────────────────────────────────────────────────────────────────────────
# Logging helper

def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers

def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def read_dat_global_header(dat_path):
    with open(dat_path, "rb") as f:
        magic = f.read(len(DAT_GLOBAL_MAGIC))
        if magic != DAT_GLOBAL_MAGIC:
            raise RuntimeError(f"Invalid .dat magic in {dat_path!r}: {magic!r}")
        major, minor, block_count, _ = GLOBAL_VER_PACK.unpack(f.read(GLOBAL_VER_PACK.size))
        return major, minor, block_count

def iter_idx_latest(idx_path: str) -> Iterator[Tuple[int, int, int]]:
    """
    Stream latest .idx entries (30 bytes/entry) without loading a huge dict.
    Yields (nid, offset, n_records) in file order.
    """
    with open(idx_path, "rb") as f:
        hdr = f.read(4)
        if len(hdr) != 4:
            raise RuntimeError("idx header too short")
        (count,) = struct.unpack("<I", hdr)

        # sanity: size may be exact; tolerate some mismatch but still stream safely
        f.seek(0, os.SEEK_END)
        size = f.tell()
        expected = 4 + count * IDX_ENTRY_SIZE
        if size != expected:
            # Still allow streaming if it looks consistent enough.
            # (Some older generators may have padding; we won’t hard-fail.)
            if size < 4 + IDX_ENTRY_SIZE:
                raise RuntimeError(f"idx size too small: file={size}, expected_at_least={4 + IDX_ENTRY_SIZE}")
        f.seek(4, os.SEEK_SET)

        for i in range(count):
            rec = f.read(IDX_ENTRY_SIZE)
            if len(rec) != IDX_ENTRY_SIZE:
                # truncated tail
                break
            nid, off, _blk_sz, nrec, _flags, _R, _C = IDX_ENTRY_PACK.unpack(rec)
            yield int(nid), int(off), int(nrec)

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
                    if nb != rb and op != '=':
                        variants.append((npos, 'X', rb, nb))
                else:
                    break
            node_pos += L
            read_pos += L
        elif op == 'I':
            ins = read_sequence[read_pos: read_pos + L].upper()
            ref_anchor = node_pos - 1 if node_pos > 0 else 0
            ref_base = node_sequence[ref_anchor].upper() if 0 <= ref_anchor < node_seq_len else "*"
            variants.append((ref_anchor, 'I', ins, ref_base))
            read_pos += L
        elif op == 'D':
            del_seq = node_sequence[node_pos: node_pos + L].upper() if node_pos + L <= node_seq_len else ""
            if del_seq:
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

    # keep track of last aligned base-quality before any deletion
    last_aligned_bq = None
    DEFAULT_DEL_ANCHOR_BQ = 50

    for L, op in read_cigar_ops_decoded:
        if op in ('M', '=', 'X'):
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

    node_len = len(node_sequence)
    aligned_read_segments = []

    try:
        worker_dat_file.seek(dat_file_offset, os.SEEK_SET)
        hdr = worker_dat_file.read(BLOCK_HDR_SIZE)
        if len(hdr) != BLOCK_HDR_SIZE:
            raise RuntimeError(f"Cannot read block header at {dat_file_offset} for node {node_id}")
        nid2, nrec2, _flags, R, C = BLOCK_HDR_PACK.unpack(hdr)
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
# JSON streaming for candidate nodes (NEW)
# Supports:
#   1) top-level list: [ {...}, {...}, ... ]
#   2) wrapper dict:  {"nodes":[ {...}, ... ], ...}
#
# This is memory-safe and works for huge JSON (even gz).
# ─────────────────────────────────────────────────────────────────────────────

def iter_node_objects_from_candidates(path: str) -> Iterator[Dict[str, Any]]:
    """
    Yield node dict objects from:
      - a top-level JSON list, OR
      - a top-level dict containing "nodes": [ ... ]
    We do a lightweight bracket/object brace-tracking parse and json.loads per object.
    """
    with open_maybe_gzip(path, "rt") as f:
        # We need to reach the array that contains node objects.
        # If file starts with '{', we find the first '[' after the "nodes" key (best-effort).
        # Otherwise, if file starts with '[', we use that array directly.

        # Read until we hit either '[' (array start) or EOF.
        ch = f.read(1)
        while ch and ch.isspace():
            ch = f.read(1)
        if not ch:
            return

        if ch == '[':
            # top-level list
            pass
        elif ch == '{':
            # wrapper dict; scan until we find the first '[' which should open the nodes array
            # (this assumes "nodes" array appears before any other big arrays; typical in your files)
            in_str = False
            esc = False
            prev = ch
            while True:
                ch = f.read(1)
                if not ch:
                    return
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
                    elif ch == '[':
                        break
            # now we're positioned just after '[' (array open)
        else:
            # unexpected top-level token
            return

        # Now extract each {...} object from the array.
        buf: List[str] = []
        depth = 0
        in_str = False
        esc = False

        while True:
            ch = f.read(1)
            if not ch:
                break

            if depth == 0:
                if ch.isspace() or ch == ',':
                    continue
                if ch == ']':
                    break
                if ch != '{':
                    continue
                buf = ['{']
                depth = 1
                in_str = False
                esc = False
                continue

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
                        obj_str = ''.join(buf)
                        buf = []
                        try:
                            obj = json.loads(obj_str)
                            if isinstance(obj, dict):
                                yield obj
                        except Exception:
                            # skip malformed object
                            pass

# ─────────────────────────────────────────────────────────────────────────────
# Main (Stage 1: tensors only) — revised to accept JSON LIST output from your
# "build_nodes_from_idx.py", and to be memory-safe on huge files.
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate variant-centered tensors from latest .dat/.idx format (R,C per block) into sharded NPYs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("dat", help=".dat alignment file (latest format)")
    parser.add_argument("idx", help=".idx index file (30B entries)")
    parser.add_argument("output", help="Base output directory (will contain prefix_XXXXX_data.npy and variant_summary.ndjson)")
    parser.add_argument("candidate_variants_json",
                        help="Candidate nodes JSON: can be {'nodes':[...]} OR a JSON LIST of node records. Supports .gz")
    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Worker processes")
    parser.add_argument("--chunksize", type=int, default=64, help="Task chunksize (not used directly; kept for compatibility)")
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
    parser.add_argument("--shard_prefix", type=str, default="shard",
                        help="Prefix for shard filenames (prefix_00000_data.npy)")

    # NEW: memory-safe wave size (how many node sequences to keep in RAM per worker pool run)
    parser.add_argument("--wave_size", type=int, default=200000,
                        help="Nodes per wave (controls memory use when candidate JSON is huge)")

    args = parser.parse_args()

    if not all(os.path.isfile(p) for p in (args.dat, args.idx, args.candidate_variants_json)):
        sys.exit("Error: dat/idx/json not found.")

    # Validate latest .dat header
    try:
        major, minor, blocks = read_dat_global_header(args.dat)
        log(f".dat header: version {major}.{minor}, blocks={blocks} (expect latest format)")
    except Exception as e:
        sys.exit(f"Error reading .dat header: {e}")

    os.makedirs(args.output, exist_ok=True)

    need_view = (args.view is not None and args.view != 0)
    shard_size = max(1, int(args.shard_size))

    # Output summary file
    summary_path = os.path.join(args.output, "variant_summary.ndjson")
    summary_f = open(summary_path, "w")

    # Shard buffers (persist across waves)
    shard_idx = 0
    buffer_tensors: List[np.ndarray] = []
    buffer_meta: List[Dict[str, Any]] = []
    total_variants = 0

    def flush_full_shard():
        nonlocal shard_idx, buffer_tensors, buffer_meta, total_variants
        while len(buffer_tensors) >= shard_size:
            chunk_tensors = buffer_tensors[:shard_size]
            chunk_meta = buffer_meta[:shard_size]
            buffer_tensors = buffer_tensors[shard_size:]
            buffer_meta = buffer_meta[shard_size:]

            xs = np.stack(chunk_tensors, axis=0)  # (N, 6, 201, 100)
            N = xs.shape[0]
            data_path = os.path.join(args.output, f"{args.shard_prefix}_{shard_idx:05d}_data.npy")
            log(f"Saving shard {shard_idx} data: {N} tensors -> {data_path}")
            np.save(data_path, xs)

            for i, m in enumerate(chunk_meta):
                m_out = dict(m)
                m_out["shard_index"] = shard_idx
                m_out["index_within_shard"] = i
                summary_f.write(json.dumps(m_out) + "\n")

            total_variants += N
            shard_idx += 1

    def flush_remainder():
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

    # ─────────────────────────────────────────────────────────────────────────
    # Core revision:
    # We stream idx sequentially, and stream the candidate JSON sequentially too.
    # If your candidate JSON was produced by your "build_nodes_from_idx.py", it
    # is in the SAME order as idx, so we can process in lockstep safely.
    #
    # We also process in waves to avoid holding all sequences for 15M nodes.
    # ─────────────────────────────────────────────────────────────────────────
    wave_size = max(1, int(args.wave_size))

    idx_iter = iter_idx_latest(args.idx)
    json_iter = iter_node_objects_from_candidates(args.candidate_variants_json)

    processed_nodes = 0
    submitted_nodes_total = 0

    t0_global = time.time()
    t_batch = time.time()
    batch_nodes = 0
    batch_variants = 0

    wave_num = 0

    while True:
        wave_num += 1

        # Build one wave: tasks + per-wave seq/af maps
        wave_tasks: List[Tuple[int, int, int, float, int, float, str]] = []
        wave_seqs: Dict[int, str] = {}
        wave_af: Dict[int, Any] = {}

        # Fill up to wave_size nodes, in idx/json lockstep
        got_any = False
        for _ in range(wave_size):
            try:
                nid, off, nrec = next(idx_iter)
            except StopIteration:
                break

            try:
                node_obj = next(json_iter)
            except StopIteration:
                # JSON ended early; stop processing.
                log("WARNING: candidate JSON ended before idx ended. Stopping.")
                break

            got_any = True

            # We expect matching node_id if JSON was generated from this idx.
            # If not present/mismatched, we still attempt to use the JSON node_id.
            nid_json_raw = node_obj.get("node_id")
            try:
                nid_json = int(str(nid_json_raw)) if nid_json_raw is not None else None
            except Exception:
                nid_json = None

            if nid_json is not None and nid_json != nid:
                # mismatch: keep going but warn occasionally
                if submitted_nodes_total < 5:
                    log(f"WARNING: idx/json node_id mismatch: idx={nid}, json={nid_json}. "
                        f"Continuing using idx offsets for idx node_id.")
                # We still use idx nid for alignment offsets and for indexing globals.
                # Sequence will come from this JSON object if it has it.
                # (This keeps behavior predictable.)

            seq = node_obj.get("sequence", "")
            if not isinstance(seq, str) or not seq:
                continue  # skip nodes without sequence

            # Optional genomead_af
            af = node_obj.get("genomead_af", [])

            wave_seqs[nid] = seq.upper()
            wave_af[nid] = af

            wave_tasks.append((nid, off, nrec,
                               args.min_af, args.min_variants,
                               args.min_allele_bq, args.variant_type))

        if not got_any:
            break
        if not wave_tasks:
            # No nodes with sequences in this wave
            continue

        # Publish wave-scoped globals BEFORE forking
        GLOBAL_NODE_SEQS.clear()
        GLOBAL_NODE_AF.clear()
        GLOBAL_NODE_SEQS.update(wave_seqs)
        GLOBAL_NODE_AF.update(wave_af)

        submitted_nodes_total += len(wave_tasks)
        log(f"[Wave {wave_num}] Submitting {len(wave_tasks)} nodes to {args.num_workers} workers "
            f"(wave_size={wave_size})")

        # Run this wave
        max_tasks_in_flight = max(1, args.num_workers * 2)
        task_iter = iter(wave_tasks)

        with ProcessPoolExecutor(max_workers=args.num_workers,
                                 initializer=init_worker,
                                 initargs=(args.dat, args.output, need_view)) as ex:
            futures = {}
            try:
                for _ in range(max_tasks_in_flight):
                    t = next(task_iter)
                    f = ex.submit(process_single_node_for_pileup, t)
                    futures[f] = True
            except StopIteration:
                pass

            while futures:
                for f in as_completed(list(futures.keys())):
                    del futures[f]
                    nid_done, view_data, node_tensors, node_meta = f.result()

                    processed_nodes += 1
                    batch_nodes += 1

                    if node_tensors:
                        buffer_tensors.extend(node_tensors)
                        buffer_meta.extend(node_meta)
                        n_new = len(node_tensors)
                        batch_variants += n_new

                    # Optional view
                    if need_view and view_data:
                        seq_for_view = GLOBAL_NODE_SEQS.get(nid_done, "")
                        display_pileup_data(view_data, str(nid_done), seq_for_view,
                                            args.max_view_reads,
                                            args.view if args.view != -1 else float('inf'))

                    # Flush full shards
                    flush_full_shard()

                    # Progress logging
                    if (processed_nodes % args.log_every_nodes == 0):
                        dt = time.time() - t_batch
                        rate_nodes = batch_nodes / dt if dt > 0 else 0.0
                        rate_vars = batch_variants / dt if dt > 0 else 0.0
                        log(f"Batch: +{batch_nodes} nodes (+{batch_variants} variants) "
                            f"in {dt:.2f}s → {rate_nodes:.1f} nodes/s, {rate_vars:.1f} vars/s. "
                            f"[total {processed_nodes} nodes, {total_variants} variants, "
                            f"{shard_idx} shards written so far]")
                        batch_nodes = 0
                        batch_variants = 0
                        t_batch = time.time()

                    # Refill
                    try:
                        while len(futures) < max_tasks_in_flight:
                            t = next(task_iter)
                            nf = ex.submit(process_single_node_for_pileup, t)
                            futures[nf] = True
                    except StopIteration:
                        pass

                    if not futures:
                        break

        # After wave completes, collect garbage-friendly
        GLOBAL_NODE_SEQS.clear()
        GLOBAL_NODE_AF.clear()

    # Final flush
    flush_remainder()
    summary_f.close()

    dt_all = time.time() - t0_global
    log("Done.")
    log(f"Total nodes processed   : {processed_nodes}")
    log(f"Total variant tensors   : {total_variants}")
    log(f"Total shards written    : {shard_idx}")
    log(f"Variant summary (NDJSON): {summary_path}")
    log(f"Wall time              : {dt_all/3600:.2f} hours")

if __name__ == "__main__":
    main()

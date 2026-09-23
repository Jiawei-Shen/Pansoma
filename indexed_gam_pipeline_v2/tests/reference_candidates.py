"""Frozen reference: support counting, windowing and tensor encoding as they were before the
VisitView/NodeReads speedup (v5 production run 363651, source sha256 recorded in its config.json).

Used only by test_reference_equivalence.py to prove the current implementation is output-identical.
Do not edit: these are verbatim copies of the previous candidates.overlap/anchor_window/make_tensor.
"""
from collections import Counter

import numpy as np

from indexed_gam_pipeline_v2.candidates import (BASES, CHANNELS, FORMAT_VERSION, OPS, ROW_ORDER, ROW_SELECTION_VERSION,
    STORAGE_VERSION, STRAND, WINDOW_ENCODING_VERSION, Column, boundary_evidence, candidate_alt_codes, encode_count,
    int8_quality, oriented_columns, split_columns)


def overlap(read, candidate, min_bq):
    """Classify one record as ('alt'|'ref'|'other', anchor visit), or None if it does not cover.

    A record counts once even when it visits the node repeatedly: ALT beats REF beats
    other, and the earliest mapping breaks ties. ALT needs an exact observation with
    quality >= min_bq. REF needs matching graph bases over the whole interval; for an
    insertion, REF needs adjacent M/X columns on both sides of the boundary.
    """
    choices = []
    for visit in read.visits:
        if visit.node != candidate.node or visit.first == visit.last:
            continue
        if candidate.kind == "INS":
            covered = visit.start <= candidate.start <= visit.end
        else:
            covered = visit.start < candidate.end and visit.end > candidate.start
        if not covered:
            continue
        exact = any(o.visit == visit.index and o.candidate == candidate and o.quality >= min_bq
                    for o in read.observations)
        cols = oriented_columns(read, visit, 1)
        support = "other"
        if exact:
            support = "alt"
        elif candidate.kind == "INS":
            left, middle, right = split_columns(cols, visit, candidate)
            if not middle and boundary_evidence(left, right, visit, candidate.start):
                support = "ref"
        else:
            local = [c for c in cols if c.visit == visit.index]
            affected = [c for c in local if not c.boundary and candidate.start <= c.pos < candidate.end]
            inserted = any(c.boundary and candidate.start < c.pos < candidate.end for c in local)
            if (len(affected) == len(candidate.ref) and not inserted
                    and all(c.op in ("M", "X") and c.read == c.ref for c in affected)):
                support = "ref"
        choices.append(({"alt": 0, "ref": 1, "other": 2}[support], visit.index, support, visit))
    if not choices:
        return None
    _, _, support, visit = min(choices, key=lambda x: x[:2])
    return support, visit


def anchor_window(read, visit, candidate, width):
    """Crop this record's columns so the candidate starts at column width // 2.

    Returns (row, cropped): `row` has exactly `width` entries (Column or None for
    missing evidence); `cropped` counts central columns lost past the right edge.
    INS rows without the insertion reserve gap ("G") slots only when both flanks
    are present; DEL rows keep one "D" column per deleted graph base.
    """
    anchor = width // 2
    cols = oriented_columns(read, visit, width)
    left, center, right = split_columns(cols, visit, candidate)
    central_length = len(center)
    if candidate.kind == "INS":
        evidence = bool(center) or boundary_evidence(left, right, visit, candidate.start)
        slots = min(len(candidate.alt), width - anchor)
        if len(center) < slots:
            gap = Column("-", "-", -1, "G", candidate.node, candidate.start, False, visit.index, True)
            center = center + [gap if evidence else None] * (slots - len(center))
        central_length = max(central_length, slots)
        missing = 0
    else:
        # Coverage starting after the anchor keeps its true offset: pad the missing start.
        covered = [c.pos for c in center if c is not None and not c.boundary]
        missing = min(width, max(0, min(covered) - candidate.start)) if covered else 0
    row = [None] * max(0, anchor - len(left)) + left[-anchor:] if anchor else []
    tail = [None] * missing + center + right
    row += tail[:width - anchor]
    row += [None] * (width - len(row))
    cropped = max(0, missing + central_length - (width - anchor))
    return row, cropped


def make_tensor(candidate, eligible, path_counts, rows=200, width=101, debug=False):
    """Encode one candidate from `eligible` = [(read, support, anchor visit)].

    Every eligible record is windowed and ranked (descending visible edit bp, then
    MAPQ, record hash, mapping index), grouped stably by its visible node path, and
    only then uniformly sampled down to `rows`. Coverage/AF use all eligible records.
    """
    if rows < 1 or width < 1:
        raise ValueError("Tensor rows and width must be positive")
    if any(not isinstance(c, (int, np.integer)) or not 0 <= c <= np.iinfo(np.int32).max
           for c in path_counts.values()):
        raise ValueError("Node path counts must be nonnegative int32 integers")
    start = width // 2
    # The candidate region is the same columns in every row: [start, start + allele length).
    alt_stripe = {start + k: code for k, code in enumerate(candidate_alt_codes(candidate)) if start + k < width}
    tensor = np.zeros((len(CHANNELS), rows, width), dtype=np.int8)
    windows = []
    for read, support, visit in eligible:
        row, overflow = anchor_window(read, visit, candidate, width)
        # Group key: the node path actually visible in this window (mapping transitions
        # only; offsets and insertion lengths are ignored).
        path, previous_visit = [], None
        for col in row:
            if col is not None and col.visit != previous_visit:
                path.append((col.node, col.reverse))
                previous_visit = col.visit
        mismatch = sum(c is not None and c.op != "G" and c.read != c.ref for c in row)
        windows.append((read, support, visit, row, tuple(path), mismatch, overflow))
    windows.sort(key=lambda w: (-w[5], -w[0].mapq, w[0].digest, w[2].index))
    path_groups = {}
    for window in windows:
        path_groups.setdefault(window[4], []).append(window)
    ordered = [window for group in path_groups.values() for window in group]
    n_selected = min(rows, len(ordered))
    if n_selected == 1:
        sample_indices = [len(ordered) // 2]
    elif n_selected:
        sample_indices = [i * (len(ordered) - 1) // (n_selected - 1) for i in range(n_selected)]
    else:
        sample_indices = []
    chosen = [ordered[i] for i in sample_indices]
    groups, omitted, details = [], [], []
    for ri, (read, support, visit, row, path_key, mismatch, overflow) in enumerate(chosen):
        path = [dict(node_id=node, reverse=reverse) for node, reverse in path_key]
        if not groups or groups[-1]["path"] != path:
            groups.append(dict(start_row=ri, end_row=ri + 1, path=path))
        else:
            groups[-1]["end_row"] = ri + 1
        strand = STRAND["reverse"] if visit.reverse else STRAND["forward"]
        for ci, col in enumerate(row):
            if col is not None:
                tensor[:, ri, ci] = [BASES.get(col.read, 5), int8_quality(col.quality), alt_stripe.get(ci, 0),
                                     int8_quality(read.mapq), OPS[col.op], BASES.get(col.ref, 5),
                                     encode_count(path_counts[col.node]), strand]
        if overflow:
            omitted.append(dict(row_index=ri, omitted_columns=overflow,
                                reason="candidate_region_extends_beyond_anchor_window"))
        if debug:
            details.append(dict(read_name=read.name, record_sha256=read.digest, support=support,
                anchor_mapping_index=visit.index, reversed_for_candidate=visit.reverse,
                window_mismatch_bp=mismatch, grouped_rank=sample_indices[ri],
                omitted_central_context_columns=overflow,
                path=[dict(node_id=v.node, start=v.start, end=v.end, reverse=v.reverse,
                           mapping_index=v.index) for v in read.visits],
                columns=[c.graph() if c is not None else None for c in row]))
    counts = Counter(s for _, s, _ in eligible)
    selected_counts = Counter(w[1] for w in chosen)
    candidate_span = max(1, len(candidate.alt) if candidate.kind == "INS" else len(candidate.ref))
    meta = dict(candidate.metadata(), tensor_storage_version=STORAGE_VERSION,
        tensor_format_version=FORMAT_VERSION, row_order=ROW_ORDER, row_groups=groups,
        row_selection_version=ROW_SELECTION_VERSION, selected_grouped_ranks=sample_indices,
        window_mismatch_bp=[w[5] for w in chosen], window_encoding_version=WINDOW_ENCODING_VERSION,
        anchor_column=start, candidate_columns=[start, min(width, start + candidate_span)],
        coverage=len(eligible), alt_count=counts["alt"], ref_count=counts["ref"],
        other_count=counts["other"], af=counts["alt"] / len(eligible) if eligible else 0,
        selected_alignments=len(chosen), omitted_context=omitted, selected_counts=dict(selected_counts))
    if debug:
        meta["rows"] = details
        meta["selection_audit"] = [dict(record_sha256=w[0].digest, mapping_quality=w[0].mapq,
            anchor_mapping_index=w[2].index, window_mismatch_bp=w[5], support=w[1],
            path=[list(p) for p in w[4]]) for w in windows]
    return tensor, meta

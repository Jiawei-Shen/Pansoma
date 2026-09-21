"""Candidate-centered GAM edit model. Coordinates are zero-based forward-node coordinates.

Identity normalization is exact node/forward interval/alleles (including reverse
complements), not repeat left-alignment or equivalence across graph nodes.
"""
from collections import Counter
from dataclasses import dataclass
import hashlib
from types import SimpleNamespace

import numpy as np

VERSION = "indexed-gam-candidate-v2"
V3_VERSION = "indexed-gam-candidate-v3"
ROW_SELECTION_VERSION = "window-edit-bp-group-uniform-v1"
ROW_ORDER = "descending visible-window substitution+insertion+deletion bp; stable visible-node-path groups in first-occurrence order; uniform ordered sampling"
CHANNELS = ["read_base", "base_quality", "event_flags", "mapping_quality",
            "alignment_operation", "row_graph_reference_base"]
V3_CHANNELS = CHANNELS + ["node_distinct_w_record_count"]
BASES = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5, "-": 6}
OPS = {"M": 1, "X": 2, "I": 3, "D": 4, "C": 5, "G": 6}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


@dataclass(frozen=True, order=True)
class Candidate:
    node: int
    start: int
    ref: str
    alt: str
    kind: str

    @property
    def end(self):
        return self.start + len(self.ref)

    def metadata(self):
        return dict(candidate_id=f"{self.node}:{self.start}:{self.kind}:{self.ref}>{self.alt}",
                    node_id=self.node, start=self.start, end=self.end, orientation="+",
                    ref=self.ref, alt=self.alt, event_type=self.kind,
                    event_length=max(len(self.ref), len(self.alt)))


@dataclass
class Column:
    __slots__ = ("read", "ref", "quality", "op", "node", "pos", "reverse", "visit", "boundary")

    read: str
    ref: str
    quality: int
    op: str
    node: int
    pos: int
    reverse: bool
    visit: int
    boundary: bool

    def flipped(self):
        return Column(rc(self.read), rc(self.ref), self.quality, self.op,
                      self.node, self.pos, not self.reverse, self.visit, self.boundary)

    def graph(self):
        return dict(node_id=self.node, **({"boundary": self.pos} if self.boundary
                    else {"offset": self.pos}), reverse=self.reverse,
                    mapping_index=self.visit)


@dataclass
class Visit:
    node: int
    start: int
    end: int
    reverse: bool
    index: int
    first: int
    last: int
    node_length: int


@dataclass
class Observation:
    candidate: Candidate
    visit: int
    quality: float


@dataclass
class Read:
    name: str
    digest: str
    mapq: int
    columns: list
    visits: list
    observations: list


def decode_alignment(alignment, sequences, max_indel=50):
    """Consume original edits, retaining every mapping and both cursors.

    Complex replacements use explicitly marked C slots (ordered zip, gaps at
    the tail), solely for context; they never become supported candidates.
    """
    columns, visits, observations, unsupported = [], [], [], []
    read_cursor = 0
    sequence = alignment.sequence.upper()
    qualities = alignment.quality
    if qualities and len(qualities) != len(sequence):
        raise ValueError("Read quality length differs from sequence length")
    for mi, mapping in enumerate(alignment.path.mapping):
        nid = mapping.position.node_id
        forward = sequences[nid]
        reverse = mapping.position.is_reverse
        reference = rc(forward) if reverse else forward
        cursor = mapping.position.offset
        initial = cursor
        first = len(columns)
        # Adjacent I or D edits encode one contiguous allele. Merge them before
        # applying the length limit, so 30I+21I cannot become supported 30/21 bp ALTs.
        edits = []
        for ei, original in enumerate(mapping.edit):
            f, t = original.from_length, original.to_length
            kind = "I" if f == 0 and t > 0 else "D" if t == 0 and f > 0 else None
            if edits and kind is not None and edits[-1].kind == kind:
                edits[-1].from_length += f
                edits[-1].to_length += t
                edits[-1].sequence += original.sequence
            else:
                edits.append(SimpleNamespace(from_length=f, to_length=t,
                    sequence=original.sequence, kind=kind, index=ei))
        for edit in edits:
            ei = edit.index
            f, t = edit.from_length, edit.to_length
            if f < 0 or t < 0 or cursor + f > len(reference) or read_cursor + t > len(sequence):
                raise ValueError(f"Invalid edit bounds: node {nid}, mapping {mi}, edit {ei}")
            ref = reference[cursor:cursor + f]
            alt = sequence[read_cursor:read_cursor + t]
            replacement = edit.sequence.upper()
            if replacement and (len(replacement) != t or replacement != alt):
                raise ValueError(f"Edit replacement disagrees with read at node {nid}")
            if not f and not t:
                continue
            op = "M" if f == t and not replacement else "X" if f == t else "I" if not f else "D" if not t else "C"
            pos = len(forward) - cursor - f if reverse else cursor
            cref, calt = (rc(ref), rc(alt)) if reverse else (ref, alt)
            qs = list(qualities[read_cursor:read_cursor+t]) if qualities else [-1] * t
            quality = sum(qs) / len(qs) if qs else -1
            # Deletions have no base quality; nearest existing read bases qualify them.
            if op == "D":
                flank = [qualities[q] for q in (read_cursor - 1, read_cursor)
                         if qualities and 0 <= q < len(qualities)]
                quality = min(flank) if flank else -1
            if op == "X":
                for k, (r, a) in enumerate(zip(cref, calt)):
                    if r != a:
                        q = qs[f - 1 - k if reverse else k]
                        observations.append(Observation(Candidate(nid, pos+k, r, a, "SNP"), mi, q))
            elif op in ("I", "D"):
                candidate = Candidate(nid, pos, cref, calt, "INS" if op == "I" else "DEL")
                if max(f, t) <= min(50, max_indel):
                    observations.append(Observation(candidate, mi, quality))
                else:
                    unsupported.append(dict(candidate.metadata(), mapping_index=mi, edit_index=ei,
                                            reason="indel_exceeds_limit"))
            elif op == "C":
                unsupported.append(dict(node_id=nid, start=pos, ref=cref, alt=calt,
                                        mapping_index=mi, edit_index=ei,
                                        reason="complex_replacement_not_supported"))
            for k in range(max(f, t)):
                is_boundary = k >= f
                p = cursor + min(k, f)
                p = len(forward) - p - (0 if is_boundary else 1) if reverse else p
                columns.append(Column(alt[k] if k < t else "-", ref[k] if k < f else "-",
                                      qs[k] if k < t else -1, op, nid, p, reverse, mi, is_boundary))
            cursor += f
            read_cursor += t
        lo, hi = (len(forward)-cursor, len(forward)-initial) if reverse else (initial, cursor)
        visits.append(Visit(nid, lo, hi, reverse, mi, first, len(columns), len(forward)))
    if read_cursor != len(sequence):
        raise ValueError("GAM edits do not consume the complete read sequence")
    return Read(alignment.name, hashlib.sha256(alignment.SerializeToString(deterministic=True)).hexdigest(),
                alignment.mapping_quality, columns, visits, observations), unsupported


def oriented_columns(read, visit, context=None):
    columns = read.columns if context is None else read.columns[
        max(0, visit.first-context):visit.last+context]
    return [c.flipped() for c in reversed(columns)] if visit.reverse else columns


def boundary_evidence(left, right, visit, boundary):
    """Immediate reference neighbors must meet this boundary, or a true node edge."""
    if not left or not right or left[-1].op not in ("M", "X") or right[0].op not in ("M", "X"):
        return False
    l, r = left[-1], right[0]
    before = l.pos == boundary-1 if l.visit == visit.index else boundary == 0
    after = r.pos == boundary if r.visit == visit.index else boundary == visit.node_length
    return before and after


def overlap(read, candidate, min_bq):
    """One count per record. Repeated visits: ALT > REF > other; first visit wins ties.

    Insertion coverage is closed-boundary overlap [start,end], including terminal
    boundaries and spanning deletions. REF needs adjacent M/X reference columns
    on both sides, with no insertion/replacement/deletion at the boundary.
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
        local = [c for c in cols if c.visit == visit.index]
        support = "other"
        if exact:
            support = "alt"
        elif candidate.kind == "INS":
            left, middle, right = split_columns(cols, visit, candidate)
            if not middle and boundary_evidence(left, right, visit, candidate.start):
                support = "ref"
        else:
            affected = [c for c in local if not c.boundary and candidate.start <= c.pos < candidate.end]
            inserted = any(c.boundary and candidate.start < c.pos < candidate.end for c in local)
            if (len(affected) == len(candidate.ref) and not inserted
                    and all(c.op in ("M", "X") and c.read == c.ref for c in affected)):
                support = "ref"
        choices.append((dict(alt=0, ref=1, other=2)[support], visit.index, support, visit))
    if not choices:
        return None
    _, _, support, visit = min(choices, key=lambda x: x[:2])
    return support, visit


def split_columns(cols, visit, candidate):
    """Separate candidate insertion boundary or reference interval from row context."""
    indices = [i for i, c in enumerate(cols) if c.visit == visit.index and
               ((c.boundary and c.pos == candidate.start) if candidate.kind == "INS" else
                (not c.boundary and candidate.start <= c.pos < candidate.end))]
    if indices:
        return cols[:indices[0]], cols[indices[0]:indices[-1]+1], cols[indices[-1]+1:]
    # Insertion with no inserted bases: place block immediately before the first
    # reference base at/after boundary, or after the last column in this visit.
    local = [i for i, c in enumerate(cols) if c.visit == visit.index]
    cut = next((i for i in local if cols[i].pos >= candidate.start),
               local[-1]+1 if local else 0)
    return cols[:cut], [], cols[cut:]


def make_tensor(candidate, eligible, rows=200, width=100, debug=False, node_walk_counts=None):
    """Select deterministically only after full-record coverage/support counting."""
    if rows < 1 or width < 1:
        raise ValueError("Tensor rows and width must be positive")
    span = max(1, len(candidate.alt) if candidate.kind == "INS" else len(candidate.ref))
    if candidate.kind == "INS":
        # Left-align insertion alleles; reserve up to the supported event limit.
        span = max([span] + [min(50, len(split_columns(oriented_columns(r, v, 1), v, candidate)[1]))
                             for r, _, v in eligible])
    if candidate.kind == "DEL":
        span += max([0] + [sum(c.boundary for c in split_columns(oriented_columns(r, v, 1), v, candidate)[1])
                           for r, _, v in eligible])
    if span > width:
        raise ValueError(f"Tensor width cannot preserve complete candidate region: "
                         f"candidate={candidate.metadata()['candidate_id']}, required_columns={span}, width={width}")
    start = (width - span) // 2
    with_walks = node_walk_counts is not None
    if with_walks and any(not isinstance(count, (int, np.integer)) or not 0 <= count <= np.iinfo(np.int32).max
                          for count in node_walk_counts.values()):
        raise ValueError("Node walk counts must be nonnegative int32 integers")
    tensor = np.zeros((7 if with_walks else 6, rows, width), dtype=np.int32 if with_walks else np.int16)
    details = []
    omitted = []
    windows = []
    for read, support, visit in eligible:
        cols = oriented_columns(read, visit, width)
        left, center, right = split_columns(cols, visit, candidate)
        overflow = 0
        if candidate.kind == "INS":
            overflow = max(0, len(center)-span)
            center = center[:span]
            # Gap slots are known aligned absence only with boundary evidence.
            evidence = bool(center) or boundary_evidence(left, right, visit, candidate.start)
            gap = Column("-", "-", -1, "G", candidate.node, candidate.start, False, visit.index, True)
            center += [gap if evidence else None] * (span - len(center))
        else:
            by_pos = {c.pos: c for c in center if not c.boundary}
            insertions = {}
            for col in center:
                if col.boundary:
                    insertions.setdefault(col.pos, []).append(col)
            center = []
            for p in range(candidate.start, candidate.end):
                center.extend(insertions.get(p, []))
                center.append(by_pos.get(p))
            # Row-specific context insertions stay on this row; unused central
            # slots are padding and never imply an insertion on another path.
            center += [None] * (span-len(center))
        row = [None] * max(0, start-len(left)) + left[-start:] if start else []
        row += center + right[:width-start-span]
        row += [None] * (width-len(row))
        # Group only the path actually visible in this candidate window. Distant
        # branches outside the window must not split otherwise identical groups.
        # Keep mapping transitions (including repeated visits), but ignore offsets
        # and insertion/gap lengths when identifying a node path.
        path = []
        previous_visit = None
        for col in row:
            if col is not None and col.visit != previous_visit:
                path.append((col.node, col.reverse))
                previous_visit = col.visit
        mismatch = sum(c is not None and c.op != "G" and c.read != c.ref for c in row)
        windows.append((read, support, visit, row, tuple(path), mismatch, overflow))
    # Rank all eligible records, then group stably BEFORE the depth limit.
    windows.sort(key=lambda w: (-w[5], -w[0].mapq, w[0].digest, w[2].index))
    path_groups = {}
    for window in windows:
        path_groups.setdefault(window[4], []).append(window)
    ordered = [window for group in path_groups.values() for window in group]
    n_selected = min(rows, len(ordered))
    if n_selected == 1:
        sample_indices = [len(ordered) // 2]
    elif n_selected:
        sample_indices = [i * (len(ordered)-1) // (n_selected-1) for i in range(n_selected)]
    else:
        sample_indices = []
    chosen = [ordered[i] for i in sample_indices]
    selected = [(w[0], w[1], w[2]) for w in chosen]
    groups = []
    for ri, (read, support, visit, row, path_key, mismatch, overflow) in enumerate(chosen):
        path = [dict(node_id=node, reverse=reverse) for node, reverse in path_key]
        if not groups or groups[-1]["path"] != path:
            groups.append(dict(start_row=ri, end_row=ri+1, path=path))
        else:
            groups[-1]["end_row"] = ri+1
        for ci, col in enumerate(row):
            if col is not None:
                tensor[:6, ri, ci] = [BASES.get(col.read, 5), col.quality,
                    int(col.op in ("I", "D", "C") or (col.op == "X" and col.read != col.ref)) | (2 if start <= ci < start+span else 0),
                    min(32767, read.mapq), OPS[col.op], BASES.get(col.ref, 5)]
                if with_walks:
                    tensor[6, ri, ci] = node_walk_counts[col.node]
            elif start <= ci < start+span:
                tensor[2, ri, ci] = 2
                if candidate.kind == "INS":
                    tensor[5, ri, ci] = BASES["-"]
        if overflow:
            omitted.append(dict(row_index=ri, omitted_columns=overflow,
                                reason="unsupported_insertion_allele_exceeds_shared_block"))
        if debug:
            details.append(dict(read_name=read.name, record_sha256=read.digest, support=support,
                anchor_mapping_index=visit.index, reversed_for_candidate=visit.reverse,
                window_mismatch_bp=mismatch, grouped_rank=sample_indices[ri],
                omitted_central_context_columns=overflow,
                path=[dict(node_id=v.node, start=v.start, end=v.end, reverse=v.reverse,
                           mapping_index=v.index) for v in read.visits],
                columns=[c.graph() if c is not None else None for c in row]))
    counts = Counter(s for _, s, _ in eligible)
    return tensor, dict(candidate.metadata(), tensor_format_version=V3_VERSION if with_walks else VERSION,
        row_order=ROW_ORDER, row_groups=groups, row_selection_version=ROW_SELECTION_VERSION,
        selected_grouped_ranks=sample_indices, window_mismatch_bp=[w[5] for w in chosen],
        candidate_columns=[start, start+span], coverage=len(eligible),
        alt_count=counts["alt"], ref_count=counts["ref"], other_count=counts["other"],
        af=counts["alt"]/len(eligible) if eligible else 0,
        selected_alignments=len(selected), omitted_context=omitted, selected_counts=dict(Counter(s for _, s, _ in selected)),
        **({"rows": details, "selection_audit": [dict(record_sha256=w[0].digest,
            mapping_quality=w[0].mapq, anchor_mapping_index=w[2].index,
            window_mismatch_bp=w[5], support=w[1],
            path=[list(p) for p in w[4]]) for w in windows]} if debug else {}))

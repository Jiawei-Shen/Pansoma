"""Candidate model: decode GAM edits, count allele support, encode one tensor per candidate.

All coordinates are zero-based on the *forward* strand of a node. Candidate
identity is the exact tuple (node, forward interval, REF, ALT, kind); there is
no left-normalization and no merging across nodes.

Tensor layout is (8, rows, width) int8. Every row is oriented to the forward strand of
the candidate node and the candidate starts at column width // 2. Cells without
evidence are 0 in every channel.
    0 read base        A=1 C=2 G=3 T=4 N=5 gap=6
    1 base quality     clip(q, -1, 127); -1 = missing quality or no read base
    2 candidate ALT    the candidate ALT base of that column (same encoding as 0; DEL columns = 6),
                       identical for every row, only inside the candidate region
    3 mapping quality  clip(MAPQ, -1, 127)
    4 operation        M=1 X=2 I=3 D=4 complex=5 aligned-no-insertion gap=6
    5 graph base       the graph reference base under this row/column, same encoding as 0
    6 path count       floor(14 * log2(distinct GBWT paths of the column's node + 1) + 0.5), max 127
    7 strand           1 = read sequenced on the candidate node's forward strand, 2 = reverse (row-level)

"Row matches REF" is channel 0 == channel 5; "row matches ALT" is channel 0 == channel 2.
"""
from collections import Counter, defaultdict
from bisect import bisect_left
from dataclasses import dataclass, field
import hashlib
import math
from types import SimpleNamespace

import numpy as np

FORMAT_VERSION = "indexed-gam-candidate-v5"
SCHEMA_VERSION = 5
STORAGE_VERSION = "int8-count-log2x14-v1"
ROW_SELECTION_VERSION = "window-edit-bp-group-uniform-v1"
WINDOW_ENCODING_VERSION = "anchor-centered-columns-v1"
ROW_ORDER = ("descending visible-window substitution+insertion+deletion bp; stable "
             "visible-node-path groups in first-occurrence order; uniform ordered sampling")
CHANNELS = ["read_base", "base_quality", "candidate_alt", "mapping_quality", "alignment_operation",
            "row_graph_reference_base", "node_distinct_gbwt_path_count", "strand"]
BASES = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5, "-": 6}
OPS = {"M": 1, "X": 2, "I": 3, "D": 4, "C": 5, "G": 6}
STRAND = {"forward": 1, "reverse": 2}
COUNT_LOG_SCALE = 14
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def int8_quality(value):
    """Keep the missing-quality sentinel -1 and saturate large qualities."""
    return max(-1, min(127, int(value)))


def encode_count(value):
    """int8 log scale for path counts: floor(14 * log2(count + 1) + 0.5), saturating at 127.

    Every count from 1 to 19 gets its own level (14, 22, 28, 33, ...); 90 -> 91, 331 -> 117,
    counts >= 524 -> 127. Invert with 2 ** (value / 14) - 1. Real nodes always have at least
    one path, so 0 only ever means "no evidence".
    """
    value = int(value)
    if value < 0:
        raise ValueError("Path counts must be nonnegative")
    return min(127, math.floor(COUNT_LOG_SCALE * math.log2(value + 1) + 0.5))


def candidate_alt_codes(candidate):
    """Channel-2 value of each candidate-region column: ALT bases, or gaps for a deletion."""
    if candidate.kind == "DEL":
        return [BASES["-"]] * len(candidate.ref)
    return [BASES.get(base, 5) for base in candidate.alt]


# --- decoded-read model --------------------------------------------------------

@dataclass(frozen=True, order=True)
class Candidate:
    node: int
    start: int
    ref: str
    alt: str
    kind: str  # SNP | INS | DEL

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
    """One alignment column: a read base (or gap) against a graph base (or gap)."""
    __slots__ = ("read", "ref", "quality", "op", "node", "pos", "reverse", "visit", "boundary")

    read: str
    ref: str
    quality: int
    op: str
    node: int
    pos: int        # forward-node offset; for boundary columns the insertion boundary
    reverse: bool
    visit: int      # mapping index within the read
    boundary: bool  # True for inserted bases (no graph base consumed)

    def flipped(self):
        return Column(rc(self.read), rc(self.ref), self.quality, self.op,
                      self.node, self.pos, not self.reverse, self.visit, self.boundary)

    def graph(self):
        return dict(node_id=self.node, **({"boundary": self.pos} if self.boundary
                    else {"offset": self.pos}), reverse=self.reverse, mapping_index=self.visit)


@dataclass
class Visit:
    """One mapping of the read to a node, in forward-node coordinates."""
    node: int
    start: int
    end: int
    reverse: bool
    index: int      # mapping index
    first: int      # first column index of this mapping
    last: int       # one past the last column index
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
    # Lazily built lookups (see visits_on / alt_quality); not part of the record's identity.
    _node_visits: dict = field(default=None, repr=False, compare=False)
    _alt_quality: dict = field(default=None, repr=False, compare=False)

    def visits_on(self, node):
        """This record's non-empty visits to `node`, in mapping order."""
        if self._node_visits is None:
            self._node_visits = defaultdict(list)
            for visit in self.visits:
                if visit.first != visit.last:
                    self._node_visits[visit.node].append(visit)
        return self._node_visits.get(node, ())

    def alt_quality(self, visit_index, candidate):
        """Best quality of an exact observation of `candidate` in one mapping, or None."""
        if self._alt_quality is None:
            best = {}
            for o in self.observations:
                key = (o.visit, o.candidate)
                if key not in best or o.quality > best[key]:
                    best[key] = o.quality
            self._alt_quality = best
        return self._alt_quality.get((visit_index, candidate))


def decode_alignment(alignment, sequences, max_indel=50, target_nodes=None):
    """Walk every edit of every mapping, producing columns, visits and candidate observations.

    Adjacent insertion (or deletion) edits within one mapping are merged before the
    length limit is applied. Substitution edits give SNP candidates per differing base;
    I/D edits up to `max_indel` give INS/DEL candidates; longer indels and complex
    replacements are returned as unsupported events (they remain as row context).
    Candidates containing N (or an insertion anchored on N) are context only.
    When `target_nodes` is given, observations are kept only for those nodes; columns,
    visits and unsupported events are unchanged.
    """
    columns, visits, observations, unsupported = [], [], [], []
    read_cursor = 0
    sequence = alignment.sequence.upper()
    qualities = alignment.quality
    if qualities and len(qualities) != len(sequence):
        raise ValueError("Read quality length differs from sequence length")
    for mi, mapping in enumerate(alignment.path.mapping):
        nid = mapping.position.node_id
        observe = target_nodes is None or nid in target_nodes
        forward = sequences[nid]
        reverse = mapping.position.is_reverse
        reference = rc(forward) if reverse else forward
        cursor = mapping.position.offset
        initial = cursor
        first = len(columns)
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
            op = ("M" if f == t and not replacement else "X" if f == t
                  else "I" if not f else "D" if not t else "C")
            pos = len(forward) - cursor - f if reverse else cursor
            cref, calt = (rc(ref), rc(alt)) if reverse else (ref, alt)
            qs = list(qualities[read_cursor:read_cursor + t]) if qualities else [-1] * t
            quality = sum(qs) / len(qs) if qs else -1
            if op == "D":  # a deletion has no bases: use the nearest read bases' quality
                flank = [qualities[q] for q in (read_cursor - 1, read_cursor)
                         if qualities and 0 <= q < len(qualities)]
                quality = min(flank) if flank else -1
            if op == "X" and observe:
                for k, (r, a) in enumerate(zip(cref, calt)):
                    if r != a and r.upper() != "N" and a.upper() != "N":
                        q = qs[f - 1 - k if reverse else k]
                        observations.append(Observation(Candidate(nid, pos + k, r, a, "SNP"), mi, q))
            elif op in ("I", "D"):
                candidate = Candidate(nid, pos, cref, calt, "INS" if op == "I" else "DEL")
                anchor = forward[max(0, pos - 1):max(0, pos - 1) + 1]  # forward-node anchor base
                ambiguous = ("N" in cref.upper() or "N" in calt.upper()
                             or (op == "I" and anchor.upper() == "N"))
                if max(f, t) > max_indel:
                    unsupported.append(dict(candidate.metadata(), mapping_index=mi, edit_index=ei,
                                            reason="indel_exceeds_limit"))
                elif observe and not ambiguous:
                    observations.append(Observation(candidate, mi, quality))
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
        lo, hi = (len(forward) - cursor, len(forward) - initial) if reverse else (initial, cursor)
        visits.append(Visit(nid, lo, hi, reverse, mi, first, len(columns), len(forward)))
    if read_cursor != len(sequence):
        raise ValueError("GAM edits do not consume the complete read sequence")
    digest = hashlib.sha256(alignment.SerializeToString(deterministic=True)).hexdigest()
    return Read(alignment.name, digest, alignment.mapping_quality, columns, visits, observations), unsupported


# --- support counting -----------------------------------------------------------

def oriented_columns(read, visit, context=None):
    """Columns of one visit (plus `context` columns each side), in candidate orientation."""
    columns = read.columns if context is None else read.columns[
        max(0, visit.first - context):visit.last + context]
    return [c.flipped() for c in reversed(columns)] if visit.reverse else columns


def boundary_evidence(left, right, visit, boundary):
    """The insertion boundary is flanked by matched/mismatched reference bases (or a node edge)."""
    if not left or not right or left[-1].op not in ("M", "X") or right[0].op not in ("M", "X"):
        return False
    l, r = left[-1], right[0]
    before = l.pos == boundary - 1 if l.visit == visit.index else boundary == 0
    after = r.pos == boundary if r.visit == visit.index else boundary == visit.node_length
    return before and after


def split_columns(cols, visit, candidate):
    """Split oriented columns into (before, candidate block, after)."""
    indices = [i for i, c in enumerate(cols) if c.visit == visit.index and
               ((c.boundary and c.pos == candidate.start) if candidate.kind == "INS" else
                (not c.boundary and candidate.start <= c.pos < candidate.end))]
    if indices:
        return cols[:indices[0]], cols[indices[0]:indices[-1] + 1], cols[indices[-1] + 1:]
    # Insertion absent from this read: cut immediately before the first reference
    # base at/after the boundary, or after the visit's last column.
    local = [i for i, c in enumerate(cols) if c.visit == visit.index]
    cut = next((i for i in local if cols[i].pos >= candidate.start), local[-1] + 1 if local else 0)
    return cols[:cut], [], cols[cut:]


class VisitView:
    """One visit of a record in candidate orientation, indexed for per-candidate lookups.

    `cols` is exactly oriented_columns(read, visit, context). The visit's own columns are
    indexed by forward position (graph bases) and by boundary (inserted bases), so REF
    checks and window cuts cost O(allele length) instead of a scan over the whole visit.
    A view is built once per (record, visit) and shared by every candidate on the node.
    If a visit ever breaks the invariants the indexes rely on (one column per graph base,
    non-decreasing positions), `indexed` is False and the plain scans are used instead.
    """
    __slots__ = ("read", "visit", "context", "cols", "base_at", "boundary_at", "local_pos", "local_idx", "indexed")

    def __init__(self, read, visit, context):
        self.read, self.visit, self.context = read, visit, context
        self.cols = oriented_columns(read, visit, context)
        self.base_at, self.boundary_at, self.local_pos, self.local_idx = {}, {}, [], []
        unique = True
        for i, c in enumerate(self.cols):
            if c.visit != visit.index:
                continue
            self.local_pos.append(c.pos)
            self.local_idx.append(i)
            if c.boundary:
                self.boundary_at[c.pos] = (self.boundary_at.get(c.pos, (i,))[0], i + 1)
            elif c.pos in self.base_at:
                unique = False
            else:
                self.base_at[c.pos] = i
        self.indexed = unique and all(a <= b for a, b in zip(self.local_pos, self.local_pos[1:]))

    def block(self, candidate):
        """(i, j): the candidate block cols[i:j], or (cut, cut) when this record has none.

        Same split as split_columns(): the block runs from the first to the last own column
        inside the candidate (boundary columns at the start for INS, graph bases otherwise);
        without one, the cut is before the first own column at/after the start.
        """
        if not self.indexed:
            left, middle, _ = split_columns(self.cols, self.visit, candidate)
            return len(left), len(left) + len(middle)
        if candidate.kind == "INS":
            if candidate.start in self.boundary_at:
                return self.boundary_at[candidate.start]
        else:
            hits = [self.base_at[p] for p in range(candidate.start, candidate.end) if p in self.base_at]
            if hits:
                return min(hits), max(hits) + 1
        k = bisect_left(self.local_pos, candidate.start)
        cut = self.local_idx[k] if k < len(self.local_idx) else (self.local_idx[-1] + 1 if self.local_idx else 0)
        return cut, cut

    def matches_reference(self, candidate):
        """REF support of this visit (no ALT): the rule documented in overlap()."""
        if candidate.kind == "INS":
            i, j = self.block(candidate)
            return i == j and boundary_evidence(self.cols[:i][-1:], self.cols[i:i + 1], self.visit, candidate.start)
        if not self.indexed:
            local = [c for c in self.cols if c.visit == self.visit.index]
            affected = [c for c in local if not c.boundary and candidate.start <= c.pos < candidate.end]
            inserted = any(c.boundary and candidate.start < c.pos < candidate.end for c in local)
            return (len(affected) == len(candidate.ref) and not inserted
                    and all(c.op in ("M", "X") and c.read == c.ref for c in affected))
        for p in range(candidate.start, candidate.end):
            i = self.base_at.get(p)
            if i is None:
                return False
            c = self.cols[i]
            if c.op not in ("M", "X") or c.read != c.ref:
                return False
        return not any(p in self.boundary_at for p in range(candidate.start + 1, candidate.end))


class NodeReads:
    """The records used for one node's candidates, with one VisitView per (record, visit).

    Views are built lazily with `context` flanking columns (the tensor width, so windows
    can reuse them) and live as long as this object: build one per node, then drop it.
    """

    def __init__(self, node, reads, context):
        self.node, self.reads, self.context = node, reads, context
        self.views = {}

    def view(self, read, visit):
        key = (id(read), visit.index)
        view = self.views.get(key)
        if view is None:
            view = self.views[key] = VisitView(read, visit, self.context)
        return view

    def classify(self, read, candidate, min_bq):
        """('alt'|'ref'|'other', anchor VisitView) for one record, or None if it does not cover.

        A record counts once even when it visits the node repeatedly: ALT beats REF beats
        other, and the earliest mapping breaks ties. ALT needs an exact observation with
        quality >= min_bq. REF needs matching graph bases over the whole interval; for an
        insertion, REF needs adjacent M/X columns on both sides of the boundary.
        """
        best = None
        for visit in read.visits_on(candidate.node):
            if candidate.kind == "INS":
                covered = visit.start <= candidate.start <= visit.end
            else:
                covered = visit.start < candidate.end and visit.end > candidate.start
            if not covered:
                continue
            quality = read.alt_quality(visit.index, candidate)
            if quality is not None and quality >= min_bq:
                return "alt", self.view(read, visit)  # visits are in mapping order: earliest ALT wins
            if best is None or best[0] == "other":  # an earlier REF already beats any later non-ALT
                view = self.view(read, visit)
                if view.matches_reference(candidate):
                    best = ("ref", view)
                elif best is None:
                    best = ("other", view)
        return best

    def eligible(self, candidate, min_bq):
        """[(record, support, anchor VisitView)] for every covering record, in record order."""
        result = []
        for read in self.reads:
            hit = self.classify(read, candidate, min_bq)
            if hit is not None:
                result.append((read, *hit))
        return result


def overlap(read, candidate, min_bq):
    """Classify one record as ('alt'|'ref'|'other', anchor visit), or None if it does not cover.

    See NodeReads.classify for the rule; this single-record form returns the Visit itself.
    """
    hit = NodeReads(candidate.node, [read], 1).classify(read, candidate, min_bq)
    return None if hit is None else (hit[0], hit[1].visit)


def alt_support_bounds(reads, candidates, min_bq):
    """Upper bound on ALT support: one vote per record per candidate (visits do not multiply)."""
    counts = Counter()
    for read in reads:
        counts.update({o.candidate for o in read.observations
                       if o.quality >= min_bq and o.candidate in candidates})
    return counts


def exact_coverage(candidates, reads):
    """{candidate: number of records overlap() would classify}, from visit intervals only.

    overlap() returns a hit iff one of the record's non-empty visits to the candidate
    node covers the candidate (closed interval for INS, half-open overlap otherwise),
    and each record counts once. That needs no column materialization, so this equals
    the eligible-record count of support counting over all records of `reads` on the
    candidate's node. One pass over all visits, then one vector test per candidate.
    """
    per_node = defaultdict(list)
    for candidate in candidates:
        per_node[candidate.node].append(candidate)
    visits = defaultdict(list)
    for ri, read in enumerate(reads):
        for v in read.visits:
            if v.first != v.last and v.node in per_node:
                visits[v.node].append((ri, v.start, v.end))
    result = {}
    for node, items in per_node.items():
        if not visits[node]:
            result.update((c, 0) for c in items)
            continue
        record, start, end = np.array(visits[node], dtype=np.int64).T
        for c in items:
            covered = (start <= c.start) & (end >= c.start) if c.kind == "INS" else (start < c.end) & (end > c.start)
            result[c] = int(np.unique(record[covered]).size)
    return result


# --- window encoding --------------------------------------------------------------

def anchor_window(read, visit, candidate, width):
    """Crop this record's columns so the candidate starts at column width // 2.

    Returns (row, cropped): `row` has exactly `width` entries (Column or None for
    missing evidence); `cropped` counts central columns lost past the right edge.
    INS rows without the insertion reserve gap ("G") slots only when both flanks
    are present; DEL rows keep one "D" column per deleted graph base.
    """
    return view_window(VisitView(read, visit, width), candidate, width)


def view_window(view, candidate, width):
    """anchor_window() on a prebuilt VisitView (needs view.context >= width)."""
    if view.context < width:
        raise ValueError("VisitView context is narrower than the window")
    anchor = width // 2
    cols, visit = view.cols, view.visit
    i, j = view.block(candidate)
    # Only the last `anchor` left columns and the first `width - anchor` right columns can
    # appear in the row (keep >= 1 on each side for the insertion flank check).
    left, center, right = cols[max(0, i - max(anchor, 1)):i], cols[i:j], cols[j:j + width - anchor]
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


BASE_CODES = np.full(256, BASES["N"], dtype=np.int64)  # any other character encodes as N, like BASES.get(x, 5)
for _base, _code in BASES.items():
    BASE_CODES[ord(_base)] = _code
OP_CODES = np.zeros(256, dtype=np.int64)
for _op, _code in OPS.items():
    OP_CODES[ord(_op)] = _code


def fill_rows(tensor, chosen, alt_stripe, path_counts, width):
    """Write the selected window rows into `tensor` (all eight channels, cells with a column only)."""
    flat = [col for w in chosen for col in w[3]]
    cells = np.flatnonzero(np.fromiter((c is not None for c in flat), dtype=bool, count=len(flat)))
    if not len(cells):
        return
    cols = [flat[k] for k in cells.tolist()]
    ri, ci = np.divmod(cells, width)

    def chars(values, table):
        raw = "".join(values).encode("latin-1", "replace")  # one character per column
        if len(raw) != len(cols):
            raise ValueError("Column bases and operations must be single characters")
        return table[np.frombuffer(raw, dtype=np.uint8)]

    ops = chars([c.op for c in cols], OP_CODES)
    if not ops.all():
        raise KeyError("Unknown alignment operation")
    nodes = np.array([c.node for c in cols], dtype=np.int64)
    used, where = np.unique(nodes, return_inverse=True)
    codes = []
    for node in used.tolist():
        count = path_counts[node]
        if not isinstance(count, (int, np.integer)) or not 0 <= count <= np.iinfo(np.int32).max:
            raise ValueError("Node path counts must be nonnegative int32 integers")
        codes.append(encode_count(count))
    stripe = np.zeros(width, dtype=np.int64)
    for column, code in alt_stripe.items():
        stripe[column] = code
    mapq = np.array([int8_quality(w[0].mapq) for w in chosen], dtype=np.int64)
    strand = np.array([STRAND["reverse"] if w[2].reverse else STRAND["forward"] for w in chosen], dtype=np.int64)
    quality = np.clip(np.array([c.quality for c in cols], dtype=np.int64), -1, 127)
    tensor[:, ri, ci] = np.stack([chars([c.read for c in cols], BASE_CODES), quality, stripe[ci], mapq[ri], ops,
                                  chars([c.ref for c in cols], BASE_CODES), np.array(codes, dtype=np.int64)[where],
                                  strand[ri]])


def make_tensor(candidate, eligible, path_counts, rows=200, width=101, debug=False):
    """Encode one candidate from `eligible` = [(read, support, anchor Visit or VisitView)].

    Every eligible record is windowed and ranked (descending visible edit bp, then
    MAPQ, record hash, mapping index), grouped stably by its visible node path, and
    only then uniformly sampled down to `rows`. Coverage/AF use all eligible records.
    """
    if rows < 1 or width < 1:
        raise ValueError("Tensor rows and width must be positive")
    start = width // 2
    # The candidate region is the same columns in every row: [start, start + allele length).
    alt_stripe = {start + k: code for k, code in enumerate(candidate_alt_codes(candidate)) if start + k < width}
    tensor = np.zeros((len(CHANNELS), rows, width), dtype=np.int8)
    windows = []
    for read, support, anchor in eligible:
        view = anchor if isinstance(anchor, VisitView) and anchor.context >= width else VisitView(
            read, anchor.visit if isinstance(anchor, VisitView) else anchor, width)
        visit = view.visit
        row, overflow = view_window(view, candidate, width)
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
    fill_rows(tensor, chosen, alt_stripe, path_counts, width)
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

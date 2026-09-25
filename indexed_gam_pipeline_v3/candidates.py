"""Candidate model: decode GAM edits, count allele support, encode one tensor per site.

All coordinates are zero-based on the *forward* strand of a node. Indels are
left-normalized while decoding (left_align_indels): every insertion/deletion moves to its
leftmost equivalent position on the node's forward strand (across short nodes of one
orientation), so forward- and reverse-strand reads report a repeat indel identically.
Candidate identity is then the exact tuple (node, forward interval, REF, ALT, kind).

Tensor layout is (8, rows, width) int8. Every row is oriented to the forward strand of
the site's node; the site (SiteLayout: insertion slots for the longest INS allele, then the
graph bases of the longest DEL allele / the SNV base) starts at column width // 2 in every
row. Cells without evidence are 0 in every channel.
    0 read base        A=1 C=2 G=3 T=4 N=5 gap=6
    1 base quality     clip(q, -1, 127); -1 = missing quality or no read base
    2 site allele      per row: the site allele this record carries, spelled over the site
                       columns (same encoding as 0; gaps = 6): an ALT's bases, or the REF
                       allele (slot gaps + graph bases); 0 for OTHER records and outside the site
    3 mapping quality  clip(MAPQ, -1, 127)
    4 operation        M=1 X=2 I=3 D=4 complex=5 aligned-no-insertion gap=6
    5 graph base       the graph reference base under this row/column, same encoding as 0
    6 path count       distinct GBWT paths of the column's node: exact up to 100, then
                       100 + ceil(log2(count - 99)) (331 -> 108), max 127
    7 strand           1 = read sequenced on the candidate node's forward strand, 2 = reverse (row-level)

Rows come in blocks A1 (best-supported ALT), A2, ..., REF, OTHER; inside a block similar
records are adjacent (similarity_order). "Row matches REF" is channel 0 == channel 5.
"""
from collections import Counter, defaultdict
from bisect import bisect_left
from dataclasses import dataclass, field
import hashlib
import math
from types import SimpleNamespace

import numpy as np

FORMAT_VERSION = "indexed-gam-candidate-v6"
SCHEMA_VERSION = 6
STORAGE_VERSION = "int8-count-linear100-log2-v1"
ROW_SELECTION_VERSION = "site-allele-blocks-uniform-similarity-v1"
WINDOW_ENCODING_VERSION = "site-layout-columns-v1"
ROW_ORDER = ("site-allele blocks A1..Ak, REF, OTHER (each by record hash), uniform ordered sampling, then "
             "average-linkage similarity order of window events inside each block (identical rows: strand, hash)")
CHANNELS = ["read_base", "base_quality", "site_allele", "mapping_quality", "alignment_operation",
            "row_graph_reference_base", "node_distinct_gbwt_path_count", "strand"]
BASES = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5, "-": 6}
OPS = {"M": 1, "X": 2, "I": 3, "D": 4, "C": 5, "G": 6}
STRAND = {"forward": 1, "reverse": 2}
COUNT_LINEAR_MAX = 100
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def int8_quality(value):
    """Keep the missing-quality sentinel -1 and saturate large qualities."""
    return max(-1, min(127, int(value)))


def encode_count(value):
    """int8 path counts: exact up to 100, then 100 + ceil(log2(count - 99)), saturating at 127.

    On HPRC v1.1 d9, 99.3 % of nodes have <= 100 paths (57 % have 85-90, near every
    haplotype), so those keep their exact count; the rest (repeat nodes revisited by
    one haplotype, max 331) are compressed: 101 -> 101, 102-103 -> 102, 104-107 -> 103,
    ..., 331 -> 108. Invert: value <= 100 is the count; value k > 100 means a count in
    (99 + 2 ** (k - 101), 99 + 2 ** (k - 100)]. Real nodes always have at least one
    path, so 0 only ever means "no evidence".
    """
    value = int(value)
    if value < 0:
        raise ValueError("Path counts must be nonnegative")
    if value <= COUNT_LINEAR_MAX:
        return value
    return min(127, COUNT_LINEAR_MAX + (value - COUNT_LINEAR_MAX).bit_length())  # = ceil(log2(value - 99))


# --- decoded-read model --------------------------------------------------------

@dataclass(frozen=True, order=True)
class Candidate:
    node: int
    start: int
    ref: str
    alt: str
    kind: str  # SNP | INS | DEL
    # A deletion continuing onto further nodes: ((node, start, end), ...) of those nodes, in
    # forward order; `ref` then spans all of them (the first node holds len(ref) - their bases).
    path: tuple = ()

    @property
    def end(self):
        return self.start + len(self.ref)

    def positions(self):
        """(node, forward offset) of every graph base of REF, in forward order."""
        first = len(self.ref) - sum(e - s for _, s, e in self.path)
        return ([(self.node, self.start + k) for k in range(first)]
                + [(node, p) for node, s, e in self.path for p in range(s, e)])

    def metadata(self):
        via = "@" + "+".join(str(node) for node, _, _ in self.path) if self.path else ""
        return dict(candidate_id=f"{self.node}:{self.start}:{self.kind}:{self.ref}>{self.alt}{via}",
                    node_id=self.node, start=self.start, end=self.end, orientation="+",
                    ref=self.ref, alt=self.alt, event_type=self.kind,
                    event_length=max(len(self.ref), len(self.alt)), path=[list(seg) for seg in self.path])


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
    # (from node, to node) of every indel that normalization moved to another node.
    moves: list = field(default_factory=list, repr=False, compare=False)
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


def decode_alignment(alignment, sequences, max_indel=50, target_nodes=None, left_align=True):
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
            if op == "X" and observe:
                for k, (r, a) in enumerate(zip(cref, calt)):
                    if r != a and r.upper() != "N" and a.upper() != "N":
                        q = qs[f - 1 - k if reverse else k]
                        observations.append(Observation(Candidate(nid, pos + k, r, a, "SNP"), mi, q))
            elif op in ("I", "D") and max(f, t) > max_indel:
                candidate = Candidate(nid, pos, cref, calt, "INS" if op == "I" else "DEL")
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
        lo, hi = (len(forward) - cursor, len(forward) - initial) if reverse else (initial, cursor)
        visits.append(Visit(nid, lo, hi, reverse, mi, first, len(columns), len(forward)))
    if read_cursor != len(sequence):
        raise ValueError("GAM edits do not consume the complete read sequence")
    moves = left_align_indels(columns, visits, max_indel) if left_align else []
    observations += indel_observations(columns, sequences, max_indel, target_nodes)
    for s, e, op in indel_runs(columns):  # joined over mappings, only the whole event is too long
        pieces = Counter(c.visit for c in columns[s:e])
        if e - s > max_indel and len(pieces) > 1 and max(pieces.values()) <= max_indel:
            first = columns[e - 1 if columns[s].reverse else s]
            unsupported.append(dict(node_id=first.node, start=first.pos, event_type="INS" if op == "I" else "DEL",
                                    event_length=e - s, mapping_index=first.visit, mappings=len(pieces),
                                    reason="indel_exceeds_limit_across_mappings"))
    digest = hashlib.sha256(alignment.SerializeToString(deterministic=True)).hexdigest()
    return Read(alignment.name, digest, alignment.mapping_quality, columns, visits, observations, moves), unsupported


# --- indel normalization --------------------------------------------------------

ACGT = frozenset("ACGT")


def indel_runs(columns):
    """(start, end, op) of every maximal run of I (inserted) or D columns.

    A run continues across mapping boundaries (same orientation): vg writes an indel
    spanning several nodes as one edit per mapping, but the read carries one event.
    """
    hi = len(columns)
    runs, i = [], 0
    while i < hi:
        c = columns[i]
        if c.op in ("I", "D"):
            j = i + 1
            while j < hi and columns[j].op == c.op and (columns[j].visit == c.visit or columns[j].reverse == c.reverse):
                j += 1
            runs.append((i, j, c.op))
            i = j
        else:
            i += 1
    return runs


def _shiftable(neighbour, reverse):
    """A matched graph base the indel may move across: M, read == graph base, A/C/G/T, same orientation."""
    return (neighbour.op == "M" and not neighbour.boundary and neighbour.reverse == reverse
            and neighbour.read == neighbour.ref and neighbour.read in ACGT)


def left_align_indels(columns, visits, max_indel):
    """Move every candidate-sized indel to its leftmost equivalent position on the forward strand.

    vg places a gap at one end of a repeat in *read* orientation, so forward- and reverse-
    strand reads put the same repeat indel at opposite ends in forward node coordinates.
    Here each insertion or deletion of at most `max_indel` bases, without N, is shifted one
    base at a time towards lower forward coordinates while the matched base it passes equals
    the indel's last forward base (the classic left-normalization). Read bases and qualities
    keep their order and every reference base keeps its column, so the read, its alignment
    span and every non-shifted column are unchanged; only which columns are I/D vs M moves.

    Shifts cross node boundaries within a run of mappings of one orientation. An indel that
    continues across mapping boundaries (vg writes one edit per node) is one event, and one that runs into another
    of the same kind merges with it (one event) and keeps moving. Finally every insertion is
    attached to the graph base that follows it on the forward strand, so an insertion
    between two nodes is always (next node, offset of its first base). `visits`
    get their column ranges updated in place; their graph intervals do not change.
    Returns [(from node, to node)] for every indel that ended on another node than it started.
    """
    runs = [r for r in indel_runs(columns) if r[1] - r[0] <= max_indel
            and not any("N" in (columns[k].read if r[2] == "I" else columns[k].ref).upper() for k in range(r[0], r[1]))]
    # Forward-strand runs move to lower read indices, reverse-strand runs to higher ones;
    # process each in the direction it moves so a run only ever meets already-final ones.
    forward = sorted((r for r in runs if not columns[r[0]].reverse), key=lambda r: r[0])
    reverse = sorted((r for r in runs if columns[r[0]].reverse), key=lambda r: -r[0])
    moves = []
    for s, e, op in forward + reverse:
        rev = columns[s].reverse
        origin = columns[s].node
        while True:
            k = e if rev else s - 1  # the column on the forward-left side
            if not 0 <= k < len(columns):
                break
            nb = columns[k]
            if nb.op == op and nb.reverse == rev:
                # An adjacent indel of the same kind is the same event: merge, then keep moving.
                j = k
                while 0 <= j + (1 if rev else -1) < len(columns):
                    c = columns[j + (1 if rev else -1)]
                    if c.op != op or c.reverse != rev:
                        break
                    j += 1 if rev else -1
                lo, hi = (s, j + 1) if rev else (j, e)
                merged = columns[lo:hi]
                if hi - lo > max_indel or any("N" in (c.read if op == "I" else c.ref).upper() for c in merged):
                    break
                s, e = lo, hi
                continue
            if not _shiftable(nb, rev):
                break
            if op == "I":
                edge = columns[s] if rev else columns[e - 1]  # the insertion's last forward base
                if nb.read != edge.read:
                    break
                bases = ([(c.read, c.quality) for c in columns[s + 1:e]] + [(nb.read, nb.quality)] if rev
                         else [(nb.read, nb.quality)] + [(c.read, c.quality) for c in columns[s:e - 1]])
                inserted = [Column(b, "-", q, "I", nb.node, nb.pos, nb.reverse, nb.visit, True) for b, q in bases]
                matched = Column(edge.read, nb.ref, edge.quality, "M", nb.node, nb.pos, nb.reverse, nb.visit, False)
                if rev:
                    columns[s:e + 1] = [matched] + inserted
                    s, e = s + 1, e + 1
                else:
                    columns[s - 1:e] = inserted + [matched]
                    s, e = s - 1, e - 1
            else:
                edge = columns[s] if rev else columns[e - 1]  # the deletion's last forward base
                if nb.ref != edge.ref:
                    break
                gap = Column("-", nb.ref, -1, "D", nb.node, nb.pos, nb.reverse, nb.visit, False)
                base = Column(nb.read, edge.ref, nb.quality, "M", edge.node, edge.pos, edge.reverse, edge.visit, False)
                columns[k] = gap
                if rev:
                    columns[s] = base
                    s, e = s + 1, e + 1
                else:
                    columns[e - 1] = base
                    s, e = s - 1, e - 1
        if op == "I":
            # Attach to the following forward base (the next read column on the forward strand);
            # at a read end without one, keep the whole insertion where its first forward base is.
            k = s - 1 if rev else e
            if 0 <= k < len(columns) and not columns[k].boundary and columns[k].reverse == rev:
                owner = columns[k]
            else:
                owner = columns[e - 1] if rev else columns[s]
            columns[s:e] = [Column(c.read, c.ref, c.quality, c.op, owner.node, owner.pos, c.reverse, owner.visit, True)
                            for c in columns[s:e]]
        if columns[s].node != origin:
            moves.append((origin, columns[s].node))
    spans = {}
    for i, c in enumerate(columns):
        spans.setdefault(c.visit, [i, i])[1] = i + 1
    edge = 0
    for visit in visits:  # mappings emptied by moving their only (inserted) columns keep an empty range
        visit.first, visit.last = spans.get(visit.index, (edge, edge))
        edge = visit.last
    return moves


def _read_base_quality(columns, i, step):
    """Quality of the nearest column holding a read base from index i in direction step (None: read edge)."""
    while 0 <= i < len(columns):
        if columns[i].read != "-":
            return None if columns[i].quality == -1 else columns[i].quality
        i += step
    return None


def _repeat_extent(columns, s, e, op, sequences):
    """Read indices of the matched bases a (leftmost) indel could still shift right across.

    Together with the indel's own columns they are the read bases of every equivalent
    placement of the event, the same for both strands' reports of one molecule.
    """
    rev = columns[s].reverse
    ordered = list(range(e - 1, s - 1, -1)) if rev else list(range(s, e))
    if op == "I":
        bases = [rc(columns[k].read) if rev else columns[k].read for k in ordered]
    else:
        bases = [sequences[columns[k].node][columns[k].pos] for k in ordered]
    passed, k = [], s - 1 if rev else e
    while 0 <= k < len(columns) and _shiftable(columns[k], rev):
        base = sequences[columns[k].node][columns[k].pos]
        if base.upper() != bases[0].upper():
            break
        bases = bases[1:] + [base]
        passed.append(k)
        k += -1 if rev else 1
    return passed


def indel_observations(columns, sequences, max_indel, target_nodes=None):
    """INS/DEL observations from the (normalized) columns: one per event, on its forward-first node.

    An event is a run of inserted or deleted columns, joined over mapping boundaries; it is
    observed when its forward-first column is on a target node. At most `max_indel` bases in
    total, no N in the alleles or in the base before an insertion (base 0 at a node start).
    A deletion over several nodes is one candidate with `path` = its further nodes.
    Quality, the same for either strand's report of the molecule: an insertion's is the mean
    over the read bases of all its equivalent placements (its bases plus the repeat bases it
    could shift across), a deletion's the lower of the nearest read bases outside that span
    (when the span reaches both read ends: the lowest read base inside it).
    """
    observations = []
    for s, e, op in indel_runs(columns):
        rev = columns[s].reverse
        first = columns[e - 1 if rev else s]  # forward-first column
        if e - s > max_indel or (target_nodes is not None and first.node not in target_nodes):
            continue
        passed = _repeat_extent(columns, s, e, op, sequences)
        forward = sequences[first.node]
        if op == "I":
            inserted = "".join(c.read for c in columns[s:e])
            alt = rc(inserted) if rev else inserted
            anchor = forward[max(0, first.pos - 1):max(0, first.pos - 1) + 1]
            if "N" in alt.upper() or anchor.upper() == "N":
                continue
            spanned = [columns[k].quality for k in range(s, e)] + [columns[k].quality for k in passed]
            quality = sum(spanned) / len(spanned)
            candidate = Candidate(first.node, first.pos, "", alt, "INS")
        else:
            ordered = columns[s:e][::-1] if rev else columns[s:e]
            ref = "".join(sequences[c.node][c.pos] for c in ordered)
            if "N" in ref.upper():
                continue
            segments = []
            for c in ordered:
                if segments and segments[-1][0] == c.node and segments[-1][2] == c.pos:
                    segments[-1][2] += 1
                else:
                    segments.append([c.node, c.pos, c.pos + 1])
            right = min(passed) - 1 if rev and passed else max(passed) + 1 if passed else (s - 1 if rev else e)
            flank = [q for q in (_read_base_quality(columns, e if rev else s - 1, 1 if rev else -1),
                                 _read_base_quality(columns, right, -1 if rev else 1)) if q is not None]
            if not flank:  # the repeat reaches both read ends: its own read bases are the evidence
                flank = [columns[k].quality for k in passed if columns[k].quality != -1]
            quality = min(flank) if flank else -1
            candidate = Candidate(first.node, segments[0][1], ref, "", "DEL", tuple(tuple(x) for x in segments[1:]))
        observations.append(Observation(candidate, first.visit, quality))
    return observations


# --- support counting -----------------------------------------------------------

def oriented_columns(read, visit, context):
    """Columns of one visit plus `context` columns each side, in candidate orientation."""
    columns = read.columns[max(0, visit.first - context):visit.last + context]
    return [c.flipped() for c in reversed(columns)] if visit.reverse else columns


def boundary_evidence(left, right, visit, boundary):
    """The insertion boundary is flanked by matched/mismatched reference bases (or a node edge)."""
    if not left or not right or left[-1].op not in ("M", "X") or right[0].op not in ("M", "X"):
        return False
    l, r = left[-1], right[0]
    before = l.pos == boundary - 1 if l.visit == visit.index else boundary == 0
    after = r.pos == boundary if r.visit == visit.index else boundary == visit.node_length
    return before and after


def aligned_across(left, right, visit, boundary):
    """The record is aligned across the boundary without inserting there: the insertion REF
    rule (boundary_evidence), or a matched/mismatched base followed by a deletion starting at it."""
    if boundary_evidence(left, right, visit, boundary):
        return True
    if not left or not right:
        return False
    l, r = left[-1], right[0]
    return (r.op == "D" and r.visit == visit.index and r.pos == boundary and l.op in ("M", "X")
            and (l.pos == boundary - 1 if l.visit == visit.index else boundary == 0))


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
        """REF support of this visit (no ALT): the rule documented in NodeReads.classify.

        SNV/DEL: the record has a plain match (M/X, read base = graph base, no inserted base) at
        every graph base of REF, over every node of a multi-node deletion, and aligned bases
        (M/X) right before and after it that are its graph neighbours (the previous/next base of
        the node, or another node at a node edge): no indel, no read end next to the site.
        INS: no inserted base at the boundary, and aligned bases on both sides of it.
        """
        cols = self.cols
        if candidate.kind == "INS":
            i, j = self.block(candidate)
            return i == j and boundary_evidence(cols[:i][-1:], cols[i:i + 1], self.visit, candidate.start)
        if self.indexed:
            i = self.base_at.get(candidate.start)
        else:
            i = next((k for k, c in enumerate(cols) if c.visit == self.visit.index and not c.boundary
                      and c.pos == candidate.start), None)
        expected = candidate.positions()
        j = None if i is None else i + len(expected)
        if i is None or i == 0 or j >= len(cols):
            return False
        own = len(expected) - sum(e - s for _, s, e in candidate.path)  # bases on this node: this visit's
        for k, (c, (node, pos)) in enumerate(zip(cols[i:j], expected)):
            if (c.boundary or (c.node, c.pos) != (node, pos) or c.op not in ("M", "X") or c.read != c.ref
                    or (k < own and c.visit != self.visit.index)):
                return False
        before, after, last = cols[i - 1], cols[j], cols[j - 1]
        if before.boundary or after.boundary or before.op not in ("M", "X") or after.op not in ("M", "X"):
            return False
        adjacent_before = (before.pos == candidate.start - 1 if before.visit == self.visit.index
                           else candidate.start == 0)
        adjacent_after = (after.pos == last.pos + 1 if after.visit == last.visit
                          else last.pos == self.read.visits[last.visit].node_length - 1)
        return adjacent_before and adjacent_after


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
        quality >= min_bq. REF needs a plain match at every graph base of REF plus aligned
        neighbouring bases on both sides (VisitView.matches_reference); for an insertion,
        REF needs adjacent M/X columns on both sides of the boundary and no inserted base.
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


def alt_support_bounds(reads, candidates, min_bq):
    """Upper bound on ALT support: one vote per record per candidate (visits do not multiply)."""
    counts = Counter()
    for read in reads:
        counts.update({o.candidate for o in read.observations
                       if o.quality >= min_bq and o.candidate in candidates})
    return counts


def exact_coverage(candidates, reads):
    """{candidate: number of records NodeReads.classify would classify}, from visit intervals only.

    NodeReads.classify returns a hit iff one of the record's non-empty visits to the candidate
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

@dataclass(frozen=True)
class SiteLayout:
    """The columns a site occupies from the anchor column on.

    `insertion_slots` inserted-base slots at boundary `start` (the longest INS allele),
    then the `span` graph bases [start, start + span) (the longest DEL allele, or the SNV
    base); `reference` is the graph sequence of that span. Every row of the site's tensor
    uses the same layout, so rows carrying alleles of different lengths stay aligned.
    """
    node: int
    start: int
    insertion_slots: int
    span: int
    reference: str

    @classmethod
    def of(cls, alleles):
        first = alleles[0]
        if any((a.node, a.start) != (first.node, first.start) for a in alleles):
            raise ValueError("Site alleles must share node and start")
        insertions = [len(a.alt) for a in alleles if a.kind == "INS"]
        spanned = max((a for a in alleles if a.kind != "INS"), key=lambda a: len(a.ref), default=None)
        return cls(first.node, first.start, max(insertions, default=0),
                   len(spanned.ref) if spanned else 0, spanned.ref if spanned else "")

    @property
    def columns(self):
        return self.insertion_slots + self.span

    def allele_codes(self, allele):
        """Channel-2 codes of one allele over the site columns (bases encoding; gap = 6).

        INS: its bases in the slots (gap-padded), then the graph bases. DEL: gaps over
        its deleted bases, graph bases for the rest of the span. SNV: the ALT base.
        None (REF): gaps in the slots, the graph bases over the span.
        """
        gap = BASES["-"]
        slots = [gap] * self.insertion_slots
        span = [BASES.get(b, 5) for b in self.reference]
        if allele is not None and allele.kind == "INS":
            slots = [BASES.get(b, 5) for b in allele.alt] + slots[len(allele.alt):]
        elif allele is not None and allele.kind == "DEL":
            span = [gap] * len(allele.ref) + span[len(allele.ref):]
        elif allele is not None:
            span = [BASES.get(b, 5) for b in allele.alt]
        return slots + span


def site_window(view, layout, width):
    """(row, cropped, cropped_inserted): this record's window over one site.

    `row` has exactly `width` entries (Column or None for missing evidence). Up to
    width // 2 left-flank columns end at the site; the site starts at column width // 2:
    the record's own inserted bases at the boundary fill the insertion slots (padded with
    "G" no-insertion gaps when both flanks are aligned M/X, else left empty; bases beyond
    the slots are cropped and counted in `cropped_inserted`; aligned_across decides "both
    flanks"), then its columns over the
    span (deleted bases keep one "D" column each; coverage starting inside the span keeps
    its true offset), then the right flank. With no insertion slots, an insertion at the
    boundary stays in the left flank. `cropped` counts site columns lost past the right edge.
    """
    if view.context < width:
        raise ValueError("VisitView context is narrower than the window")
    anchor = width // 2
    cols, visit = view.cols, view.visit
    i, j = view.block(Candidate(layout.node, layout.start, "", "", "INS"))
    ri, rj = (view.block(Candidate(layout.node, layout.start, layout.reference, "", "DEL")) if layout.span
              else (j, j))
    left_end = i if layout.insertion_slots else ri
    # Only the last `anchor` left columns and the first `width - anchor` right columns can
    # appear in the row (keep >= 1 on each side for the insertion flank check).
    left = cols[max(0, left_end - max(anchor, 1)):left_end]
    inserted, cropped_inserted, extent = [], 0, 0
    if layout.insertion_slots:
        own = cols[i:j]
        evidence = bool(own) or aligned_across(left, cols[j:j + 1], visit, layout.start)
        kept = own[:layout.insertion_slots]  # bases beyond the site's longest insertion are cropped
        cropped_inserted = len(own) - len(kept)
        gap = Column("-", "-", -1, "G", layout.node, layout.start, False, visit.index, True)
        inserted = kept + [gap if evidence else None] * (layout.insertion_slots - len(kept))
        extent = max(len(kept), min(layout.insertion_slots, width - anchor))  # columns that count as cropped
    spanned, missing = [], 0
    if layout.span:
        spanned = cols[ri:rj]
        # Coverage starting after the anchor keeps its true offset: pad the missing start.
        covered = [c.pos for c in spanned if not c.boundary]
        missing = min(width, max(0, min(covered) - layout.start)) if covered else 0
    right = cols[rj:rj + width - anchor]
    row = [None] * max(0, anchor - len(left)) + left[-anchor:] if anchor else []
    tail = inserted + [None] * missing + spanned + right
    row += tail[:width - anchor]
    row += [None] * (width - len(row))
    cropped = max(0, extent + missing + len(spanned) - (width - anchor))
    return row, cropped, cropped_inserted


BASE_CODES = np.full(256, BASES["N"], dtype=np.int64)  # any other character encodes as N, like BASES.get(x, 5)
for _base, _code in BASES.items():
    BASE_CODES[ord(_base)] = _code
OP_CODES = np.zeros(256, dtype=np.int64)
for _op, _code in OPS.items():
    OP_CODES[ord(_op)] = _code


def fill_rows(tensor, rows, stripes, path_counts, width):
    """Write window rows into `tensor` (all eight channels, cells with a column only).

    rows = [(read, anchor Visit, row columns)]; stripes = (len(rows), width) channel-2 codes.
    """
    flat = [col for _, _, row in rows for col in row]
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
    mapq = np.array([int8_quality(read.mapq) for read, _, _ in rows], dtype=np.int64)
    strand = np.array([STRAND["reverse"] if visit.reverse else STRAND["forward"] for _, visit, _ in rows],
                      dtype=np.int64)
    quality = np.clip(np.array([c.quality for c in cols], dtype=np.int64), -1, 127)
    tensor[:, ri, ci] = np.stack([chars([c.read for c in cols], BASE_CODES), quality, stripes[ri, ci], mapq[ri], ops,
                                  chars([c.ref for c in cols], BASE_CODES), np.array(codes, dtype=np.int64)[where],
                                  strand[ri]])


# --- row order ------------------------------------------------------------------

def window_events(row):
    """({event: first column}, (first, last) covered column) of one window row.

    Events are what distinguishes records in graph coordinates: every column that is not
    a plain match (mismatch base, inserted base, deletion, complex), keyed by (node,
    position, boundary, op, read base), and every node transition of the visible path.
    """
    events, previous, covered = {}, None, [i for i, c in enumerate(row) if c is not None]
    for ci, col in enumerate(row):
        if col is None:
            continue
        if previous is not None and (col.node, col.reverse) != previous[:2] and col.visit != previous[2]:
            events.setdefault(("path", previous[0], previous[1], col.node, col.reverse), ci)
        previous = (col.node, col.reverse, col.visit)
        if not (col.op in ("M", "G") and col.read == col.ref):
            events.setdefault((col.node, col.pos, col.boundary, col.op, col.read), ci)
    return events, (covered[0], covered[-1]) if covered else (0, -1)


def average_linkage(distance):
    """UPGMA merges of a square distance matrix as rows (a, b, height, size), numbered like scipy's
    linkage (leaves 0..n-1, the cluster made at step i is n + i); ties go to the first pair in
    row-major order of the current clusters. The fallback for scipy trees that are not valid."""
    n = len(distance)
    d = np.array(distance, dtype=float)
    np.fill_diagonal(d, np.inf)
    ids, sizes, tree = list(range(n)), [1] * n, []
    for step in range(n - 1):
        i, j = np.unravel_index(int(np.argmin(d)), d.shape)
        i, j = min(i, j), max(i, j)
        height = d[i, j]
        merged = (sizes[i] * d[i] + sizes[j] * d[j]) / (sizes[i] + sizes[j])
        tree.append((min(ids[i], ids[j]), max(ids[i], ids[j]), height, sizes[i] + sizes[j]))
        d[i, :] = d[:, i] = merged
        d[j, :] = d[:, j] = np.inf
        d[i, i] = np.inf
        ids[i], sizes[i] = n + step, sizes[i] + sizes[j]
    return np.array(tree, dtype=float).reshape(-1, 4)


def similarity_order(items):
    """Order window rows so similar records are adjacent; `items` = [(events, covered, tie_key)].

    Distance = number of events one record has and the other lacks although it covers that
    column (both ways). Uncovered columns are unknown, not different, and an event only one
    record has shifts that record's distance to all others equally, so it does not decide
    who is nearest. Average-linkage clustering; rows follow the tree's leaves with the larger
    subtree first (ties: smaller tie_key). Runs of rows with identical events are finally
    ordered by tie_key (strand, then record hash). Deterministic; no threshold.
    """
    n = len(items)
    if n < 3:
        order = sorted(range(n), key=lambda k: items[k][2])
    else:
        from scipy.cluster.hierarchy import is_valid_linkage, linkage
        vocabulary = {}
        for events, _, _ in items:
            for key, ci in events.items():
                vocabulary[key] = min(ci, vocabulary.get(key, ci))
        keys = list(vocabulary)
        index = {key: f for f, key in enumerate(keys)}
        has = np.zeros((n, len(keys)))
        for r, (events, _, _) in enumerate(items):
            has[r, [index[k] for k in events]] = 1
        where = np.array([vocabulary[k] for k in keys])
        lo = np.array([c[0] for _, c, _ in items])[:, None]
        hi = np.array([c[1] for _, c, _ in items])[:, None]
        covers = np.maximum(((lo <= where) & (where <= hi)).astype(float), has)
        distance = has @ covers.T + covers @ has.T - 2 * has @ has.T
        condensed = distance[np.triu_indices(n, 1)]
        tree = linkage(condensed, method="average")
        if not is_valid_linkage(tree):  # scipy 1.16 can merge a cluster with itself on tied distances
            tree = average_linkage(distance)
        rank = sorted(range(n), key=lambda k: items[k][2])
        first = {leaf: position for position, leaf in enumerate(rank)}
        members = {k: [k] for k in range(n)}
        for c, (a, b, _, _) in enumerate(tree, start=n):
            a, b = members.pop(int(a)), members.pop(int(b))
            if (len(b), -min(first[k] for k in b)) > (len(a), -min(first[k] for k in a)):
                a, b = b, a
            members[c] = a + b
        order = members[2 * n - 2]
    ordered, run = [], []
    for k in order + [None]:
        if run and (k is None or set(items[k][0]) != set(items[run[0]][0])):
            ordered += sorted(run, key=lambda r: items[r][2])
            run = []
        if k is not None:
            run.append(k)
    return ordered


# --- site tensors ---------------------------------------------------------------

REF_LABEL, OTHER_LABEL = "REF", "OTHER"


def site_labels(alleles, eligible_by_allele):
    """[(read, label, anchor)] in record order: one label per record over all site alleles.

    A record is allele A<k> for the first (best-supported) allele it carries, REF if it
    matches the reference for every allele it covers, otherwise OTHER. The anchor visit is
    that allele's anchor (REF/OTHER: the first covered allele's).
    """
    seen, order = {}, []
    for k, eligible in enumerate(eligible_by_allele):
        occurrence = Counter()
        for read, support, anchor in eligible:
            key = (id(read), occurrence[id(read)])  # the same record listed twice is two records
            occurrence[id(read)] += 1
            if key not in seen:
                seen[key] = (read, [])
                order.append(key)
            seen[key][1].append((k, support, anchor))
    result = []
    for key in order:
        read, hits = seen[key]
        alt = next(((k, anchor) for k, support, anchor in sorted(hits, key=lambda h: h[0]) if support == "alt"), None)
        if alt is not None:
            label, anchor = f"A{alt[0] + 1}", alt[1]
        else:
            label = REF_LABEL if all(support == "ref" for _, support, _ in hits) else OTHER_LABEL
            anchor = min(hits, key=lambda h: h[0])[2]
        result.append((read, label, anchor))
    return result


def make_site_tensor(alleles, eligible_by_allele, path_counts, rows=200, width=101, debug=False):
    """Encode one site: `alleles` in rank order (tensor representative first), each with its
    `eligible` = [(read, support, anchor Visit or VisitView)].

    Rows: every record eligible for any allele is labeled (site_labels), windowed over the
    site layout and put in blocks A1, A2, ..., REF, OTHER, each ordered by record hash;
    the concatenation is uniformly sampled down to `rows` (keeps each block's share), then
    each block is ordered by similarity (similarity_order). Channel 2 of a row spells the
    site allele that row carries (SiteLayout.allele_codes; OTHER rows 0).
    """
    if rows < 1 or width < 1:
        raise ValueError("Tensor rows and width must be positive")
    layout = SiteLayout.of(alleles)
    anchor_column = width // 2
    labels = [f"A{k + 1}" for k in range(len(alleles))]
    codes = {label: layout.allele_codes(allele) for label, allele in zip(labels, alleles)}
    codes[REF_LABEL] = layout.allele_codes(None)
    blocks = labels + [REF_LABEL, OTHER_LABEL]
    windows = []
    for read, label, anchor in site_labels(alleles, eligible_by_allele):
        view = anchor if isinstance(anchor, VisitView) and anchor.context >= width else VisitView(
            read, anchor.visit if isinstance(anchor, VisitView) else anchor, width)
        row, cropped, cropped_inserted = site_window(view, layout, width)
        mismatch = sum(c is not None and c.op != "G" and c.read != c.ref for c in row)
        windows.append(dict(read=read, label=label, visit=view.visit, row=row, cropped=cropped,
                            cropped_inserted=cropped_inserted, mismatch=mismatch,
                            tie=(view.visit.reverse, read.digest, view.visit.index)))
    ranked = sorted(windows, key=lambda w: (blocks.index(w["label"]), w["tie"][1], w["tie"][2]))
    n_selected = min(rows, len(ranked))
    if n_selected == 1:
        sample_indices = [len(ranked) // 2]
    elif n_selected:
        sample_indices = [i * (len(ranked) - 1) // (n_selected - 1) for i in range(n_selected)]
    else:
        sample_indices = []
    sampled = [(i, ranked[i]) for i in sample_indices]
    chosen, ranks, groups = [], [], []
    for label in blocks:
        block = [(i, w) for i, w in sampled if w["label"] == label]
        if not block:
            continue
        order = similarity_order([(*window_events(w["row"]), w["tie"]) for _, w in block])
        groups.append(dict(start_row=len(chosen), end_row=len(chosen) + len(block), allele=label))
        chosen += [block[k][1] for k in order]
        ranks += [block[k][0] for k in order]
    tensor = np.zeros((len(CHANNELS), rows, width), dtype=np.int8)
    stripes = np.zeros((len(chosen), width), dtype=np.int64)
    for ri, w in enumerate(chosen):
        for k, code in enumerate(codes.get(w["label"], ())):
            if anchor_column + k < width:
                stripes[ri, anchor_column + k] = code
    fill_rows(tensor, [(w["read"], w["visit"], w["row"]) for w in chosen], stripes, path_counts, width)
    omitted = [dict(row_index=ri, omitted_columns=w["cropped"], cropped_inserted_bases=w["cropped_inserted"],
                    reason="site_extends_beyond_window_or_insertion_longer_than_slots")
               for ri, w in enumerate(chosen) if w["cropped"] or w["cropped_inserted"]]
    representative = alleles[0]
    first = [(r, s, a) for r, s, a in eligible_by_allele[0]]
    counts = Counter(s for _, s, _ in first)
    site_counts = Counter(w["label"] for w in windows)
    meta = dict(representative.metadata(), tensor_storage_version=STORAGE_VERSION,
        tensor_format_version=FORMAT_VERSION, row_order=ROW_ORDER, row_groups=groups,
        row_selection_version=ROW_SELECTION_VERSION, selected_ranks=ranks,
        window_mismatch_bp=[w["mismatch"] for w in chosen], window_encoding_version=WINDOW_ENCODING_VERSION,
        anchor_column=anchor_column, candidate_columns=[anchor_column, min(width, anchor_column + layout.columns)],
        site_layout=dict(insertion_slots=layout.insertion_slots, span=layout.span, reference=layout.reference),
        allele_labels={label: allele.metadata()["candidate_id"] for label, allele in zip(labels, alleles)},
        coverage=len(first), alt_count=counts["alt"], ref_count=counts["ref"],
        other_count=counts["other"], af=counts["alt"] / len(first) if first else 0,
        site_coverage=len(windows), site_counts={b: site_counts[b] for b in blocks if site_counts[b]},
        selected_alignments=len(chosen), omitted_context=omitted,
        selected_counts={b: sum(w["label"] == b for w in chosen) for b in blocks if any(w["label"] == b for w in chosen)})
    if debug:
        meta["rows"] = [dict(read_name=w["read"].name, record_sha256=w["read"].digest, site_allele=w["label"],
            anchor_mapping_index=w["visit"].index, reversed_for_candidate=w["visit"].reverse,
            window_mismatch_bp=w["mismatch"], selected_rank=ranks[ri],
            omitted_central_context_columns=w["cropped"], cropped_inserted_bases=w["cropped_inserted"],
            path=[dict(node_id=v.node, start=v.start, end=v.end, reverse=v.reverse,
                       mapping_index=v.index) for v in w["read"].visits],
            columns=[c.graph() if c is not None else None for c in w["row"]]) for ri, w in enumerate(chosen)]
        meta["selection_audit"] = [dict(record_sha256=w["read"].digest, mapping_quality=w["read"].mapq,
            anchor_mapping_index=w["visit"].index, reversed_for_candidate=w["visit"].reverse,
            window_mismatch_bp=w["mismatch"], site_allele=w["label"]) for w in windows]
    return tensor, meta

"""Indel left-normalization: same molecule on either strand -> same candidates, read and alignment preserved."""
import copy
import random
import unittest

from fixtures import spec_alignment  # noqa: F401  (sets sys.path)
from indexed_gam_pipeline_v2.candidates import (Candidate, decode_alignment, indel_runs, left_align_indels, overlap,
                                                rc)


def mirror(specs, sequences):
    """The same molecule as the other strand would report it: reversed path, flipped mappings, mirrored edits."""
    mirrored = []
    for node, offset, reverse, edits in reversed(specs):
        span = sum(f for f, _, _ in edits)
        mirrored.append((node, len(sequences[node]) - offset - span, not reverse,
                         [(f, t, rc(seq) if seq else seq) for f, t, seq in reversed(edits)]))
    return mirrored


def decode(specs, sequences, **kwargs):
    return decode_alignment(spec_alignment(specs, sequences), sequences, **kwargs)[0]


def candidates(read):
    return sorted({o.candidate for o in read.observations})


def columns_of(read):
    return [(c.node, c.pos, c.boundary, c.op, c.read, c.ref, c.reverse) for c in read.columns]


def mirrored_columns(read):
    """The columns the other strand's report must have: reversed read order, each column flipped."""
    return [(c.node, c.pos, c.boundary, c.op, rc(c.read), rc(c.ref), not c.reverse) for c in reversed(read.columns)]


class LeftAlignTest(unittest.TestCase):
    def test_homopolymer_and_repeat_insertions_and_deletions(self):
        sequences = {1: "GAAAAC", 2: "CATATG", 3: "GTTTTC"}
        cases = [
            ((1, 0, [(5, 5, ""), (0, 1, "A"), (1, 1, "")]), Candidate(1, 1, "", "A", "INS")),       # after last A
            ((2, 0, [(5, 5, ""), (0, 2, "AT"), (1, 1, "")]), Candidate(2, 1, "", "AT", "INS")),     # rotates AT/TA
            ((3, 0, [(4, 4, ""), (1, 0, ""), (1, 1, "")]), Candidate(3, 1, "T", "", "DEL")),        # last T deleted
            ((3, 0, [(1, 1, ""), (0, 1, "T"), (5, 5, "")]), Candidate(3, 1, "", "T", "INS")),       # already leftmost
        ]
        for (node, offset, edits), expected in cases:
            for reverse in (False, True):
                specs = [(node, offset, False, edits)]
                if reverse:
                    specs = mirror(specs, sequences)
                with self.subTest(expected=expected, reverse=reverse):
                    read = decode(specs, sequences)
                    self.assertEqual(candidates(read), [expected])
                    self.assertEqual(overlap(read, expected, 10)[0], "alt")

    def test_insertion_crosses_short_repeat_nodes_and_attaches_to_the_next_node(self):
        # A TA repeat chopped into short nodes: C | AT | AT | G. +AT after the second AT node
        # (forward reads) and before the first (reverse reads) is one event: (node 2, offset 0).
        sequences = {1: "C", 2: "AT", 3: "AT", 4: "G"}
        specs = [(1, 0, False, [(1, 1, "")]), (2, 0, False, [(2, 2, "")]), (3, 0, False, [(2, 2, ""), (0, 2, "AT")]),
                 (4, 0, False, [(1, 1, "")])]
        for strand_specs in (specs, mirror(specs, sequences)):
            read = decode(strand_specs, sequences)
            self.assertEqual(candidates(read), [Candidate(2, 0, "", "AT", "INS")])
        # The previous node's end is the same boundary: it is attached to the next node's start.
        specs = [(1, 0, False, [(1, 1, ""), (0, 1, "G")]), (2, 0, False, [(2, 2, "")])]
        self.assertEqual(candidates(decode(specs, {1: "C", 2: "AT"})), [Candidate(2, 0, "", "G", "INS")])

    def test_blockers_stop_the_shift(self):
        sequences = {1: "GAAAAC", 2: "AA", 3: "AAC"}
        # a mismatch in the run
        read = decode([(1, 0, False, [(2, 2, ""), (1, 1, "T"), (2, 2, ""), (0, 1, "A"), (1, 1, "")])], sequences)
        self.assertIn(Candidate(1, 3, "", "A", "INS"), candidates(read))
        # another indel directly left
        read = decode([(1, 0, False, [(3, 3, ""), (1, 0, ""), (0, 1, "A"), (2, 2, "")])], sequences)
        self.assertEqual(candidates(read), [Candidate(1, 1, "A", "", "DEL"), Candidate(1, 2, "", "A", "INS")])
        # a 2-bp deletion stops at a node boundary; a 1-bp deletion crosses it
        read = decode([(2, 0, False, [(2, 2, "")]), (3, 0, False, [(2, 0, ""), (1, 1, "")])], sequences)
        self.assertEqual(candidates(read), [Candidate(3, 0, "AA", "", "DEL")])
        read = decode([(2, 0, False, [(2, 2, "")]), (3, 0, False, [(1, 1, ""), (1, 0, ""), (1, 1, "")])], sequences)
        self.assertEqual(candidates(read), [Candidate(2, 0, "A", "", "DEL")])
        # an orientation change stops it
        read = decode([(2, 0, True, [(2, 2, "")]), (3, 0, False, [(0, 1, "A"), (3, 3, "")])], sequences)
        self.assertEqual(candidates(read), [Candidate(3, 0, "", "A", "INS")])
        # N bases are never moved, nor moved across
        read = decode([(1, 0, False, [(5, 5, ""), (0, 1, "N"), (1, 1, "")])], sequences)
        self.assertEqual(candidates(read), [])
        self.assertEqual([c.pos for c in read.columns if c.boundary], [5])

    def test_disabled_keeps_vg_placement(self):
        sequences = {1: "GAAAAC"}
        read = decode([(1, 0, False, [(5, 5, ""), (0, 1, "A"), (1, 1, "")])], sequences, left_align=False)
        self.assertEqual(candidates(read), [Candidate(1, 5, "", "A", "INS")])

    def test_random_strand_symmetry_read_preservation_idempotence_and_leftmost(self):
        rng = random.Random(20260923)
        for case in range(1500):
            sequences = {n: "".join(rng.choice("ACGT" if rng.random() < .5 else "AT") for _ in range(rng.randint(2, 12)))
                         for n in range(1, 5)}
            reverse, mixed = rng.random() < 0.3, rng.random() < 0.3
            specs, nodes = [], sorted(sequences)
            for i in range(rng.randint(1, 4)):  # a walk (offset only on the first mapping), sometimes mixed-orientation
                node = nodes[i]
                offset = rng.randint(0, len(sequences[node]) - 1) if i == 0 else 0
                orientation = rng.random() < 0.5 if mixed else reverse
                specs.append((node, offset, orientation, random_edits(rng, len(sequences[node]), offset)))
            a = spec_alignment(specs, sequences)
            with self.subTest(case=case, specs=specs):
                read = decode_alignment(a, sequences)[0]
                raw = decode_alignment(a, sequences, left_align=False)[0]
                # the read, its qualities in read order and its graph span are unchanged
                self.assertEqual("".join(c.read for c in read.columns if c.read != "-"), a.sequence.upper())
                self.assertEqual([c.quality for c in read.columns if c.read != "-"],
                                 [c.quality for c in raw.columns if c.read != "-"])
                self.assertEqual([(c.node, c.pos) for c in read.columns if not c.boundary],
                                 [(c.node, c.pos) for c in raw.columns if not c.boundary])
                self.assertEqual([(v.node, v.start, v.end) for v in read.visits],
                                 [(v.node, v.start, v.end) for v in raw.visits])
                self.assertTrue(all(c.read == c.ref for c in read.columns if c.op == "M"))
                for v in read.visits:
                    self.assertTrue(all(c.visit == v.index for c in read.columns[v.first:v.last]))
                # idempotent
                again = copy.deepcopy(read.columns)
                left_align_indels(again, copy.deepcopy(read.visits), 50)
                self.assertEqual([vars_(c) for c in again], [vars_(c) for c in read.columns])
                # leftmost: no normalized indel can move further
                for s, e, op in indel_runs(read.columns):
                    rev = read.columns[s].reverse
                    k = e if rev else s - 1
                    if not 0 <= k < len(read.columns):
                        continue
                    nb, edge = read.columns[k], read.columns[s] if rev else read.columns[e - 1]
                    movable = (nb.op == "M" and not nb.boundary and nb.reverse == rev and nb.read in "ACGT"
                               and (nb.read == edge.read if op == "I" else nb.ref == edge.ref and
                                    (e - s == 1 or nb.visit == edge.visit)))
                    self.assertFalse(movable and "N" not in "".join(c.read + c.ref for c in read.columns[s:e]),
                                     (s, e, op))
                # the other strand's report of the same molecule gives the same candidates and columns
                other = decode(mirror(specs, sequences), sequences)
                self.assertEqual(candidates(other), candidates(read))
                self.assertEqual(mirrored_columns(other), columns_of(read))


def random_edits(rng, length, offset):
    """M/X/I/D edits (no N, no complex replacements) consuming [offset, end) of a node."""
    edits, cursor = [], offset
    for _ in range(rng.randint(1, 7)):
        room = length - cursor
        kind = rng.choice("MMXIIDD")
        if kind == "I":
            inserted = "".join(rng.choice("ACGT" if rng.random() < .5 else "AT") for _ in range(rng.randint(1, 4)))
            edits.append((0, len(inserted), inserted))
        elif room == 0:
            continue
        elif kind == "M":
            f = rng.randint(1, room)
            edits.append((f, f, ""))
        elif kind == "X":
            f = rng.randint(1, min(room, 2))
            edits.append((f, f, "".join(rng.choice("ACGT") for _ in range(f))))
        else:
            edits.append((rng.randint(1, min(room, 3)), 0, ""))
        cursor += edits[-1][0]
    if not edits:
        edits.append((0, 1, "A"))
    return edits


def vars_(column):
    return tuple(getattr(column, k) for k in column.__slots__)


if __name__ == "__main__":
    unittest.main()

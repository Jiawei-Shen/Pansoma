"""NodeReads / VisitView on random multi-mapping records.

Random records (both orientations, repeated node visits, M/X/I/D/complex edits, N bases, partial
node coverage, missing or varied qualities) are decoded, then every observed and probed candidate is
classified record by record (fixtures.overlap) and by NodeReads at several widths:
- NodeReads.eligible equals the per-record classification at every width (views are context independent);
- every one-node SNV/DEL REF call passes aligned_neighbours, an independent check of the REF rule's
  aligned-neighbour condition;
- site tensors do not depend on record order, and selected_alignments == min(rows, site_coverage);
- the forced scan fallback (VisitView.indexed = False) gives identical eligibility and tensors;
- a prebuilt narrow view is rebuilt for wide windows.

These are the parts of v2's test_reference_equivalence.py that do not need the frozen v5 oracle
(random_edits, random_reads, probes and aligned_neighbours are copied from it); byte identity with v2
is held by the goldens.
"""
import random
import unittest
from unittest import mock

import numpy as np

from .. import candidates
from ..candidates import Candidate, NodeReads, decode_alignment, oriented_columns
from .fixtures import make_tensor, overlap, spec_alignment

BASES = "ACGT"
LONG_PROBE_STRIDE = 5  # long nodes: every observed candidate plus every 5th synthetic probe (keeps the file near 15 s)


def random_edits(rng, reference_length, offset):
    """Edits consuming reference [offset, end) with at least one column."""
    edits, cursor = [], offset
    for _ in range(rng.randint(1, 6)):
        room = reference_length - cursor
        kind = rng.choice("MMXXIDC")
        if kind == "I":
            inserted = "".join(rng.choice(BASES + "N") for _ in range(rng.randint(1, 4)))
            edits.append((0, len(inserted), inserted))
        elif room == 0:
            continue
        elif kind == "M":
            f = rng.randint(1, room)
            edits.append((f, f, ""))
        elif kind == "X":
            f = rng.randint(1, min(room, 3))
            edits.append((f, f, "".join(rng.choice(BASES + "N") for _ in range(f))))
        elif kind == "D":
            edits.append((rng.randint(1, min(room, 4)), 0, ""))
        else:
            f = rng.randint(1, min(room, 3))
            t = rng.choice([x for x in range(1, 5) if x != f])
            edits.append((f, t, "".join(rng.choice(BASES) for _ in range(t))))
        cursor += edits[-1][0]
    if not edits:
        edits.append((0, 1, rng.choice(BASES)))
    return edits


def random_reads(rng, sequences):
    reads = []
    for ri in range(rng.randint(1, 7)):
        specs = []
        for _ in range(rng.randint(1, 4)):
            node = rng.choice(sorted(sequences))
            offset = rng.randint(0, len(sequences[node]))
            specs.append((node, offset, rng.random() < 0.5, random_edits(rng, len(sequences[node]), offset)))
        a = spec_alignment(specs, sequences, name=f"r{ri}", mapq=rng.choice([0, 20, 60, 255]))
        choice = rng.random()
        a.quality = (b"" if choice < 0.15 else bytes(rng.randint(0, 255) for _ in a.sequence) if choice < 0.6
                     else bytes([30] * len(a.sequence)))
        reads.append(decode_alignment(a, sequences, max_indel=rng.choice([2, 50]))[0])
    return reads


def probes(reads, sequences):
    found = {o.candidate for r in reads for o in r.observations}
    for node, s in sequences.items():
        for pos in range(len(s) + 1):
            found.add(Candidate(node, pos, "", "G", "INS"))
            for length in (1, 2, 3):
                if pos + length <= len(s):
                    found.add(Candidate(node, pos, s[pos:pos + length], "", "DEL"))
            if pos < len(s):
                found.add(Candidate(node, pos, s[pos], "C" if s[pos] != "C" else "A", "SNP"))
    return sorted(found)


def aligned_neighbours(read, visit, candidate):
    """Independent check of v6's extra REF condition for a one-node SNV/DEL: the columns right before
    and after the interval (candidate orientation) are M/X graph neighbours, no inserted base at either end."""
    cols = oriented_columns(read, visit, 3 + len(candidate.ref))
    own = [k for k, c in enumerate(cols) if c.visit == visit.index and not c.boundary
           and candidate.start <= c.pos < candidate.end]
    if not own or own[0] == 0 or own[-1] + 1 >= len(cols):
        return False
    before, after = cols[own[0] - 1], cols[own[-1] + 1]
    if before.boundary or after.boundary or before.op not in "MX" or after.op not in "MX":
        return False
    ok_before = before.pos == candidate.start - 1 if before.visit == visit.index else candidate.start == 0
    ok_after = after.pos == candidate.end if after.visit == visit.index else candidate.end == visit.node_length
    return ok_before and ok_after


def short_sequences(rng):
    return {n: "".join(rng.choice(BASES + ("N" if rng.random() < 0.2 else "")) for _ in range(rng.randint(1, 9)))
            for n in range(1, rng.randint(2, 4))}


class ViewTest(unittest.TestCase):
    def check_scenario(self, rng, sequences, widths, probe_stride=1):
        """Assert the view properties on one random scenario; returns what was computed (for comparisons).

        Every candidate a record carries is checked; of the other probes, every probe_stride-th.
        """
        reads = random_reads(rng, sequences)
        counts = {n: rng.choice([1, 2, 90, 600]) for n in sequences}
        observed = {o.candidate for r in reads for o in r.observations}
        computed, ref_calls = [], 0
        for k, candidate in enumerate(probes(reads, sequences)):
            if k % probe_stride and candidate not in observed:
                continue
            for min_bq in (0, 10, 40):
                classified = [(r, overlap(r, candidate, min_bq)) for r in reads]
                expected = [(r.name, *hit) for r, hit in classified if hit]
                for r, hit in classified:
                    if hit and hit[0] == "ref" and candidate.kind != "INS" and not candidate.path:
                        self.assertTrue(aligned_neighbours(r, hit[1], candidate), candidate)
                        ref_calls += 1
                for width in widths:
                    actual = NodeReads(candidate.node, reads, width).eligible(candidate, min_bq)
                    found = [(r.name, s, v.visit) for r, s, v in actual]
                    self.assertEqual(found, [(name, s, v) for name, s, v in expected], candidate)
                    computed.append((candidate, min_bq, width, found))
                    if min_bq == 10 and width == widths[-1]:
                        for rows in (1, 5, 200):  # tensors are well-formed and independent of record order
                            x, m = make_tensor(candidate, actual, counts, rows=rows, width=width, debug=True)
                            y, _ = make_tensor(candidate, actual[::-1], counts, rows=rows, width=width, debug=True)
                            np.testing.assert_array_equal(x, y, err_msg=str(candidate))
                            self.assertEqual(m["selected_alignments"], min(rows, m["site_coverage"]))
                            computed.append((candidate, rows, x.tobytes(), m))
        return computed, ref_calls

    def short_scenarios(self, seed, cases):
        rng, computed, ref_calls = random.Random(seed), [], 0
        for case in range(cases):
            sequences = short_sequences(rng)
            with self.subTest(case=case):
                found, refs = self.check_scenario(rng, sequences, widths=(1, 2, 3, 11))
                computed.append(found)
                ref_calls += refs
        return computed, ref_calls

    def test_random_short_nodes(self):
        _, ref_calls = self.short_scenarios(20260923, 60)
        self.assertGreater(ref_calls, 100)  # the aligned-neighbour oracle was exercised

    def test_scan_fallback_gives_identical_eligibility_and_tensors(self):
        """Views that fail their invariants use plain scans; forcing that path everywhere changes nothing."""
        indexed, _ = self.short_scenarios(11, 20)
        build, views = candidates.VisitView.__init__, []

        def unindexed(view, *args):
            build(view, *args)
            view.indexed = False
            views.append(view)
        with mock.patch.object(candidates.VisitView, "__init__", unindexed):
            fallback, _ = self.short_scenarios(11, 20)
        self.assertTrue(views)
        self.assertFalse(any(view.indexed for view in views))
        self.assertEqual(len(fallback), len(indexed))
        for case, (a, b) in enumerate(zip(indexed, fallback)):
            self.assertEqual(a, b, f"case {case}")

    def test_random_long_nodes_and_default_width(self):
        rng = random.Random(7)
        for case in range(4):
            sequences = {n: "".join(rng.choice(BASES) for _ in range(rng.randint(40, 160))) for n in (1, 2)}
            with self.subTest(case=case):
                self.check_scenario(rng, sequences, widths=(11, 101), probe_stride=LONG_PROBE_STRIDE)

    def test_prebuilt_narrow_view_is_rebuilt_for_wide_windows(self):
        rng = random.Random(3)
        sequences = {1: "ACGTACGTAC", 2: "TTGCA"}
        reads = random_reads(rng, sequences)
        for candidate in probes(reads, sequences):
            narrow = NodeReads(candidate.node, reads, 1).eligible(candidate, 10)
            wide = [(r, s, v.visit) for r, s, v in narrow]
            x, m = make_tensor(candidate, narrow, {1: 5, 2: 5}, rows=4, width=11, debug=True)
            y, n = make_tensor(candidate, wide, {1: 5, 2: 5}, rows=4, width=11, debug=True)
            np.testing.assert_array_equal(x, y)
            self.assertEqual(m, n)


if __name__ == "__main__":
    unittest.main()

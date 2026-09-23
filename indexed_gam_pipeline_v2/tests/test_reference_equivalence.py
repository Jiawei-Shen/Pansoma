"""The indexed VisitView/NodeReads path equals the frozen reference implementation exactly.

Random multi-mapping records (both orientations, repeated node visits, M/X/I/D/complex
edits, N bases, partial node coverage, missing or varied qualities) are decoded, then every
observed and probed candidate is classified, windowed and encoded by both implementations
at several widths and row counts. Support, anchor visit, rows, tensors and metadata must match.
"""
import random
import unittest
from unittest import mock

import numpy as np

from fixtures import spec_alignment  # noqa: F401  (sets sys.path)
import reference_candidates as reference
from indexed_gam_pipeline_v2 import candidates
from indexed_gam_pipeline_v2.candidates import (Candidate, NodeReads, anchor_window, decode_alignment, make_tensor,
                                                overlap)

BASES = "ACGT"


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


class ReferenceEquivalenceTest(unittest.TestCase):
    def check_scenario(self, rng, sequences, widths=(1, 2, 3, 11, 101)):
        reads = random_reads(rng, sequences)
        counts = {n: rng.choice([1, 2, 90, 600]) for n in sequences}
        for candidate in probes(reads, sequences):
            for min_bq in (0, 10, 40):
                expected = [(r, *reference.overlap(r, candidate, min_bq)) for r in reads
                            if reference.overlap(r, candidate, min_bq)]
                for r in reads:
                    new, old = overlap(r, candidate, min_bq), reference.overlap(r, candidate, min_bq)
                    self.assertEqual(new is None, old is None, candidate)
                    if old is not None:
                        self.assertEqual(new[0], old[0], candidate)
                        self.assertIs(new[1], old[1], candidate)
                for width in widths:
                    node_reads = NodeReads(candidate.node, reads, width)
                    actual = node_reads.eligible(candidate, min_bq)
                    self.assertEqual([(r.name, s, v.visit) for r, s, v in actual],
                                     [(r.name, s, v) for r, s, v in expected], candidate)
                    for r, _, visit in expected:
                        self.assertEqual(anchor_window(r, visit, candidate, width),
                                         reference.anchor_window(r, visit, candidate, width), candidate)
                    if min_bq != 10:
                        continue
                    for rows in (1, 2, 5, 200):
                        x, m = make_tensor(candidate, actual, counts, rows=rows, width=width, debug=True)
                        y, n = reference.make_tensor(candidate, expected, counts, rows=rows, width=width, debug=True)
                        np.testing.assert_array_equal(x, y, err_msg=str(candidate))
                        self.assertEqual(m, n, candidate)

    def short_scenarios(self, seed, cases):
        rng = random.Random(seed)
        for case in range(cases):
            sequences = {n: "".join(rng.choice(BASES + ("N" if rng.random() < 0.2 else "")) for _ in range(rng.randint(1, 9)))
                         for n in range(1, rng.randint(2, 4))}
            with self.subTest(case=case):
                self.check_scenario(rng, sequences, widths=(1, 2, 3, 11))

    def test_random_short_nodes(self):
        self.short_scenarios(20260923, 160)

    def test_scan_fallback_is_also_identical(self):
        """Views that fail their invariants use plain scans; force that path everywhere."""
        build = candidates.VisitView.__init__

        def unindexed(view, *args):
            build(view, *args)
            view.indexed = False
        with mock.patch.object(candidates.VisitView, "__init__", unindexed):
            self.short_scenarios(11, 40)

    def test_random_long_nodes_and_default_width(self):
        rng = random.Random(7)
        for case in range(8):
            sequences = {n: "".join(rng.choice(BASES) for _ in range(rng.randint(40, 160))) for n in (1, 2)}
            with self.subTest(case=case):
                self.check_scenario(rng, sequences, widths=(11, 101))

    def test_prebuilt_narrow_view_is_rebuilt_for_wide_windows(self):
        rng = random.Random(3)
        sequences = {1: "ACGTACGTAC", 2: "TTGCA"}
        reads = random_reads(rng, sequences)
        for candidate in probes(reads, sequences):
            narrow = NodeReads(candidate.node, reads, 1).eligible(candidate, 10)
            wide = [(r, s, v.visit) for r, s, v in narrow]
            x, m = make_tensor(candidate, narrow, {1: 5, 2: 5}, rows=4, width=11, debug=True)
            y, n = reference.make_tensor(candidate, wide, {1: 5, 2: 5}, rows=4, width=11, debug=True)
            np.testing.assert_array_equal(x, y)
            self.assertEqual(m, n)


if __name__ == "__main__":
    unittest.main()

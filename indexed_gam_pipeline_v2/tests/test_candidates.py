"""Edit decoding, candidate identity, support counting, window encoding and row selection."""
import unittest

import numpy as np

from fixtures import spec_alignment  # noqa: F401  (sets sys.path)
from indexed_gam_pipeline_v2.candidates import (BASES, OPS, Candidate, alt_support_bounds, decode_alignment,
                                                encode_count, make_tensor, overlap, rc)

COUNTS = {n: 90 for n in range(1, 10)}


def read(specs, sequences, **kwargs):
    return decode_alignment(spec_alignment(specs, sequences, **kwargs), sequences)[0]


def eligible(candidate, reads, min_bq=10):
    return [(r, *overlap(r, candidate, min_bq)) for r in reads if overlap(r, candidate, min_bq)]


def tensor(candidate, reads, **kwargs):
    kwargs.setdefault("path_counts", COUNTS)
    return make_tensor(candidate, eligible(candidate, reads), **kwargs)


class DecodeTest(unittest.TestCase):
    def test_n_filter_keeps_context_and_is_strand_independent(self):
        cases = [
            ("NAC", [(3, 3, "ANT")], [Candidate(1, 2, "C", "T", "SNP")]),
            ("AC", [(1, 1, ""), (0, 2, "GN"), (1, 1, "T")], [Candidate(1, 1, "C", "T", "SNP")]),
            ("ANC", [(1, 1, ""), (1, 0, ""), (1, 1, "T")], [Candidate(1, 2, "C", "T", "SNP")]),
            ("NC", [(1, 1, ""), (0, 1, "G"), (1, 1, "T")], [Candidate(1, 1, "C", "T", "SNP")]),
            ("NC", [(0, 1, "G"), (2, 2, "")], []),  # insertion anchored on N at node start
            ("AC", [(0, 1, "G"), (2, 2, "")], [Candidate(1, 0, "", "G", "INS")]),
            ("AN", [(1, 1, ""), (0, 1, "G"), (1, 1, "")], [Candidate(1, 1, "", "G", "INS")]),
            ("NA", [(2, 2, ""), (0, 1, "G")], [Candidate(1, 2, "", "G", "INS")]),
        ]
        for reference, edits, expected in cases:
            for reverse in (False, True):
                oriented = [(f, t, rc(s)) for f, t, s in reversed(edits)] if reverse else edits
                with self.subTest(reference=reference, reverse=reverse, edits=edits):
                    seq = {1: reference}
                    r, unsupported = decode_alignment(spec_alignment([(1, 0, reverse, oriented)], seq), seq)
                    self.assertEqual([o.candidate for o in r.observations], expected)
                    self.assertFalse(unsupported)
                    self.assertEqual(len(r.columns), sum(max(f, t) for f, t, _ in edits))
                    self.assertEqual((r.visits[0].start, r.visits[0].end), (0, len(reference)))

    def test_edits_indel_limits_and_unsupported_events(self):
        seq = {1: "A" * 110}
        a = spec_alignment([(1, 0, False, [(1, 1, "T"), (0, 50, "C" * 50), (50, 0, ""),
                                           (0, 51, "G" * 51), (51, 0, ""), (2, 3, "TGC")])], seq)
        r, unsupported = decode_alignment(a, seq)
        self.assertEqual([(o.candidate.kind, max(len(o.candidate.ref), len(o.candidate.alt)))
                          for o in r.observations], [("SNP", 1), ("INS", 50), ("DEL", 50)])
        self.assertEqual([e["reason"] for e in unsupported],
                         ["indel_exceeds_limit", "indel_exceeds_limit", "complex_replacement_not_supported"])
        for obs in r.observations:
            _, meta = tensor(obs.candidate, [r])
            self.assertGreaterEqual(meta["candidate_columns"][1] - meta["candidate_columns"][0],
                                    max(len(obs.candidate.ref), len(obs.candidate.alt)))
        # A GAM match edit is never rescanned into SNPs.
        a = spec_alignment([(1, 0, False, [(4, 4, "")])], seq)
        a.sequence = "TTTT"
        self.assertFalse(decode_alignment(a, seq)[0].observations)
        # Adjacent insertions merge before the length limit: 30I + 21I is one 51-bp event.
        seq = {1: "A" * 60}
        a = spec_alignment([(1, 0, False, [(1, 1, ""), (0, 30, "T" * 30), (0, 21, "G" * 21), (30, 0, ""), (21, 0, "")])], seq)
        r, rejected = decode_alignment(a, seq)
        self.assertFalse(r.observations)
        self.assertEqual([e["event_length"] for e in rejected], [51, 51])

    def test_target_nodes_restrict_observations_only(self):
        seq = {1: "ACGTAC", 2: "ACGTAC", 3: "A" * 60}
        for reverse in (False, True):
            a = spec_alignment([(1, 0, reverse, [(1, 1, "T"), (1, 1, ""), (0, 2, "GG"), (2, 0, ""), (2, 2, "")]),
                                (2, 0, reverse, [(1, 1, "T"), (1, 1, ""), (0, 2, "GG"), (2, 0, ""), (2, 2, "")]),
                                (1, 0, reverse, [(6, 6, "")]),
                                (3, 0, reverse, [(0, 51, "C" * 51), (60, 60, "")])], seq)
            full, rejected = decode_alignment(a, seq)
            for targets in (set(), {1}, {2}, {1, 2, 3}):
                new, new_rejected = decode_alignment(a, seq, target_nodes=targets)
                self.assertEqual((full.columns, full.visits, full.digest, rejected),
                                 (new.columns, new.visits, new.digest, new_rejected))
                self.assertEqual(new.observations, [o for o in full.observations if o.candidate.node in targets])
                for o in new.observations:
                    x, m = make_tensor(o.candidate, [(full, *overlap(full, o.candidate, 10))], COUNTS, rows=2, debug=True)
                    y, n = make_tensor(o.candidate, [(new, *overlap(new, o.candidate, 10))], COUNTS, rows=2, debug=True)
                    np.testing.assert_array_equal(x, y)
                    self.assertEqual(m, n)

    def test_invalid_alignments_are_rejected(self):
        seq = {1: "ACGT"}
        a = spec_alignment([(1, 0, False, [(4, 4, "")])], seq)
        a.quality = bytes([30])
        with self.assertRaisesRegex(ValueError, "quality length"):
            decode_alignment(a, seq)
        a = spec_alignment([(1, 0, False, [(4, 4, "")])], seq)
        a.sequence += "A"
        a.quality += bytes([30])
        with self.assertRaisesRegex(ValueError, "complete read"):
            decode_alignment(a, seq)
        a = spec_alignment([(1, 2, False, [(4, 4, "")])], seq)
        with self.assertRaisesRegex(ValueError, "edit bounds"):
            decode_alignment(a, seq)


class SupportTest(unittest.TestCase):
    def test_position_coverage_deletion_spans_and_row_limits(self):
        seq = {1: "ACGTAC"}
        c = Candidate(1, 2, "G", "T", "SNP")
        before = read([(1, 0, False, [(2, 2, "")])], seq)
        deletion = read([(1, 0, False, [(1, 1, ""), (4, 0, ""), (1, 1, "")])], seq)
        alt = read([(1, 0, False, [(2, 2, ""), (1, 1, "T"), (3, 3, "")])], seq)
        ref = read([(1, 0, False, [(6, 6, "")])], seq)
        e = eligible(c, [before, deletion, alt, alt, ref])
        self.assertEqual(len(e), 4)
        _, small = make_tensor(c, e, COUNTS, rows=1)
        _, large = make_tensor(c, e, COUNTS, rows=10)
        for k in ("coverage", "af", "alt_count", "ref_count", "other_count"):
            self.assertEqual(small[k], large[k])
        self.assertEqual((small["af"], small["selected_alignments"], large["selected_alignments"]), (.5, 1, 4))
        first, _ = make_tensor(c, e, COUNTS, rows=3)
        shuffled, _ = make_tensor(c, list(reversed(e)), COUNTS, rows=3)
        np.testing.assert_array_equal(first, shuffled)
        self.assertEqual(overlap(deletion, Candidate(1, 2, "", "T", "INS"), 10)[0], "other")
        self.assertEqual(overlap(before, Candidate(1, 1, "CGTA", "", "DEL"), 10)[0], "other")

    def test_insertion_boundaries_and_node_edges(self):
        seq = {1: "AC", 2: "GT", 3: "TA"}
        c = Candidate(1, 2, "", "T", "INS")
        crossing = read([(1, 0, False, [(2, 2, "")]), (2, 0, False, [(2, 2, "")])], seq)
        terminal = read([(1, 0, False, [(2, 2, "")])], seq)
        self.assertEqual(overlap(crossing, c, 10)[0], "ref")
        self.assertEqual(overlap(terminal, c, 10)[0], "other")
        partial = read([(1, 0, False, [(1, 1, "")]), (2, 0, False, [(2, 2, "")])], seq)
        self.assertEqual(overlap(partial, Candidate(1, 1, "", "T", "INS"), 10)[0], "other")
        self.assertEqual(overlap(crossing, Candidate(2, 0, "", "T", "INS"), 10)[0], "ref")
        c = Candidate(1, 1, "C", "T", "SNP")
        a = read([(1, 0, False, [(1, 1, ""), (1, 1, "T")]), (2, 0, False, [(1, 1, ""), (0, 1, "A"), (1, 1, "")])], seq)
        b = read([(1, 0, False, [(2, 2, "")]), (3, 0, False, [(2, 2, "")])], seq)
        x, m = tensor(c, [a, b], width=9)
        lo = m["candidate_columns"][0]
        self.assertEqual(x[5, 0, lo + 1:lo + 4].tolist(), [BASES[v] for v in "G-T"])
        self.assertEqual(x[5, 1, lo + 1:lo + 3].tolist(), [BASES[v] for v in "TA"])

    def test_repeated_visits_count_once_and_keep_path(self):
        seq = {1: "ACGT", 2: "T"}
        r = read([(1, 0, False, [(4, 4, "")]), (2, 0, False, [(1, 1, "")]),
                  (1, 0, False, [(1, 1, ""), (1, 1, "T"), (2, 2, "")])], seq)
        c = r.observations[0].candidate
        _, m = tensor(c, [r, r], debug=True)
        self.assertEqual((m["coverage"], m["alt_count"]), (2, 2))
        self.assertEqual(m["rows"][0]["anchor_mapping_index"], 2)
        self.assertEqual([n["node_id"] for n in m["row_groups"][0]["path"]], [1, 2, 1])

    def test_low_quality_alt_is_other(self):
        seq = {1: "C"}
        a = spec_alignment([(1, 0, True, [(1, 1, "T")])], seq)
        a.quality = bytes([3])
        r, _ = decode_alignment(a, seq)
        c = r.observations[0].candidate
        self.assertEqual(c, Candidate(1, 0, "C", "A", "SNP"))
        self.assertEqual(overlap(r, c, 10)[0], "other")
        self.assertEqual(overlap(r, c, 0)[0], "alt")
        _, m = tensor(c, [r])
        self.assertEqual((m["coverage"], m["other_count"]), (1, 1))

    def test_alt_support_upper_bound(self):
        seq = {1: "ACGT", 2: "TT"}
        a = spec_alignment([(1, 0, False, [(1, 1, ""), (1, 1, "T"), (2, 2, "")]), (2, 0, False, [(2, 2, "")]),
                            (1, 0, False, [(1, 1, ""), (1, 1, "T"), (2, 2, "")])], seq)
        b = spec_alignment([(1, 0, True, [(2, 2, ""), (1, 1, "A"), (1, 1, "")])], seq, name="reverse")
        b.quality = bytes([2] * 4)
        c = spec_alignment([(1, 0, False, [(1, 1, ""), (2, 0, ""), (1, 1, "")])], seq, name="del")
        d = spec_alignment([(1, 0, False, [(2, 2, ""), (0, 2, "TT"), (2, 2, "")])], seq, name="ins")
        reads = [decode_alignment(x, seq)[0] for x in (a, a, b, c, d)]
        candidates = {o.candidate for r in reads for o in r.observations}
        counts = alt_support_bounds(reads, candidates, 10)
        snp = Candidate(1, 1, "C", "T", "SNP")
        self.assertEqual(counts[snp], 2)  # duplicate records count twice, repeated visits once
        self.assertEqual(alt_support_bounds(reads, candidates, 0)[snp], 3)
        for candidate in candidates:
            actual = sum(overlap(r, candidate, 10)[0] == "alt" for r in reads if overlap(r, candidate, 10))
            self.assertLessEqual(actual, counts[candidate])


class WindowTest(unittest.TestCase):
    def test_default_101_columns_has_50_flanking_columns(self):
        seq = {1: "A" * 200}
        r = read([(1, 0, False, [(100, 100, ""), (1, 1, "T"), (99, 99, "")])], seq)
        c = r.observations[0].candidate
        x, m = tensor(c, [r], debug=True)
        self.assertEqual(x.shape, (8, 200, 101))
        self.assertEqual((m["anchor_column"], m["candidate_columns"]), (50, [50, 51]))
        self.assertEqual([p["offset"] for p in m["rows"][0]["columns"]], list(range(50, 151)))

    def test_long_deletion_and_background_insertion(self):
        seq = {1: "A" * 200}
        candidate = Candidate(1, 60, "A" * 40, "", "DEL")
        alt = read([(1, 0, False, [(60, 60, ""), (40, 0, ""), (100, 100, "")])], seq)
        other = read([(1, 0, False, [(65, 65, ""), (0, 150, "T" * 150), (135, 135, "")])], seq)
        x, m = tensor(candidate, [alt, other], rows=3, width=100, debug=True)
        self.assertEqual((x.shape, m["anchor_column"], m["coverage"], m["alt_count"]), ((8, 3, 100), 50, 2, 1))
        for ri, row in enumerate(m["rows"]):
            self.assertEqual(row["columns"][50]["offset"], 60)
            if row["record_sha256"] == other.digest:
                self.assertTrue(np.all(x[4, ri, 55:] == OPS["I"]))
                self.assertEqual(row["window_mismatch_bp"], 45)
                self.assertGreater(row["omitted_central_context_columns"], 0)
            else:
                self.assertTrue(np.all(x[4, ri, 50:90] == OPS["D"]))
                self.assertTrue(np.all(x[0, ri, 50:90] == BASES["-"]))
        self.assertFalse(x[:, 2].any())

    def test_long_insertion_reverse_equivalence_and_partial_coverage(self):
        seq = {1: "A" * 120}
        f = read([(1, 0, False, [(60, 60, ""), (0, 150, "C" * 150), (60, 60, "")])], seq)
        r = read([(1, 0, True, [(60, 60, ""), (0, 150, "G" * 150), (60, 60, "")])], seq)
        candidate = Candidate(1, 60, "", "C" * 30, "INS")
        xf, m = tensor(candidate, [f], width=100)
        xr, _ = tensor(candidate, [r], width=100)
        np.testing.assert_array_equal(xf[:7], xr[:7])  # everything but strand is orientation-invariant
        self.assertTrue(np.all(xf[7, 0][xf[0, 0] != 0] == 1) and np.all(xr[7, 0][xr[0, 0] != 0] == 2))
        self.assertTrue(np.all(xf[4, 0, 50:] == OPS["I"]))
        self.assertEqual(m["window_mismatch_bp"], [50])
        partial = read([(1, 65, False, [(20, 20, "")])], seq)
        x, m = tensor(Candidate(1, 60, "A" * 20, "", "DEL"), [partial], width=100, debug=True)
        self.assertFalse(x[:, 0, 50:55].any())
        self.assertEqual(m["rows"][0]["columns"][55]["offset"], 65)

    def test_window_edit_bp_counts_only_visible_edits(self):
        seq = {1: "A" * 60}
        candidate = Candidate(1, 30, "A", "T", "SNP")
        specs = [[(10, 10, "T" * 10), (20, 20, ""), (1, 1, "T"), (29, 29, "")],
                 [(28, 28, ""), (0, 2, "GG"), (2, 2, ""), (1, 1, "T"), (29, 29, "")],
                 [(28, 28, ""), (2, 0, ""), (1, 1, "T"), (29, 29, "")],
                 [(29, 29, ""), (1, 1, "C"), (1, 1, "T"), (29, 29, "")]]
        reads = [read([(1, 0, False, edits)], seq, name=str(i)) for i, edits in enumerate(specs)]
        _, meta = tensor(candidate, reads, width=11, debug=True)
        self.assertEqual({row["read_name"]: row["window_mismatch_bp"] for row in meta["rows"]},
                         {"0": 1, "1": 3, "2": 3, "3": 2})
        self.assertEqual(meta["window_mismatch_bp"], [3, 3, 2, 1])

    def test_group_before_uniform_sampling(self):
        seq = {1: "AC", 2: "GG", 3: "TT"}
        reads = []
        for i in range(401):
            branch = 3 if i % 2 else 2
            first = [(1, 1, "A"), (1, 1, "")] if branch == 3 else [(2, 2, "")]
            reads.append(read([(branch, 0, False, first), (1, 0, False, [(1, 1, ""), (1, 1, "T")])], seq, name=f"read-{i:04d}"))
        c = Candidate(1, 1, "C", "T", "SNP")
        e = eligible(c, reads)
        counts = {1: 90, 2: 3, 3: 7}
        full, all_meta = make_tensor(c, e, counts, rows=401, width=9, debug=True)
        x, meta = make_tensor(c, list(reversed(e)), counts, rows=200, width=9, debug=True)
        indices = [i * 400 // 199 for i in range(200)]
        np.testing.assert_array_equal(x, full[:, indices, :])
        self.assertEqual(meta["selected_grouped_ranks"], indices)
        self.assertEqual([g["path"][0]["node_id"] for g in meta["row_groups"]], [3, 2])
        self.assertEqual((meta["coverage"], meta["alt_count"], meta["selected_alignments"]), (401, 401, 200))
        self.assertEqual([r["record_sha256"] for r in meta["rows"]],
                         [all_meta["rows"][i]["record_sha256"] for i in indices])

    def test_multinode_branches_and_per_row_reference(self):
        seq = {1: "AACG", 2: "TT", 3: "GC", 4: "AT"}
        a = read([(2, 0, False, [(2, 2, "")]), (1, 1, False, [(1, 1, ""), (1, 1, "T"), (1, 1, "")]), (4, 0, False, [(2, 2, "")])], seq)
        b = read([(3, 0, False, [(2, 2, "")]), (1, 1, False, [(3, 3, "")]), (4, 0, False, [(2, 2, "")])], seq)
        c = a.observations[0].candidate
        x, meta = tensor(c, [a, b], rows=3, width=9, debug=True)
        self.assertEqual((meta["coverage"], meta["alt_count"], meta["ref_count"]), (2, 1, 1))
        self.assertEqual(x[5, 0, 1:8].tolist(), [BASES[v] for v in "TTACGAT"])
        self.assertEqual(x[5, 1, 1:8].tolist(), [BASES[v] for v in "GCACGAT"])
        self.assertFalse(x[:, 2].any())
        self.assertEqual([c["node_id"] for c in meta["rows"][0]["columns"][1:8]], [2, 2, 1, 1, 1, 4, 4])

    def test_reverse_equivalence_and_quality_orientation(self):
        seq = {1: "AACG", 2: "TC"}
        forward = read([(2, 0, False, [(2, 2, "")]), (1, 0, False, [(1, 1, ""), (1, 1, "T"), (2, 2, "")])], seq)
        a = spec_alignment([(1, 0, True, [(2, 2, ""), (1, 1, "A"), (1, 1, "")]), (2, 0, True, [(2, 2, "")])], seq)
        a.quality = bytes([10, 11, 12, 13, 14, 15])
        reverse = decode_alignment(a, seq)[0]
        c = forward.observations[0].candidate
        self.assertEqual(c, reverse.observations[0].candidate)
        x, meta = tensor(c, [reverse], width=9, debug=True)
        start = meta["candidate_columns"][0]
        self.assertEqual((x[0, 0, start], x[1, 0, start]), (BASES["T"], 12))
        self.assertEqual(x[5, 0, start - 3:start + 3].tolist(), [BASES[v] for v in "TCAACG"])
        self.assertTrue(meta["rows"][0]["reversed_for_candidate"])
        for edits, rev_edits in (([(2, 2, ""), (0, 2, "TG"), (2, 2, "")], [(2, 2, ""), (0, 2, "CA"), (2, 2, "")]),
                                 ([(1, 1, ""), (2, 0, ""), (1, 1, "")], [(1, 1, ""), (2, 0, ""), (1, 1, "")])):
            f = read([(1, 0, False, edits)], seq)
            r = read([(1, 0, True, rev_edits)], seq)
            self.assertEqual(f.observations[0].candidate, r.observations[0].candidate)
            xf, _ = tensor(f.observations[0].candidate, [f], width=9)
            xr, _ = tensor(r.observations[0].candidate, [r], width=9)
            np.testing.assert_array_equal(xf[:7], xr[:7])
            self.assertEqual((int(xf[7, 0, 4]), int(xr[7, 0, 4])), (1, 2))

    def test_insertion_gap_padding_and_alleles(self):
        seq = {1: "ACGT"}
        a = read([(1, 0, False, [(2, 2, ""), (0, 2, "TA"), (2, 2, "")])], seq)
        b = read([(1, 0, False, [(4, 4, "")])], seq)
        d = read([(1, 0, False, [(2, 2, ""), (0, 1, "G"), (2, 2, "")])], seq)
        terminal = read([(1, 0, False, [(2, 2, "")])], seq)
        c = a.observations[0].candidate
        x, m = tensor(c, [a, b, d, terminal], width=8, debug=True)
        lo, hi = m["candidate_columns"]
        self.assertEqual((m["coverage"], m["alt_count"], m["ref_count"], m["other_count"], hi - lo), (4, 1, 1, 2, 2))
        self.assertEqual(x[0, 0, lo:hi].tolist(), [BASES["T"], BASES["A"]])
        for row, detail in enumerate(m["rows"]):
            if detail["record_sha256"] == b.digest:
                self.assertEqual(x[0, row, lo:hi].tolist(), [6, 6])
                self.assertEqual(x[1, row, lo], -1)
            if detail["record_sha256"] == terminal.digest:
                self.assertFalse(x[:, row, lo:hi].any())
            else:
                self.assertTrue(np.all(x[5, row, lo:hi] == 6))
                self.assertEqual(x[2, row, lo:hi].tolist(), [BASES["T"], BASES["A"]])
            if detail["record_sha256"] == d.digest:
                self.assertEqual(x[0, row, lo:hi].tolist(), [BASES["G"], 6])

    def test_central_context_insertions_within_deletion(self):
        seq = {1: "ACGTAC"}
        alt = read([(1, 0, False, [(1, 1, ""), (4, 0, ""), (1, 1, "")])], seq)
        other = read([(1, 0, False, [(2, 2, ""), (0, 2, "TT"), (4, 4, "")])], seq)
        c = alt.observations[0].candidate
        x, m = tensor(c, [alt, other], width=12, debug=True)
        lo, hi = m["candidate_columns"]
        self.assertEqual(hi - lo, 4)
        self.assertEqual(x[0, 1, lo:lo + 6].tolist(), [BASES[b] for b in "CTTGTA"])
        self.assertEqual(x[5, 1, lo:lo + 6].tolist(), [BASES[b] for b in "C--GTA"])
        self.assertEqual(m["other_count"], 1)

    def test_node_path_grouping_preserves_counts_and_channels(self):
        seq = {1: "AC", 2: "GG", 3: "TT"}
        rows = []
        for branch, base in ((3, "T"), (2, "C"), (2, "T"), (3, "C")):
            edits = [(1, 1, ""), (1, 1, "T")] if base == "T" else [(2, 2, "")]
            rows.append(read([(branch, 0, False, [(2, 2, "")]), (1, 0, False, edits)], seq))
        c = Candidate(1, 1, "C", "T", "SNP")
        e = eligible(c, rows)
        x, m = make_tensor(c, e, COUNTS, rows=6, width=9, debug=True)
        self.assertEqual([r["support"] for r in m["rows"]], ["alt", "ref", "alt", "ref"])
        self.assertEqual([(g["start_row"], g["end_row"]) for g in m["row_groups"]], [(0, 2), (2, 4)])
        self.assertEqual((m["coverage"], m["alt_count"], m["ref_count"], m["af"]), (4, 2, 2, .5))
        self.assertFalse(x[:, 4:].any())
        for ri, detail in enumerate(m["rows"]):
            original = next(item for item in e if item[0].digest == detail["record_sha256"])
            single, meta = make_tensor(c, [original], COUNTS, rows=1, width=9, debug=True)
            np.testing.assert_array_equal(x[:, ri], single[:, 0])
            self.assertEqual(detail["columns"], meta["rows"][0]["columns"])
        _, small = make_tensor(c, e, COUNTS, rows=2, width=9, debug=True)
        self.assertEqual((small["selected_counts"], small["selected_grouped_ranks"], small["af"]), ({"alt": 1, "ref": 1}, [0, 3], .5))
        shuffled, _ = make_tensor(c, list(reversed(e)), COUNTS, rows=6, width=9)
        np.testing.assert_array_equal(x, shuffled)

    def test_group_key_ignores_distant_branches_and_normalizes_strand(self):
        seq = {1: "AC", 2: "G" * 20, 3: "TT", 4: "AA"}
        a = read([(3, 0, False, [(2, 2, "")]), (2, 0, False, [(20, 20, "")]), (1, 0, False, [(1, 1, ""), (1, 1, "T")])], seq)
        b = read([(1, 0, True, [(1, 1, "A"), (1, 1, "")]), (2, 0, True, [(20, 20, "")]), (4, 0, True, [(2, 2, "")])], seq)
        c = a.observations[0].candidate
        _, m = tensor(c, [a, b], width=9, debug=True)
        self.assertEqual(m["row_groups"], [dict(start_row=0, end_row=2, path=[dict(node_id=2, reverse=False), dict(node_id=1, reverse=False)])])


class CandidateAltAndStrandTest(unittest.TestCase):
    def test_snp_stripe_encodes_the_question_and_strand_the_read_orientation(self):
        seq = {1: "ACGTAC"}
        c = Candidate(1, 2, "G", "T", "SNP")
        alt = read([(1, 0, False, [(2, 2, ""), (1, 1, "T"), (3, 3, "")])], seq, name="alt")
        ref = read([(1, 0, False, [(6, 6, "")])], seq, name="ref")
        other = read([(1, 0, False, [(2, 2, ""), (1, 1, "A"), (3, 3, "")])], seq, name="other")
        # Reverse-strand read carrying the same G>T: on its own strand it reads A instead of C.
        rev = read([(1, 0, True, [(3, 3, ""), (1, 1, "A"), (2, 2, "")])], seq, name="rev")
        self.assertEqual(rev.observations[0].candidate, c)
        x, m = tensor(c, [alt, ref, other, rev], rows=6, width=9, debug=True)
        anchor = m["anchor_column"]
        self.assertEqual(m["candidate_columns"], [anchor, anchor + 1])
        by_name = {row["read_name"]: ri for ri, row in enumerate(m["rows"])}
        for name, ri in by_name.items():
            stripe = [0] * 9
            stripe[anchor] = BASES["T"]
            self.assertEqual(x[2, ri].tolist(), stripe, name)  # identical question in every row
            evidence = x[0, ri] != 0
            self.assertTrue(np.all(x[7, ri][evidence] == (2 if name == "rev" else 1)), name)
            self.assertFalse(x[7, ri][~evidence].any(), name)
        matches_alt = {name: bool(x[0, ri, anchor] == x[2, ri, anchor]) for name, ri in by_name.items()}
        matches_ref = {name: bool(x[0, ri, anchor] == x[5, ri, anchor]) for name, ri in by_name.items()}
        self.assertEqual(matches_alt, dict(alt=True, rev=True, ref=False, other=False))
        self.assertEqual(matches_ref, dict(alt=False, rev=False, ref=True, other=False))
        self.assertEqual({row["read_name"]: row["support"] for row in m["rows"]},
                         dict(alt="alt", rev="alt", ref="ref", other="other"))
        self.assertFalse(x[:, 4:].any())  # unused rows stay zero in all eight channels

    def test_insertion_and_deletion_stripes(self):
        seq = {1: "ACGT"}
        a = read([(1, 0, False, [(2, 2, ""), (0, 2, "TA"), (2, 2, "")])], seq, name="ins")
        b = read([(1, 0, False, [(4, 4, "")])], seq, name="ref")
        d = read([(1, 0, False, [(2, 2, ""), (0, 1, "G"), (2, 2, "")])], seq, name="other-ins")
        terminal = read([(1, 0, False, [(2, 2, "")])], seq, name="terminal")
        c = a.observations[0].candidate
        self.assertEqual(c, Candidate(1, 2, "", "TA", "INS"))
        x, m = tensor(c, [a, b, d, terminal], width=8, debug=True)
        lo, hi = m["candidate_columns"]
        self.assertEqual(hi - lo, 2)
        for ri, row in enumerate(m["rows"]):
            self.assertFalse(x[2, ri, :lo].any() or x[2, ri, hi:].any(), row["read_name"])
            if row["read_name"] == "terminal":  # no boundary evidence: no cells, so no stripe
                self.assertFalse(x[:, ri, lo:hi].any())
                continue
            self.assertEqual(x[2, ri, lo:hi].tolist(), [BASES["T"], BASES["A"]], row["read_name"])
            exact = bool(np.all(x[0, ri, lo:hi] == x[2, ri, lo:hi]))
            self.assertEqual(exact, row["read_name"] == "ins", row["read_name"])
        seq = {1: "ACGTAC"}
        alt = read([(1, 0, False, [(1, 1, ""), (4, 0, ""), (1, 1, "")])], seq, name="del")
        ref = read([(1, 0, False, [(6, 6, "")])], seq, name="ref")
        c = alt.observations[0].candidate
        self.assertEqual(c, Candidate(1, 1, "CGTA", "", "DEL"))
        x, m = tensor(c, [alt, ref], width=12, debug=True)
        lo, hi = m["candidate_columns"]
        self.assertEqual(hi - lo, 4)
        for ri, row in enumerate(m["rows"]):
            self.assertEqual(x[2, ri, lo:hi].tolist(), [6, 6, 6, 6])
            self.assertEqual(x[5, ri, lo:hi].tolist(), [BASES[b] for b in "CGTA"])  # REF stays readable in channel 5
            self.assertEqual(bool(np.all(x[0, ri, lo:hi] == x[2, ri, lo:hi])), row["read_name"] == "del")
        # A long deletion is cropped at the window edge together with its stripe.
        seq = {1: "A" * 200}
        c = Candidate(1, 60, "A" * 60, "", "DEL")
        r = decode_alignment(spec_alignment([(1, 0, False, [(60, 60, ""), (60, 0, ""), (80, 80, "")])], seq), seq, max_indel=60)[0]
        x, m = tensor(c, [r], width=100)
        self.assertEqual(m["candidate_columns"], [50, 100])
        self.assertTrue(np.all(x[2, 0, 50:] == 6) and not x[2, 0, :50].any())


class StorageTest(unittest.TestCase):
    def test_path_count_log_encoding_and_padding(self):
        seq = {1: "ACA"}
        r = read([(1, 0, False, [(1, 1, ""), (1, 1, "T"), (1, 1, "")])], seq)
        c = r.observations[0].candidate
        expected = {0: 0, 1: 14, 2: 22, 3: 28, 4: 33, 5: 36, 7: 42, 8: 44, 15: 56, 16: 57, 45: 77,
                    89: 91, 90: 91, 127: 98, 331: 117, 523: 126, 524: 127, 550: 127, 40000: 127}
        for count, code in expected.items():
            with self.subTest(count=count):
                self.assertEqual(encode_count(count), code)
                x, m = tensor(c, [r], rows=2, width=5, path_counts={1: count})
                self.assertEqual((x.dtype, x.nbytes), (np.int8, 8 * 2 * 5))
                np.testing.assert_array_equal(x[6, 0], [0, code, code, code, 0])
                self.assertFalse(x[:, 1].any())
                self.assertEqual((m["coverage"], m["alt_count"], m["af"]), (1, 1, 1.0))
        self.assertEqual(len({encode_count(c) for c in range(1, 20)}), 19)  # small counts stay distinct
        for count in range(1, 3000):
            self.assertLessEqual(encode_count(count), encode_count(count + 1))  # monotonic
        with self.assertRaisesRegex(ValueError, "int32"):
            tensor(c, [r], path_counts={1: -1})

    def test_quality_and_mapq_clip_without_changing_raw_values(self):
        seq = {1: "ACA"}
        for quality, expected in ((b"", -1), (bytes([0] * 3), 0), (bytes([127] * 3), 127), (bytes([255] * 3), 127)):
            for mapq in (0, 127, 255, 40000):
                with self.subTest(quality=quality, mapq=mapq):
                    a = spec_alignment([(1, 0, False, [(1, 1, ""), (1, 1, "T"), (1, 1, "")])], seq)
                    a.quality, a.mapping_quality = quality, mapq
                    r, _ = decode_alignment(a, seq)
                    x, _ = make_tensor(r.observations[0].candidate, eligible(r.observations[0].candidate, [r], 0),
                                       {1: 331}, rows=1, width=3)
                    self.assertTrue(np.all(x[1] == expected))
                    self.assertTrue(np.all(x[3] == min(mapq, 127)))
                    self.assertEqual(r.mapq, mapq)


if __name__ == "__main__":
    unittest.main()

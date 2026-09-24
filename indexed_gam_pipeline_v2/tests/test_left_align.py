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
        left_end = [(1, 0, False, [(1, 1, "")]), (2, 0, False, [(0, 2, "AT"), (2, 2, "")]), (3, 0, False, [(2, 2, "")]),
                    (4, 0, False, [(1, 1, "")])]
        for strand_specs, moved in ((specs, [(3, 2)]), (mirror(specs, sequences), [(3, 2)]),
                                    (mirror(left_end, sequences), [])):
            read = decode(strand_specs, sequences)
            self.assertEqual(candidates(read), [Candidate(2, 0, "", "AT", "INS")])
            self.assertEqual(read.moves, moved)  # placed on node 3 -> moved to node 2; placed on node 2 -> stays
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
        # deletions of any length cross node boundaries
        read = decode([(2, 0, False, [(2, 2, "")]), (3, 0, False, [(2, 0, ""), (1, 1, "")])], sequences)
        self.assertEqual(candidates(read), [Candidate(2, 0, "AA", "", "DEL")])
        read = decode([(2, 0, False, [(2, 2, "")]), (3, 0, False, [(1, 1, ""), (1, 0, ""), (1, 1, "")])], sequences)
        self.assertEqual(candidates(read), [Candidate(2, 0, "A", "", "DEL")])
        # an orientation change stops it
        read = decode([(2, 0, True, [(2, 2, "")]), (3, 0, False, [(0, 1, "A"), (3, 3, "")])], sequences)
        self.assertEqual(candidates(read), [Candidate(3, 0, "", "A", "INS")])
        # N bases are never moved, nor moved across
        read = decode([(1, 0, False, [(5, 5, ""), (0, 1, "N"), (1, 1, "")])], sequences)
        self.assertEqual(candidates(read), [])
        self.assertEqual([c.pos for c in read.columns if c.boundary], [5])

    def test_indel_quality_is_strand_symmetric_under_read_direction_quality_patterns(self):
        # HiFi-like pattern: the last base of a homopolymer in sequencing direction has a low QV.
        sequences = {1: "GAAAAC"}
        forward = [(1, 0, False, [(5, 5, ""), (0, 1, "A"), (1, 1, "")])]  # G AAAA +A C
        deletion = [(1, 0, False, [(4, 4, ""), (1, 0, ""), (1, 1, "")])]  # G AAA -A C
        for specs, low in ((forward, 5), (deletion, 4)):
            reads = {}
            for strand, strand_specs in (("forward", specs), ("reverse", mirror(specs, sequences))):
                a = spec_alignment(strand_specs, sequences)
                quality = [40] * len(a.sequence)
                quality[low] = 5  # read index `low` = the run's last base in read order
                a.quality = bytes(quality)
                reads[strand] = decode_alignment(a, sequences)[0]
            (f,), (r,) = reads["forward"].observations, reads["reverse"].observations
            self.assertEqual(f.candidate, r.candidate)
            self.assertEqual(f.quality, r.quality)  # both see the whole repeat, low base included
        # outside a repeat the quality is the plain one: inserted bases / nearest flanks
        plain = decode_alignment(spec_alignment([(1, 0, False, [(1, 1, ""), (0, 1, "T"), (5, 5, "")])], sequences),
                                 sequences)[0]
        self.assertEqual([o.quality for o in plain.observations], [30])

    def test_deletion_across_nodes_is_one_candidate_with_its_path(self):
        from indexed_gam_pipeline_v2.candidates import make_tensor, NodeReads
        sequences = {1: "GCA", 2: "TG", 3: "GAT"}
        # vg writes it as three edits: A at the end of node 1, all of node 2, G at the start of node 3
        spec = [(1, 0, False, [(2, 2, ""), (1, 0, "")]), (2, 0, False, [(2, 0, "")]), (3, 0, False, [(1, 0, ""), (2, 2, "")])]
        expected = Candidate(1, 2, "ATGG", "", "DEL", ((2, 0, 2), (3, 0, 1)))
        self.assertEqual(expected.metadata()["candidate_id"], "1:2:DEL:ATGG>@2+3")
        self.assertEqual(expected.positions(), [(1, 2), (2, 0), (2, 1), (3, 0)])
        for strand in (spec, mirror(spec, sequences)):
            read = decode(strand, sequences)
            self.assertEqual(candidates(read), [expected])
            self.assertEqual(overlap(read, expected, 10)[0], "alt")
        ref = decode([(1, 0, False, [(3, 3, "")]), (2, 0, False, [(2, 2, "")]), (3, 0, False, [(3, 3, "")])], sequences)
        self.assertEqual(overlap(ref, expected, 10)[0], "ref")  # matched over all three nodes, aligned neighbours
        branch = decode([(1, 0, False, [(3, 3, "")]), (3, 0, False, [(3, 3, "")])], {**sequences})
        self.assertEqual(overlap(branch, expected, 10)[0], "other")  # skips node 2: not the reference path
        alt = decode(spec, sequences)
        eligible = NodeReads(1, [alt, ref], 11).eligible(expected, 10)
        x, m = make_tensor(expected, eligible, {1: 5, 2: 5, 3: 5}, rows=2, width=11, debug=True)
        self.assertEqual((m["candidate_columns"], m["site_layout"]["reference"], m["path"]), ([5, 9], "ATGG", [[2, 0, 2], [3, 0, 1]]))
        rows = {r["site_allele"]: i for i, r in enumerate(m["rows"])}
        self.assertEqual(x[2, rows["A1"], 5:9].tolist(), [6, 6, 6, 6])
        self.assertEqual(x[2, rows["REF"], 5:9].tolist(), [1, 4, 3, 3])  # A T G G

    def test_deletion_longer_than_the_limit_only_in_total_is_not_a_candidate(self):
        sequences = {1: "GCA", 2: "TG", 3: "GAT"}
        spec = [(1, 0, False, [(2, 2, ""), (1, 0, "")]), (2, 0, False, [(2, 0, "")]), (3, 0, False, [(1, 0, ""), (2, 2, "")])]
        read, unsupported = decode_alignment(spec_alignment(spec, sequences), sequences, max_indel=3)
        self.assertEqual(candidates(read), [])  # pieces of 1, 2 and 1 bp, 4 bp in total
        self.assertEqual([(u["reason"], u["event_length"], u["mappings"]) for u in unsupported],
                         [("indel_exceeds_limit_across_mappings", 4, 3)])

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
                for s, e, op in indel_runs(read.columns, join=True):
                    rev = read.columns[s].reverse
                    k = e if rev else s - 1
                    if not 0 <= k < len(read.columns):
                        continue
                    nb, edge = read.columns[k], read.columns[s] if rev else read.columns[e - 1]
                    movable = (nb.op == "M" and not nb.boundary and nb.reverse == rev and nb.read in "ACGT"
                               and (nb.read == edge.read if op == "I" else nb.ref == edge.ref))
                    self.assertFalse(movable and "N" not in "".join(c.read + c.ref for c in read.columns[s:e]),
                                     (s, e, op))
                # the other strand's report of the same molecule gives the same candidates and columns,
                # and every observation the same quality (qualities travel with their bases)
                other = decode(mirror(specs, sequences), sequences)
                self.assertEqual(candidates(other), candidates(read))
                self.assertEqual(mirrored_columns(other), columns_of(read))
                a.quality = bytes(rng.randint(0, 60) for _ in a.sequence)
                b = spec_alignment(mirror(specs, sequences), sequences)
                b.quality = a.quality[::-1]
                observed = lambda x: sorted((o.candidate, o.quality) for o in decode_alignment(x, sequences)[0].observations)
                self.assertEqual(observed(b), observed(a))


class SupplementRunTest(unittest.TestCase):
    """An indel normalized off the target nodes is listed, and a supplement run over the list builds it."""

    def test_displaced_indel_is_listed_and_recovered_with_both_strands(self):
        from contextlib import redirect_stdout
        import io
        import json
        from pathlib import Path
        import tempfile
        from fixtures import build_args, graph_fixture, write_gam
        from indexed_gam_pipeline_v2.build import build
        from indexed_gam_pipeline_v2 import orchestrate
        sequences = {1: "C", 2: "AT", 3: "AT", 4: "G"}
        # vg-like placement: forward reads at the repeat's right end (node 3), reverse reads at its
        # left end (node 2 in forward coordinates).
        forward = [(1, 0, False, [(1, 1, "")]), (2, 0, False, [(2, 2, "")]),
                   (3, 0, False, [(2, 2, ""), (0, 2, "AT")]), (4, 0, False, [(1, 1, "")])]
        left_end = [(1, 0, False, [(1, 1, "")]), (2, 0, False, [(0, 2, "AT"), (2, 2, "")]),
                    (3, 0, False, [(2, 2, "")]), (4, 0, False, [(1, 1, "")])]
        ref = [(n, 0, False, [(len(sequences[n]), len(sequences[n]), "")]) for n in (1, 2, 3, 4)]
        rows = [spec_alignment(forward, sequences, name=f"f{i}") for i in range(3)]
        rows += [spec_alignment(mirror(left_end, sequences), sequences, name=f"r{i}") for i in range(3)]
        rows += [spec_alignment(ref, sequences, name=f"ref{i}") for i in range(4)]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam = write_gam(root / "g.gam", rows)
            graph = graph_fixture(root / "graph.sqlite", [(n, s, 5) for n, s in sequences.items()])
            (root / "main.txt").write_text("3\n")  # vg's forward placement made only node 3 a target
            common = dict(gam=str(gam), graph_index=str(graph), min_variants=1, min_af=.05, rows=12, width=11,
                          max_indel_len=5, batch_nodes=4)
            with redirect_stdout(io.StringIO()):
                build(build_args(**common, nodes=str(root / "main.txt"), output=str(root / "main")))
            self.assertEqual((root / "main/displaced_nodes.tsv").read_text(), "2\t3\n")
            self.assertFalse((root / "main/variant_summary.ndjson").read_text())  # the site is not on node 3
            (root / "supplement.txt").write_text("2\n")
            with redirect_stdout(io.StringIO()):
                build(build_args(**common, nodes=str(root / "supplement.txt"), output=str(root / "supplement")))
            (site,) = [json.loads(line) for line in (root / "supplement/variant_summary.ndjson").read_text().splitlines()]
            self.assertEqual((site["site_id"], site["alt_count"], site["coverage"]), ("2:0:INDEL", 6, 10))
            self.assertEqual(site["site_counts"], {"A1": 6, "REF": 4})
            self.assertEqual((root / "supplement/displaced_nodes.tsv").read_text(), "")

    def test_multi_node_deletion_build_passes_the_independent_audit(self):
        from contextlib import redirect_stdout
        import io
        import json
        from pathlib import Path
        import tempfile
        from fixtures import build_args, graph_fixture, write_gam
        from indexed_gam_pipeline_v2.build import build
        from indexed_gam_pipeline_v2.validate_examples import validate
        sequences = {1: "GCA", 2: "TG", 3: "GAT"}
        spec = [(1, 0, False, [(2, 2, ""), (1, 0, "")]), (2, 0, False, [(2, 0, "")]), (3, 0, False, [(1, 0, ""), (2, 2, "")])]
        ref = [(n, 0, False, [(len(s), len(s), "")]) for n, s in sequences.items()]
        rows = [spec_alignment(spec, sequences, name=f"a{i}") for i in range(3)]
        rows += [spec_alignment(mirror(spec, sequences), sequences, name=f"b{i}") for i in range(2)]
        rows += [spec_alignment(ref, sequences, name=f"r{i}") for i in range(4)]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam = write_gam(root / "g.gam", rows)
            graph = graph_fixture(root / "graph.sqlite", [(n, s, 5) for n, s in sequences.items()])
            (root / "nodes.txt").write_text("1\n2\n3\n")
            with redirect_stdout(io.StringIO()):
                build(build_args(gam=str(gam), graph_index=str(graph), nodes=str(root / "nodes.txt"), output=str(root / "o"),
                                 min_variants=1, width=11, max_indel_len=5, rows=12, debug_rows=True))
            (site,) = [json.loads(l) for l in (root / "o/variant_summary.ndjson").read_text().splitlines()]
            self.assertEqual((site["candidate_id"], site["alt_count"], site["ref_count"]), ("1:2:DEL:ATGG>@2+3", 5, 4))
            with redirect_stdout(io.StringIO()):
                self.assertTrue(validate(root / "o", str(gam), None, str(graph))["passed"])

    def test_displaced_command_collects_tasks_and_drops_covered_nodes(self):
        import io
        import json
        from contextlib import redirect_stdout
        from pathlib import Path
        import tempfile
        from indexed_gam_pipeline_v2 import orchestrate
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "targets.txt").write_text("5\n9\n")
            (root / "earlier.txt").write_text("7\n")
            config = dict(tasks=2, tensors=str(root / "tensors"), variant_outputs=dict(SNV=.06, INDEL=.08),
                          inputs=dict(nodes=dict(path=str(root / "targets.txt"))))
            (root / "config.json").write_text(json.dumps(config))
            for i, table in enumerate(("4\t2\n7\t1\n", "4\t1\n9\t3\n8\t1\n")):
                folder = root / "tensors/shared" / f"task_{i:04d}"
                folder.mkdir(parents=True)
                (folder / "displaced_nodes.tsv").write_text(table)
            with redirect_stdout(io.StringIO()):
                summary = orchestrate.displaced_nodes(root, root / "supplement.txt", [root / "earlier.txt"], 2)
            self.assertEqual((root / "supplement.txt").read_text(), "4\n")  # 7 excluded, 9 a target, 8 < 2 records
            self.assertEqual((summary["nodes"], summary["records"], summary["already_covered"], summary["below_min_records"]),
                             (1, 3, 1, 2))


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

"""Discovery: counts equal a plain-Python reference of the rule (records decoded and left-normalized
with candidates.decode_alignment), do not depend on the process count, node_stats.json keeps its
exact format, and segments tile the GAM."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import random
import tempfile
import unittest

import numpy as np
import pysam

from .fixtures import af_gam, build_index, encode_varint, graph_fixture, spec_alignment
from .. import discovery, native
from ..candidates import decode_alignment
from ..common import write_json
from ..gam_reader import group

MODULE, INFO = native.load()
needs_module = unittest.skipIf(MODULE is None, f"native module unavailable: {INFO.get('reason')}")


def write_grouped_gam(path, alignments, per_group):
    """Several records per GAM group (as vg writes them), one group per BGZF flush, plus its GAI."""
    with pysam.BGZFile(str(path), "wb") as stream:
        for start in range(0, len(alignments), per_group):
            raws = [a.SerializeToString() for a in alignments[start:start + per_group]]
            stream.write(encode_varint(len(raws) + 1) + encode_varint(3) + b"GAM"
                         + b"".join(encode_varint(len(r)) + r for r in raws))
            stream.flush()
    build_index(path, str(path) + ".gai")
    return path


def synthetic_world(seed, cases):
    """Records of native.synthetic_case with node IDs moved apart (one shared graph), valid records only."""
    rng = random.Random(seed)
    sequences, records = {}, []
    for index in range(cases):
        a, seqs, _, _, _ = native.synthetic_case(rng, index)
        if index % 97 == 0:
            continue  # the generator's malformed-quality records
        shift = 1000 * (index + 1)
        for m in a.path.mapping:
            m.position.node_id += shift
        sequences.update({n + shift: s for n, s in seqs.items()})
        records.append(a)
    return sequences, records


def reference_counts(records, sequences, min_mapq, max_indel, normalized):
    """The rules in plain Python (normalized: candidates.decode_alignment with left-normalization)."""
    stats, fallbacks = {}, 0

    def count(node, imperfect, length):
        s = stats.setdefault(node, [0, 0, 0])
        s[int(imperfect)] += 1
        s[2] = max(s[2], length)

    for a in records:
        if a.mapping_quality <= min_mapq:
            continue
        if normalized:
            try:
                read, _ = decode_alignment(a, sequences, max_indel)
            except Exception:  # noqa: BLE001 -- counted from the edits as written, like the native code
                fallbacks += 1
            else:
                edited = [False] * len(read.visits)
                for c in read.columns:
                    if c.op != "M":
                        edited[c.visit] = True
                for v in read.visits:
                    if v.node:
                        count(v.node, edited[v.index], len(a.sequence))
                continue
        for m in a.path.mapping:
            if m.position.node_id:
                count(m.position.node_id, any(e.from_length != e.to_length or e.sequence for e in m.edit), len(a.sequence))
    nodes = list(stats)
    return ([np.array(nodes, dtype=np.int64)] + [np.array([stats[n][k] for n in nodes], dtype=np.int64) for k in range(3)],
            fallbacks)


def assert_counts(test, actual, expected):
    for a, b in zip(actual[:4], expected):
        np.testing.assert_array_equal(a, b)


class NodeStatsFormatTest(unittest.TestCase):
    def test_streamed_writer_equals_write_json(self):
        rng = random.Random(3)
        for n in (0, 1, 7, 2500):
            nodes = np.array(rng.sample(range(1, 10**9), n), dtype=np.int64)
            cols = [np.array([rng.randint(0, 10**6) for _ in range(n)], dtype=np.int64) for _ in range(3)]
            with tempfile.TemporaryDirectory() as tmp:
                a, b = Path(tmp) / "a.json", Path(tmp) / "b.json"
                old_chunk = discovery.CHUNK
                discovery.CHUNK = 1000  # exercise the chunk joints
                try:
                    discovery.write_node_stats(a, nodes, *cols)
                finally:
                    discovery.CHUNK = old_chunk
                write_json(b, {str(x): dict(perfect=p, not_perfect=q, max_read_length=r)
                               for x, p, q, r in zip(nodes.tolist(), *(c.tolist() for c in cols))})
                self.assertEqual(a.read_bytes(), b.read_bytes())


def world_gam(tmp, seed, cases, per_group):
    """A grouped GAM of synthetic records and a graph index holding their node sequences."""
    sequences, records = synthetic_world(seed, cases)
    gam = write_grouped_gam(Path(tmp) / f"world{seed}.gam", records, per_group)
    graph = graph_fixture(Path(tmp) / f"world{seed}.sqlite", [(n, s, 1) for n, s in sequences.items()])
    return gam, graph, sequences, records


@needs_module
class ParallelScanTest(unittest.TestCase):
    def test_counts_equal_the_reference_for_any_process_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            gam, graph, sequences, records = world_gam(tmp, 11, 400, 7)
            expected, fallbacks = reference_counts(records, sequences, 5, 50, True)
            for processes in (1, 2, 3):
                with self.subTest(processes=processes):
                    actual = discovery.native_counts(gam, None, processes, str(graph), 5, 50)
                    assert_counts(self, actual, expected)
                    self.assertEqual(actual[4]["fallbacks"], fallbacks)
                    self.assertEqual(actual[4]["alignments"], len(records))

    def test_segments_tile_the_gam(self):
        with tempfile.TemporaryDirectory() as tmp:
            gam, graph, _, _ = world_gam(tmp, 5, 300, 3)
            with pysam.BGZFile(str(gam), "rb") as stream:
                total = 0
                while (messages := group(stream)) is not None:
                    total += len(messages)
            for count in (1, 2, 5, 40, 1000):
                cuts = discovery.segments(gam, None, count)
                self.assertEqual(cuts[0][0], 0)
                self.assertIsNone(cuts[-1][1])
                self.assertTrue(all(a[1] == b[0] for a, b in zip(cuts, cuts[1:])))
                seen = sum(discovery._scan_segment((str(gam), c, -1, 50, str(graph)))[4]["alignments"] for c in cuts)
                self.assertEqual(seen, total)

    def test_discover_command_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            gam, graph, _, _ = world_gam(tmp, 3, 200, 5)
            outputs = {}
            for processes in (1, 3):
                out = Path(tmp) / f"p{processes}"
                args = type("Args", (), dict(gam=str(gam), output=str(out), min_mapq=5, node_alt=0.05,
                                             graph_index=str(graph), processes=processes))()
                with redirect_stdout(io.StringIO()):
                    discovery.discover(args)
                outputs[processes] = {f: (out / f).read_bytes() for f in ("node_stats.json", "target_nodes.txt")}
                report = json.loads((out / "discovery_report.json").read_text())
                self.assertEqual((report["rule"], report["processes"]), ("normalized", processes))
            self.assertEqual(outputs[1], outputs[3])
            self.assertTrue(outputs[1]["target_nodes.txt"])


@needs_module
class NormalizedRuleTest(unittest.TestCase):
    def test_native_normalized_counts_equal_the_python_decoder(self):
        sequences, records = synthetic_world(7, 600)
        counter = MODULE.Discovery()
        counter.add_normalized([a.SerializeToString() for a in records], sequences, 5, 50)
        actual = counter.result()
        expected, fallbacks = reference_counts(records, sequences, 5, 50, True)
        assert_counts(self, actual, expected)
        self.assertEqual(actual[4]["fallbacks"], fallbacks)
        self.assertGreater(sum(actual[2]), 0)

    def test_normalized_scan_of_a_gam_with_a_graph_index(self):
        sequences, records = synthetic_world(9, 300)
        with tempfile.TemporaryDirectory() as tmp:
            gam = write_grouped_gam(Path(tmp) / "g.gam", records, 5)
            graph = graph_fixture(Path(tmp) / "graph.sqlite", [(n, s, 1) for n, s in sequences.items()])
            expected, fallbacks = reference_counts(records, sequences, 5, 50, True)
            for processes in (1, 3):
                actual = discovery.native_counts(gam, None, processes, str(graph), 5, 50)
                assert_counts(self, actual, expected)
                self.assertEqual(actual[4]["fallbacks"], fallbacks)

    def test_a_homopolymer_deletion_counts_on_the_node_it_normalizes_to(self):
        sequences = {1: "CA", 2: "A", 3: "AG"}  # forward walk 1 -> 2 -> 3 over CA|A|AG
        # vg put the 1-bp deletion of the A run on node 3; left-normalized it belongs to node 1
        a = spec_alignment([(1, 0, False, [(2, 2, "")]), (2, 0, False, [(1, 1, "")]),
                            (3, 0, False, [(1, 0, ""), (1, 1, "")])], sequences)
        normalized = MODULE.Discovery()
        normalized.add_normalized([a.SerializeToString()], sequences, 5, 50)
        as_dict = lambda r: {int(n): (int(p), int(q)) for n, p, q in zip(r[0], r[1], r[2])}  # noqa: E731
        self.assertEqual(as_dict(reference_counts([a], sequences, 5, 50, False)[0]), {1: (1, 0), 2: (1, 0), 3: (0, 1)})
        self.assertEqual(as_dict(normalized.result()), {1: (0, 1), 2: (1, 0), 3: (1, 0)})


class DiscoverArgumentsTest(unittest.TestCase):
    def test_discover_needs_a_graph_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            gam, _ = af_gam(tmp)
            args = type("Args", (), dict(gam=str(gam), output=str(Path(tmp) / "d"), min_mapq=5, node_alt=0.05))()
            with self.assertRaisesRegex(ValueError, "--graph-index"):
                discovery.discover(args)


if __name__ == "__main__":
    unittest.main()

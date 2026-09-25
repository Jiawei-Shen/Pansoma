"""Discovery: the parallel native scan equals the sequential Python scan (raw rule), node_stats.json
keeps its exact format, segments tile the GAM, and the normalized rule equals the Python decoder."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import random
import tempfile
import unittest

import numpy as np
import pysam

from .fixtures import af_gam, build_index, encode_varint, graph_fixture, spec_alignment, tiny_gam
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
            except Exception:  # noqa: BLE001 -- counted with the raw rule, like the native code
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


@needs_module
class ParallelScanTest(unittest.TestCase):
    def test_parallel_native_equals_sequential_python(self):
        sequences, records = synthetic_world(11, 400)
        with tempfile.TemporaryDirectory() as tmp:
            gams = [tiny_gam(tmp)[0], af_gam(tmp)[0],
                    write_grouped_gam(Path(tmp) / "grouped.gam", records, 7)]
            for gam in gams:
                expected = discovery.python_counts(gam, 5)
                for processes in (1, 2, 3):
                    with self.subTest(gam=gam.name, processes=processes):
                        actual = discovery.native_counts(gam, None, processes, min_mapq=5)
                        assert_counts(self, actual, expected[:4])
                        self.assertEqual(actual[4]["alignments"], expected[4]["alignments"])
                        self.assertEqual(actual[4]["used"], expected[4]["used"])

    def test_segments_tile_the_gam(self):
        _, records = synthetic_world(5, 300)
        with tempfile.TemporaryDirectory() as tmp:
            gam = write_grouped_gam(Path(tmp) / "g.gam", records, 3)
            with pysam.BGZFile(str(gam), "rb") as stream:
                total = 0
                while (messages := group(stream)) is not None:
                    total += len(messages)
            for count in (1, 2, 5, 40, 1000):
                cuts = discovery.segments(gam, None, count)
                self.assertEqual(cuts[0][0], 0)
                self.assertIsNone(cuts[-1][1])
                self.assertTrue(all(a[1] == b[0] for a, b in zip(cuts, cuts[1:])))
                seen = sum(discovery._scan_segment((str(gam), c, False, -1, 50, None))[4]["alignments"] for c in cuts)
                self.assertEqual(seen, total)

    def test_discover_command_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            gam, _ = af_gam(tmp)
            outputs = {}
            for name, extra in (("python", dict(max_alignments=10**9)), ("native", dict(processes=2))):
                fields = dict(gam=str(gam), output=str(Path(tmp) / name), min_mapq=5, node_alt=0.05,
                              max_alignments=None, processes=None, raw=True)
                args = type("Args", (), dict(fields, **extra))()
                with redirect_stdout(io.StringIO()):
                    discovery.discover(args)
                outputs[name] = {f: (Path(tmp) / name / f).read_bytes() for f in ("node_stats.json", "target_nodes.txt")}
                report = json.loads((Path(tmp) / name / "discovery_report.json").read_text())
                self.assertEqual(report["engine"], name)
            self.assertEqual(outputs["python"], outputs["native"])


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
                actual = discovery.native_counts(gam, None, processes, True, 5, 50, str(graph))
                assert_counts(self, actual, expected)
                self.assertEqual(actual[4]["fallbacks"], fallbacks)

    def test_a_homopolymer_deletion_counts_on_the_node_it_normalizes_to(self):
        sequences = {1: "CA", 2: "A", 3: "AG"}  # forward walk 1 -> 2 -> 3 over CA|A|AG
        # vg put the 1-bp deletion of the A run on node 3; left-normalized it belongs to node 1
        a = spec_alignment([(1, 0, False, [(2, 2, "")]), (2, 0, False, [(1, 1, "")]),
                            (3, 0, False, [(1, 0, ""), (1, 1, "")])], sequences)
        raw, normalized = MODULE.Discovery(), MODULE.Discovery()
        raw.add_raw([a.SerializeToString()], 5)
        normalized.add_normalized([a.SerializeToString()], sequences, 5, 50)
        as_dict = lambda r: {int(n): (int(p), int(q)) for n, p, q in zip(r[0], r[1], r[2])}  # noqa: E731
        self.assertEqual(as_dict(raw.result()), {1: (1, 0), 2: (1, 0), 3: (0, 1)})
        self.assertEqual(as_dict(normalized.result()), {1: (0, 1), 2: (1, 0), 3: (1, 0)})


class NormalizedDefaultTest(unittest.TestCase):
    def test_discover_defaults_to_normalized_and_needs_a_graph_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            gam, _ = af_gam(tmp)
            args = type("Args", (), dict(gam=str(gam), output=str(Path(tmp) / "d"), min_mapq=5, node_alt=0.05))()
            with self.assertRaisesRegex(ValueError, "--graph-index"):
                discovery.discover(args)
            if MODULE is not None:
                graph = graph_fixture(Path(tmp) / "graph.sqlite", [(n, "AAAAAA", 1) for n in (10, 20, 30, 40, 50)])
                args.graph_index, args.processes = str(graph), 2
                with redirect_stdout(io.StringIO()):
                    discovery.discover(args)
                report = json.loads((Path(tmp) / "d" / "discovery_report.json").read_text())
                self.assertEqual((report["rule"], report["engine"]), ("normalized", "native"))


if __name__ == "__main__":
    unittest.main()

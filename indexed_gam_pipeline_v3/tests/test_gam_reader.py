"""GAI v1 reading (and the fixture GAI writer), indexed retrieval and the CLI helpers around them."""
from collections import Counter
from contextlib import redirect_stdout
import gzip
import io
import json
from pathlib import Path
import struct
import tempfile
import unittest

from .fixtures import build_args, build_index, encode_varint, simple_alignment, tiny_gam, write_gam
from ..common import batches, load_nodes
from ..gam_reader import IndexedGam, scan_gam, varint
from ..run import discover


def multiset(alignments):
    return Counter(a.SerializeToString() for a in alignments)


class GamReaderTest(unittest.TestCase):
    def test_indexed_fetch_matches_full_scan_for_every_node_set(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = tiny_gam(directory)
            reader = IndexedGam(path)
            self.assertEqual(reader.version, 1)
            for nodes in ({10}, {30}, {1000}, {10, 20, 30}, {999999}, {10, 20, 30, 1000}):
                expected = multiset(a for a in scan_gam(path)
                                    if any(m.position.node_id in nodes for m in a.path.mapping))
                self.assertEqual(multiset(reader.fetch(nodes)), expected, nodes)

    def test_group_cache_preserves_order_duplicates_and_is_bounded(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = tiny_gam(directory)
            for limit in (1, 2048, 64 * 1024 * 1024):  # 1 byte: no group fits, every fetch decodes afresh
                reader = IndexedGam(path, cache_bytes=limit)
                for nodes in ({10, 30}, {20, 1000}, {10}, {10, 20, 30, 1000}, {10, 30}):
                    # file order, each record once
                    expected = [a.SerializeToString() for a in scan_gam(path)
                                if any(m.position.node_id in nodes for m in a.path.mapping)]
                    actual = list(reader.fetch(nodes))
                    self.assertEqual([a.SerializeToString() for a in actual], expected)
                    if actual:  # callers must not be able to mutate a cached record
                        actual[0].path.mapping[0].position.node_id = 99999
                    self.assertEqual([a.SerializeToString() for a in reader.fetch(nodes)], expected)
                self.assertLessEqual(reader.cache_stats["peak_accounted_bytes"], limit)
                if limit == 64 * 1024 * 1024:
                    self.assertGreater(reader.cache_stats["group_hits"], 0)
                    self.assertEqual(reader.cache_stats["indexed_records"], 6)
            for limit in (0, -1):  # the uncached mode is gone; --gam-cache-mb must be >= 1
                with self.assertRaisesRegex(ValueError, "positive"):
                    IndexedGam(path, cache_bytes=limit)
            reader = IndexedGam(path)
            with path.open("ab") as stream:
                stream.write(b"changed")
            with self.assertRaisesRegex(ValueError, "GAM changed"):
                list(reader.fetch({10}))

    def test_build_index_reports_counts_and_reproduces_retrieval(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = tiny_gam(directory)
            index = Path(directory) / "rebuilt.gai"
            report = build_index(path, index)
            self.assertEqual((report["alignments"], report["groups"]), (6, 6))
            nodes = {10, 30, 1000}
            expected = multiset(a for a in scan_gam(path) if any(m.position.node_id in nodes for m in a.path.mapping))
            self.assertEqual(multiset(IndexedGam(path, index).fetch(nodes)), expected)

    def test_overlapping_and_duplicate_index_runs_are_merged(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = tiny_gam(directory)
            reader = IndexedGam(path)
            original = multiset(reader.fetch({10, 20, 30}))
            reader.bins += list(reader.bins)
            self.assertEqual(multiset(reader.fetch({10, 20, 30})), original)

    def test_equivalent_bgzf_eof_offset_is_accepted_but_mid_group_end_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = tiny_gam(directory)
            reader = IndexedGam(path)
            expected = multiset(reader.fetch({30}))
            data = path.read_bytes()
            blocks, offset = [], 0
            while offset < len(data):
                size = struct.unpack_from("<H", data, offset + 16)[0] + 1
                blocks.append((offset, struct.unpack_from("<I", data, offset + size - 4)[0]))
                offset += size
            (eof_block, eof_size), (data_block, data_size) = blocks[-1], blocks[-2]
            self.assertEqual(eof_size, 0)
            canonical, alternate = eof_block << 16, (data_block << 16) | data_size
            changed = False
            for _, _, runs in reader.bins:
                for i, (start, end) in enumerate(runs):
                    if end == canonical:
                        runs[i], changed = (start, alternate), True
            self.assertTrue(changed)
            self.assertEqual(multiset(reader.fetch({30})), expected)
            reader = IndexedGam(path)
            start, end = reader.bins[0][2][0]
            reader.bins[0][2][0] = (start, start + 1)
            with self.assertRaisesRegex(ValueError, "group boundary"):
                list(reader.fetch({10, 20, 30, 1000}))

    def test_foreign_tagged_groups_are_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            rows = [simple_alignment((10,)), simple_alignment((10, 30)), simple_alignment((30,)), simple_alignment((1000,))]
            path = write_gam(Path(directory) / "tagged.gam", rows, foreign_groups=(0, 2))
            self.assertEqual(multiset(scan_gam(path)), multiset(rows))
            reader = IndexedGam(path)
            for nodes in ({10}, {30}, {1000}, {10, 30, 1000}):
                expected = multiset(a for a in rows if any(m.position.node_id in nodes for m in a.path.mapping))
                self.assertEqual(multiset(reader.fetch(nodes)), expected, nodes)
            self.assertEqual(multiset(IndexedGam(path, cache_bytes=1).fetch({10, 30})), multiset(rows[:3]))

    def test_bad_index_and_truncated_varint_fail(self):
        with self.assertRaises(ValueError):
            varint(io.BytesIO(b"\x80"))
        self.assertEqual(varint(io.BytesIO(encode_varint(300))), 300)
        with tempfile.TemporaryDirectory() as directory:
            path, _ = tiny_gam(directory)
            with gzip.open(str(path) + ".gai", "wb") as stream:
                stream.write(b"GAI!" + encode_varint(99))
            with self.assertRaisesRegex(ValueError, "Unsupported"):
                IndexedGam(path)
            with gzip.open(str(path) + ".gai", "wb") as stream:  # GAI v0: no magic, no version (0 bins, 0 windows)
                stream.write(encode_varint(0) + encode_varint(0))
            with self.assertRaisesRegex(ValueError, "Unsupported GAI"):
                IndexedGam(path)

    def test_helpers(self):
        self.assertEqual(list(batches([10, 20, 30, 500], 2, 100)), [[10, 20], [30], [500]])
        with tempfile.TemporaryDirectory() as directory:
            nodes = Path(directory) / "nodes.txt"
            nodes.write_text("30\n10\n10\n\n20\n")
            self.assertEqual(load_nodes(nodes), [10, 20, 30])
            nodes.write_text("0\n")
            with self.assertRaises(ValueError):
                load_nodes(nodes)

    def test_discover_selects_imperfect_nodes(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = tiny_gam(directory)
            output = Path(directory) / "discovery"
            args = build_args(gam=str(path), output=str(output), min_mapq=5, node_alt=0.05, max_alignments=None, raw=True)
            with redirect_stdout(io.StringIO()):
                discover(args)
            stats = json.loads((output / "node_stats.json").read_text())
            self.assertEqual(stats["10"], dict(perfect=0, not_perfect=3, max_read_length=12))
            self.assertEqual((output / "target_nodes.txt").read_text().split(), ["10", "20", "30", "1000"])
            report = json.loads((output / "discovery_report.json").read_text())
            self.assertEqual((report["alignments_scanned"], report["nodes_selected"], report["exploratory"]),
                             (6, 4, False))
            self.assertNotIn("max_nodes", report)


if __name__ == "__main__":
    unittest.main()

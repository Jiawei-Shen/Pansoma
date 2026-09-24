"""Synthetic regression checks; no external graph or vg executable needed."""
from collections import Counter
import argparse
from contextlib import redirect_stdout
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import struct
import sqlite3
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
import numpy as np
import pysam
from indexed_gam_pipeline import vg_pb2
from indexed_gam_pipeline.gam_reader import IndexedGam, build_index, scan_gam, varint
from indexed_gam_pipeline.segments import on_chromosome, raw_segments, orient
from indexed_gam_pipeline.run import batches, build, discover, node_records
from pangenome_ml_data_generation.tensors.builders import build_tensors_from_segments


def vi(value):
    out = bytearray()
    while value > 127:
        out.append((value & 127) | 128)
        value >>= 7
    out.append(value)
    return bytes(out)


def alignment(nodes=(10,), reverse=False, mapq=60, mutation=True):
    a = vg_pb2.Alignment(name="same-name", mapping_quality=mapq)
    for nid in nodes:
        m = a.path.mapping.add()
        m.position.node_id = nid
        m.position.is_reverse = reverse
        if mutation:
            m.edit.add(from_length=2, to_length=2)
            m.edit.add(from_length=1, to_length=1, sequence="T")
            m.edit.add(from_length=3, to_length=3)
            a.sequence += "AATAAA"
        else:
            m.edit.add(from_length=6, to_length=6)
            a.sequence += "AAAAAA"
    a.quality = bytes([30] * len(a.sequence))
    return a


def fixture(directory):
    path = Path(directory) / "tiny.gam"
    # Includes distinct records sharing a name, exact duplicates, a mapping
    # whose non-initial node is queried, and a wide cross-bin alignment.
    rows = [alignment((10,)), alignment((10,)), alignment((10, 30)),
            alignment((20,)), alignment((20, 1000)), alignment((30,))]
    bins = {}
    with pysam.BGZFile(str(path), "wb") as stream:
        for a in rows:
            start = stream.tell()
            raw = a.SerializeToString()
            stream.write(vi(2) + vi(3) + b"GAM" + vi(len(raw)) + raw)
            stream.flush()  # Store canonical BGZF offsets at block boundaries.
            end = stream.tell()
            ids = [m.position.node_id for m in a.path.mapping]
            shift = max(1, (min(ids) ^ max(ids)).bit_length())
            number = (min(ids) >> shift) + ((1 << (64 - shift)) - 1)
            bins.setdefault(number, []).append((start, end))
    payload = b"GAI!" + vi(1) + vi(len(bins))
    for number, runs in bins.items():
        payload += vi(number) + vi(len(runs))
        for start, end in runs:
            payload += vi(start) + vi(end)
    payload += vi(0)  # Query implementation does not require window speedup.
    with gzip.open(str(path) + ".gai", "wb") as stream:
        stream.write(payload)
    return path, rows


class PipelineTest(unittest.TestCase):
    def test_group_cache_preserves_order_duplicates_and_is_bounded(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = fixture(directory)
            for limit in (0, 1, 2048, 64 * 1024 * 1024):
                reader = IndexedGam(path, cache_bytes=limit)
                uncached = IndexedGam(path, cache_bytes=0)
                for nodes in ({10,30}, {20,1000}, {10}, {10,20,30,1000}, {10,30}):
                    expected = [a.SerializeToString() for a in uncached.fetch(nodes)]
                    actual = list(reader.fetch(nodes))
                    self.assertEqual([a.SerializeToString() for a in actual], expected)
                    # Callers cannot mutate a cached record or its node index.
                    if actual:
                        actual[0].path.mapping[0].position.node_id = 99999
                    self.assertEqual([a.SerializeToString() for a in reader.fetch(nodes)],expected)
                self.assertLessEqual(reader.cache_stats['peak_accounted_bytes'],limit)
                if limit == 64 * 1024 * 1024:
                    self.assertGreater(reader.cache_stats['group_hits'],0)
                    self.assertEqual(reader.cache_stats['indexed_records'],6)
            with self.assertRaisesRegex(ValueError,'nonnegative'):
                IndexedGam(path,cache_bytes=-1)
            reader = IndexedGam(path)
            with path.open('ab') as stream:
                stream.write(b'changed')
            with self.assertRaisesRegex(ValueError,'GAM changed'):
                list(reader.fetch({10}))

    def test_graph_sqlite_batched_queries_and_missing_nodes(self):
        with tempfile.TemporaryDirectory() as directory:
            db = Path(directory)/'graph.sqlite'
            with sqlite3.connect(db) as c:
                c.execute('CREATE TABLE nodes (node_id INTEGER PRIMARY KEY, seq TEXT)')
                c.executemany('INSERT INTO nodes VALUES (?,?)',[(n,'ac') for n in range(1,1002)])
            args = argparse.Namespace(node_json=None,node_sqlite=str(db),gfa=None)
            records = node_records(args,set(range(1,1002)))
            self.assertEqual(len(records),1001)
            self.assertTrue(all(r['sequence']=='AC' for r in records.values()))
            with self.assertRaisesRegex(ValueError,'Missing graph sequences'):
                node_records(args,{1,1002})

    def test_discovery_preserves_writer_statistics_schema(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = fixture(directory)
            output = Path(directory) / "discovery"
            args = argparse.Namespace(gam=str(path), output=str(output),
                                      chr="", chr_nodes=None, node_alt=.05,
                                      max_alignments=None, max_nodes=None)
            with redirect_stdout(io.StringIO()):
                discover(args)
            stats = json.loads((output / "node_stats.json").read_text())
            self.assertEqual(stats["10"], {"perfect": 0, "not_perfect": 3,
                                           "max_read_length": 6,
                                           "max_cigar_length": 6})
            self.assertIn("10", (output / "target_nodes.txt").read_text().splitlines())

    def test_bundled_hg008_sample_channels_and_rows(self):
        folder = ROOT / "indexed_gam_pipeline" / "samples"
        manifest = json.loads((folder / "manifest.json").read_text())
        with gzip.open(folder / "hg008_3nodes_data.npy.gz", "rb") as stream:
            tensors = np.load(io.BytesIO(stream.read()))
        summary = [json.loads(line) for line in
                   (folder / "hg008_3nodes_summary.ndjson").read_text().splitlines()]
        self.assertEqual(tensors.shape, tuple(manifest["shape"]))
        self.assertEqual(str(tensors.dtype), manifest["dtype"])
        self.assertEqual(len(summary), len(tensors))
        for tensor, meta, expected in zip(tensors, summary, manifest["variants"]):
            self.assertEqual((meta["node_id"], meta["variant_key"]),
                             (expected["node_id"], expected["variant_key"]))
            rows = sorted(tensor[:, row, :].tobytes() for row in range(1, 201))
            self.assertEqual(hashlib.sha256(b"".join(rows)).hexdigest(),
                             expected["canonical_read_rows_sha256"])
            self.assertEqual(hashlib.sha256(tensor[:, 0, :].tobytes()).hexdigest(),
                             expected["reference_row_sha256"])
            # SNV center: reference base, ALT read base, BQ, mismatch flag,
            # MAPQ and X CIGAR code occupy the expected five channels.
            center = 50
            encoding = {"A": 20, "C": 30, "G": 50, "T": 70}
            self.assertEqual(int(tensor[0, 0, center]), encoding[meta["v_ref"]])
            self.assertEqual(int(tensor[0, 1, center]), encoding[meta["v_alt"]])
            self.assertGreater(int(tensor[1, 1, center]), 0)
            self.assertEqual(int(tensor[2, 1, center]), 5)
            self.assertGreater(int(tensor[3, 1, center]), 0)
            self.assertEqual(int(tensor[4, 1, center]), 90)

    def test_rebuilt_index_matches_full_scan(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = fixture(directory)
            index = Path(directory) / "rebuilt.gai"
            report = build_index(path, index)
            self.assertEqual(report["alignments"], 6)
            nodes = {10, 30, 1000}
            expected = Counter(a.SerializeToString() for a in scan_gam(path)
                               if any(m.position.node_id in nodes for m in a.path.mapping))
            actual = Counter(a.SerializeToString() for a in IndexedGam(path, index).fetch(nodes))
            self.assertEqual(actual, expected)

    def test_graph_sqlite_supplies_sequences_and_checks_json(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            db = folder / "graph.sqlite"
            with sqlite3.connect(db) as connection:
                connection.execute("CREATE TABLE nodes (node_id TEXT PRIMARY KEY, seq TEXT)")
                connection.execute("INSERT INTO nodes VALUES ('10', 'AAAAAA')")
            args = argparse.Namespace(node_json=None, node_sqlite=str(db), gfa=None)
            self.assertEqual(node_records(args, {10})[10]["sequence"], "AAAAAA")
            (folder / "nodes.json").write_text('[{"node_id": "10", "sequence": "CCCCCC"}]')
            args.node_json = str(folder / "nodes.json")
            with self.assertRaisesRegex(ValueError, "SQLite/JSON sequence mismatch"):
                node_records(args, {10})

    def test_index_matches_full_scan_with_duplicates_and_cross_node_reads(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = fixture(directory)
            reader = IndexedGam(path)
            for nodes in ({10}, {30}, {1000}, {10, 20, 30}, {999999}):
                expected = Counter(a.SerializeToString() for a in scan_gam(path)
                                   if any(m.position.node_id in nodes for m in a.path.mapping))
                got = Counter(a.SerializeToString() for a in reader.fetch(nodes))
                self.assertEqual(expected, got)

    def test_overlapping_index_runs_are_merged(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = fixture(directory)
            reader = IndexedGam(path)
            original = list(reader.fetch({10,20,30}))
            # Redundant and overlapping runs in separate bins must not multiply records.
            reader.bins += list(reader.bins)
            self.assertEqual(Counter(a.SerializeToString() for a in reader.fetch({10,20,30})),
                             Counter(a.SerializeToString() for a in original))

    def test_equivalent_bgzf_eof_offset_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = fixture(directory)
            reader = IndexedGam(path)
            expected = Counter(a.SerializeToString() for a in reader.fetch({30}))

            # The same logical EOF can be encoded at the end of the final data
            # block or at offset zero in the following BGZF EOF block.
            data = path.read_bytes()
            blocks = []
            offset = 0
            while offset < len(data):
                size = struct.unpack_from("<H", data, offset + 16)[0] + 1
                uncompressed_size = struct.unpack_from("<I", data, offset + size - 4)[0]
                blocks.append((offset, size, uncompressed_size))
                offset += size
            eof_block, _, eof_size = blocks[-1]
            data_block, _, data_size = blocks[-2]
            self.assertEqual(eof_size, 0)
            canonical_eof = eof_block << 16
            alternate_eof = (data_block << 16) | data_size

            changed = False
            for _, _, runs in reader.bins:
                for i, (start, end) in enumerate(runs):
                    if end == canonical_eof:
                        runs[i] = (start, alternate_eof)
                        changed = True
            self.assertTrue(changed)
            self.assertNotEqual(canonical_eof, alternate_eof)
            self.assertEqual(Counter(a.SerializeToString() for a in reader.fetch({30})),
                             expected)

    def test_true_mid_group_index_end_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = fixture(directory)
            reader = IndexedGam(path)
            changed = False
            for _, _, runs in reader.bins:
                for i, (start, end) in enumerate(runs):
                    if start < end:
                        runs[i] = (start, start + 1)
                        changed = True
                        break
                if changed:
                    break
            self.assertTrue(changed)
            with self.assertRaisesRegex(ValueError, "group boundary"):
                list(reader.fetch({10, 20, 30, 1000}))

    def test_bad_index_and_truncated_varint_fail(self):
        with self.assertRaises(ValueError):
            varint(io.BytesIO(b"\x80"))
        with tempfile.TemporaryDirectory() as directory:
            path, _ = fixture(directory)
            with gzip.open(str(path) + ".gai", "wb") as stream:
                stream.write(b"GAI!" + vi(99))
            with self.assertRaisesRegex(ValueError, "Unsupported"):
                IndexedGam(path)

    def test_segment_read_cursor_mapq_zero_quality_and_reverse(self):
        a = alignment((10, 20), reverse=True)
        self.assertFalse(on_chromosome(a, "chr1"))
        a.refpos.add(name="chr1")
        self.assertTrue(on_chromosome(a, "chr1"))
        self.assertFalse(on_chromosome(a, "chr2"))
        a.sequence = "CCCCCCAATAAA"
        a.quality = bytes([30] * 11 + [0])
        nid, segment = list(raw_segments(a, {20}))[0]
        self.assertEqual(nid, 20)
        self.assertEqual(segment["read_sequence"], "AATAAA")
        self.assertEqual(segment["processed_quality_values"][-1], 0)
        forward = orient(segment, 10)
        self.assertEqual(forward["offset_on_node"], 4)
        self.assertEqual(forward["read_sequence"], "TTTATT")
        self.assertEqual(forward["processed_quality_values"][0], 0)
        self.assertEqual(list(raw_segments(alignment(mapq=10), {10})), [])
        self.assertEqual(len(list(raw_segments(alignment(mapq=11), {10}))), 1)
        self.assertEqual(list(raw_segments(alignment(mapq=11), {10}, 12)), [])

    def test_batch_ownership(self):
        result = list(batches([10, 20, 30, 500], 2, 100))
        self.assertEqual(result, [[10, 20], [30], [500]])
        a = alignment((10, 20, 30))
        got = [nid for batch in result for nid, _ in raw_segments(a, set(batch))]
        self.assertEqual(got, [10, 20, 30])

    def test_synthetic_shard_and_summary_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = fixture(directory)
            folder = Path(directory)
            (folder / "nodes.txt").write_text("10\n20\n30\n1000\n")
            (folder / "nodes.json").write_text(json.dumps([
                {"node_id": n, "sequence": "aaaaaa"} for n in (10, 20, 30, 1000)]))
            args = argparse.Namespace(command="build", format="legacy", gam=str(path), index=None,
                nodes=str(folder / "nodes.txt"), node_json=str(folder / "nodes.json"),
                gfa=None, output=str(folder / "out"), batch_nodes=2, max_node_span=100,
                max_batch_segments=100, shard_size=2, min_mapq=10, min_af=.05,
                min_variants=1, min_allele_bq=10., variant_type="all", max_indel_len=50)
            with redirect_stdout(io.StringIO()):
                build(args)
            report = json.loads((folder / "out/run_report.json").read_text())
            self.assertEqual(report["status"], "complete")
            self.assertEqual(report["tensors"], 4)
            self.assertEqual(report["shards"], 2)
            summary = [json.loads(line) for line in
                       (folder / "out/variant_summary.ndjson").read_text().splitlines()]
            for i, meta in enumerate(summary):
                self.assertEqual(meta["shard_index"], i // 2)
                self.assertEqual(meta["index_within_shard"], i % 2)
                tensor = np.load(folder / f"out/shard_{i // 2:05d}_data.npy")
                self.assertEqual(tensor.shape, (2, 5, 201, 100))

    def test_shared_core_matches_unchanged_legacy_dat_reader(self):
        spec = importlib.util.spec_from_file_location("legacy_tensor", ROOT / "scripts/generate_testing_tensors.py")
        legacy = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(legacy)
        raw = []
        for reverse in (False, True):
            for mutation in (False, True):
                raw.extend(segment for _, segment in raw_segments(
                    alignment(reverse=reverse, mutation=mutation), {10}))
        # Also cover insertion and deletion CIGARs in both source paths.
        for ops, seq in (([(2, "M"), (1, "I"), (4, "M")], "AATAAAA"),
                         ([(2, "M"), (1, "D"), (3, "M")], "AAAAA")):
            item = dict(raw[0], cigar_ops=ops, read_sequence=seq,
                        original_cigar_str="".join(f"{n}{op}" for n, op in ops),
                        processed_quality_values=[30] * len(seq))
            raw.append(item)
        raw = raw * 3
        R = max(len(s["read_sequence"]) for s in raw)
        C = max(len(s["original_cigar_str"]) for s in raw)
        layout = legacy.make_record_struct(R, C)
        packed = legacy.BLOCK_HDR_PACK.pack(10, len(raw), 0, R, C)
        for s in raw:
            packed += layout.pack(s["offset_on_node"], s["read_sequence"].encode(),
                                  bytes(s["processed_quality_values"]),
                                  s["original_cigar_str"].encode(),
                                  s["mapping_quality"], s["strand"].encode())
        legacy.worker_dat_file = io.BytesIO(packed)
        legacy.GLOBAL_NODE_SEQS = {10: "AAAAAA"}
        expected = legacy.process_single_node_for_pileup((10, 0, len(raw), .05, 3, 10., "all", 10, 50))
        actual = build_tensors_from_segments(10, "AAAAAA", [orient(s, 6) for s in raw])
        self.assertGreater(len(actual[2]), 0)
        self.assertEqual(actual[3], expected[3])
        np.testing.assert_array_equal(np.stack(actual[2]), np.stack(expected[2]))
        self.assertEqual(actual[2][0].shape, (5, 201, 100))
        self.assertEqual(actual[2][0].dtype, np.int8)


if __name__ == "__main__":
    unittest.main()

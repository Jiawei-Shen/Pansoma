"""The record path (gam_record_index: <gam>.gri or PANSOMA_GAM_RECORD_INDEX=memory) and the GAI walk on BgzfReader
against the reader of commit 0886f70 (gam_reader_oracle.py) on vg-style GAMs: groups of 1-50 records, several
groups per block and groups across blocks, empty blocks, PARAMS_JSON and empty groups, a GAI with vg's bin rule
(contiguous groups of one bin merged into a run), gamsorted and unsorted record orders, node-ID jumps that force
coarse bins, Struct annotations written out of key order (raw bytes != deterministic serialization), identical
records in different groups, MAPQ around the threshold, repeat visits and zero-edit mappings, GAI offsets in
another form of the same position. Error parity: a GAM changed after opening, GAI runs the walk cannot read
(the build refuses them), stale or truncated indexes, corrupt blocks, bad GAIs."""
from contextlib import redirect_stdout
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import random
import struct
import tempfile
import unittest
from unittest.mock import patch
import zlib

import pysam

from . import gam_reader_oracle as oracle
from .fixtures import encode_varint, graph_fixture, tiny_gam
from .test_bgzf_reader import bgzf_bytes, random_cuts
from .. import gam_record_index, orchestrate, vg_pb2
from ..gam_reader import IndexedGam, decode, group
from ..gam_record_index import ENVIRONMENT, Refused

# sha256 of indexed_gam_pipeline_v4/gam_reader.py at 0886f70: the oracle with its import line put back
ORACLE_SHA256 = "080c3e3e999e90b20491fa510e7a4c3699c385b26a259df4536a006d7f45716b"
BLOCKS = ((1, 300), (70_000, 70_200), (5_000_000, 5_000_100))  # node-ID blocks: jumps between them
PROBE = 9_000_000  # a node no record visits


def struct_field(rng):
    """Alignment.annotation (field 100, a Struct) with its map entries in descending key order."""
    entries = []
    for key in sorted(rng.sample(["mapq_uncapped", "fragment", "secondary", "proper", "time", "tag", "z"],
                                 rng.randint(2, 6)), reverse=True):
        kind = rng.random()
        if kind < 0.4:
            value = b"\x11" + struct.pack("<d", rng.random() * 100)
        elif kind < 0.6:
            value = b"\x20" + encode_varint(rng.randint(0, 1))
        else:
            text = "".join(rng.choice("abc") for _ in range(rng.randint(0, 5))).encode()
            value = b"\x1a" + encode_varint(len(text)) + text
        entry = b"\x0a" + encode_varint(len(key)) + key.encode() + b"\x12" + encode_varint(len(value)) + value
        entries.append(b"\x0a" + encode_varint(len(entry)) + entry)
    body = b"".join(entries)
    return encode_varint((100 << 3) | 2) + encode_varint(len(body)) + body


def random_record(rng, index):
    """Raw bytes of a vg-like Alignment: 1-12 node visits in a node-ID block, now and then jumping to another
    block, repeat visits, zero-edit mappings, either strand; 5 % unmapped; an out-of-order annotation."""
    a = vg_pb2.Alignment(name=f"read{index}", mapping_quality=rng.choice([0, 5, 9, 10, 10, 11, 30, 60]))
    if rng.random() > 0.05:
        lo, hi = rng.choice(BLOCKS)
        node = rng.randint(lo, hi)
        reverse = rng.random() < 0.3
        for _ in range(rng.randint(1, 12)):
            m = a.path.mapping.add()
            m.position.node_id, m.position.is_reverse = node, reverse
            m.position.offset = rng.choice([0, 0, 3])
            if rng.random() > 0.15:  # else a mapping without edits
                m.edit.add(from_length=4, to_length=4)
                if rng.random() < 0.3:
                    m.edit.add(from_length=1, to_length=1, sequence="T")
            step = rng.random()
            if step < 0.08:
                lo, hi = rng.choice(BLOCKS)
                node = rng.randint(lo, hi)
            elif step > 0.12:  # else visit the node again
                node += -1 if reverse and node > 1 else 1
        a.sequence = "ACGT" * len(a.path.mapping)
        a.quality = bytes(rng.randint(2, 40) for _ in a.sequence)
    raw = a.SerializeToString(deterministic=True)
    return raw + (struct_field(rng) if rng.random() < 0.7 else b"")


def gamsort_key(raw):
    """vg gamsort's order: the smallest (node_id, is_reverse, offset) of the path; unmapped records first."""
    mappings = decode(raw).path.mapping
    return min(((m.position.node_id, m.position.is_reverse, m.position.offset) for m in mappings), default=(0,))


def common_bin(lo, hi):
    shift = max(1, (lo ^ hi).bit_length())
    return (lo >> shift) + (1 << (64 - shift)) - 1


def write_gai(path, bins):
    """A GAI ('GAI!', format number 1) of {bin number: [(start, end), ...]}, no window table."""
    with gzip.open(path, "wb") as stream:
        stream.write(b"GAI!" + encode_varint(1) + encode_varint(len(bins)))
        for number, runs in sorted(bins.items()):
            stream.write(encode_varint(number) + encode_varint(len(runs)))
            for start, end in runs:
                stream.write(encode_varint(start) + encode_varint(end))
        stream.write(encode_varint(0))


class World:
    """A vg-style GAM with its GAI; groups = [(start tell, end tell, kind, raw messages)] as pysam reads them."""

    def __init__(self, directory, name, seed, records=300, ordered=True, foreign=0.08, index_foreign=False,
                 alternate_starts=0.0, alternate_eof=False, raws=None):
        rng = random.Random(seed)
        self.path = Path(directory) / f"{name}.gam"
        raws = list(raws) if raws is not None else [random_record(rng, i) for i in range(records)]
        for _ in range(len(raws) // 25):  # identical records, in other groups
            raws.insert(rng.randrange(len(raws) + 1), rng.choice(raws))
        if ordered:
            raws.sort(key=gamsort_key)
        else:
            rng.shuffle(raws)
        payload, marks, kinds, at = bytearray(), [], [], 0
        while at < len(raws) or not kinds:
            marks.append(len(payload))
            kind = rng.random()
            if kind < foreign:
                text = b'{"parameters": true}'
                payload += encode_varint(2) + encode_varint(11) + b"PARAMS_JSON" + encode_varint(len(text)) + text
                kinds.append("foreign")
            elif kind < foreign + 0.02:
                payload += encode_varint(0) if rng.random() < 0.5 else encode_varint(1) + encode_varint(3) + b"GAM"
                kinds.append("empty")
            else:
                chunk = raws[at:at + rng.choice([1, 2, 5, rng.randint(1, 50)])]
                at += len(chunk)
                payload += encode_varint(len(chunk) + 1) + encode_varint(3) + b"GAM"
                for raw in chunk:
                    marks.append(len(payload))
                    payload += encode_varint(len(raw)) + raw
                kinds.append("GAM")
        data, self.blocks = bgzf_bytes(random_cuts(rng, bytes(payload), boundaries=marks, empty=0.08))
        self.path.write_bytes(data)
        self.groups = []
        with pysam.BGZFile(str(self.path), "rb") as stream:
            for kind in kinds:
                start = stream.tell()
                messages = group(stream)
                self.groups.append((start, stream.tell(), kind, messages))
            assert group(stream) is None
        self.raws = raws
        self.nodes = sorted({m.position.node_id for raw in raws for m in decode(raw).path.mapping})
        self.ordered, self.rng = ordered, rng
        self.write_index(index_foreign, alternate_starts, alternate_eof)

    def other_form(self, offset):
        """The end-of-previous-block form of a block-start offset, if the previous block holds data."""
        if offset & 0xFFFF:
            return None
        sizes = {address: len(chunk) for address, chunk in self.blocks}
        previous = [address for address, _ in self.blocks if address < offset >> 16]
        if not previous or not sizes[previous[-1]]:
            return None
        return (previous[-1] << 16) | sizes[previous[-1]]

    def write_index(self, index_foreign=False, alternate_starts=0.0, alternate_eof=False, edit=None):
        """The GAI by vg's rule: each group in common_bin(min, max) of its node IDs (unmapped records and indexed
        foreign groups count as node 0), contiguous groups of one bin in one run. `alternate_starts`: that share
        of the run starts in their other form; `alternate_eof`: the run ending at EOF ends at the other form;
        edit(bins): a last change."""
        bins = {}
        for start, end, kind, messages in self.groups:
            if kind == "GAM" and messages:
                ids = []
                for raw in messages:
                    ids.extend([m.position.node_id for m in decode(raw).path.mapping] or [0])
            elif kind == "foreign" and index_foreign:
                ids = [0]
            else:
                continue
            runs = bins.setdefault(common_bin(min(ids), max(ids)), [])
            if runs and runs[-1][1] == start:
                runs[-1] = (runs[-1][0], end)
            else:
                runs.append((start, end))
        eof = self.groups[-1][1]
        for runs in bins.values():
            for i, (start, end) in enumerate(runs):
                other = self.other_form(start)
                if other is not None and self.rng.random() < alternate_starts:
                    start = other
                other = self.other_form(end)
                if alternate_eof and end == eof and other is not None:
                    end = other
                runs[i] = (start, end)
        if edit is not None:
            edit(bins)
        write_gai(str(self.path) + ".gai", bins)
        self.bins = bins
        return bins


def outcome(reader, nodes, cap, method="capped"):
    """What a fetch gives: [(deterministic bytes, plain bytes, left_out)], sampled, returned; or the error."""
    metrics = {}
    try:
        if method == "capped":
            fetched = reader.fetch_capped(nodes, cap, metrics)
            return ([(a.SerializeToString(deterministic=True), a.SerializeToString(), out) for a, out in fetched],
                    metrics["sampled"], metrics["returned_alignments"])
        return [a.SerializeToString(deterministic=True) for a in reader.fetch(nodes, metrics)], \
            metrics["returned_alignments"]
    except Exception as error:  # noqa: BLE001 -- compared, type and message
        context = error.__context__
        return ("error", type(error).__name__, str(error), type(context).__name__ if context else None,
                str(context) if context else None)


def opened(path, mode, **options):
    with patch.dict(os.environ, {ENVIRONMENT: mode}), redirect_stdout(io.StringIO()), \
            patch("sys.stderr", io.StringIO()):
        return IndexedGam(path, **options)


def quiet_build(path, **options):
    with patch("sys.stderr", io.StringIO()):
        return gam_record_index.build(path, force=True, **options)


def batches(rng, nodes):
    """A sequence of node sets: consecutive slices of the targets (as the builder forms them), random sets,
    nodes no record visits."""
    found = []
    while len(found) < rng.randint(2, 6):
        kind = rng.random()
        if kind < 0.5:
            i = rng.randrange(len(nodes))
            found.append(nodes[i:i + rng.choice([1, 3, 10, 40, 200])])
        elif kind < 0.8:
            found.append(rng.sample(nodes, min(len(nodes), rng.randint(1, 30))))
        else:
            found.append([rng.choice([2_000_000, 400, 70_150]), *rng.sample(nodes, 2)])
    return found


class RecordPathTest(unittest.TestCase):
    """~500 random batch sequences: oracle == GAI walk == record path from <gam>.gri == record path in memory."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        directory = Path(cls.tmp.name)
        cls.worlds = []
        settings = [dict(ordered=True), dict(ordered=False), dict(ordered=True, foreign=0.2, index_foreign=True),
                    dict(ordered=True, alternate_starts=0.5, alternate_eof=True), dict(ordered=False, records=80),
                    dict(ordered=True, records=600, alternate_starts=0.3)]
        layouts = [dict(), dict(stride=3, bits=(2, 26), fraction=0.3), dict(stride=1, bits=(2, 26), fraction=0.01),
                   dict(stride=7, bits=(12, 26), fraction=0.5)]
        for i, setting in enumerate(settings):
            for j, layout in enumerate(layouts[:2] if i % 2 else layouts):
                (directory / f"w{i}_{j}").mkdir()
                world = World(directory / f"w{i}_{j}", "world", 100 * i + j, **setting)
                with patch.object(gam_record_index, "CHECKPOINT_STRIDE", layout.get("stride", 4096)), \
                        patch.object(gam_record_index, "W_BITS", layout.get("bits", (12, 26))), \
                        patch.object(gam_record_index, "LONG_FRACTION", layout.get("fraction", 0.01)):
                    world.summary = quiet_build(world.path, processes=1 + (i + j) % 3)
                cls.worlds.append(world)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_oracle_is_the_reader_of_0886f70(self):
        text = Path(oracle.__file__).read_text().replace("from .. import vg_pb2\n", "from . import vg_pb2\n", 1)
        self.assertEqual(hashlib.sha256(text.encode()).hexdigest(), ORACLE_SHA256)

    def test_worlds_cover_the_cases(self):
        raws = [raw for w in self.worlds for raw in w.raws]
        self.assertTrue(any(raw != decode(raw).SerializeToString(deterministic=True) for raw in raws))
        self.assertTrue(any(len(w.raws) != len(set(w.raws)) for w in self.worlds))
        self.assertTrue(any(s["long_records"] for s in (w.summary for w in self.worlds)))
        self.assertTrue(any(not s["sorted_by_lo1"] for s in (w.summary for w in self.worlds)))
        self.assertTrue(any(s["sorted_by_lo1"] for s in (w.summary for w in self.worlds)))
        kinds = {kind for w in self.worlds for _, _, kind, _ in w.groups}
        self.assertEqual(kinds, {"GAM", "foreign", "empty"})
        self.assertTrue(any(len(m) > 20 for w in self.worlds for _, _, kind, m in w.groups if kind == "GAM"))
        self.assertTrue(any(e >> 16 != s >> 16 for w in self.worlds for s, e, _, _ in w.groups))  # across blocks
        self.assertTrue(any(not b for w in self.worlds for _, b in w.blocks[:-1]))  # empty blocks inside
        for w in self.worlds:  # GAI offsets in another form than the walk's tell() (both mapped)
            index = gam_record_index.load(gam_record_index.gri_path(w.path))
            exact = set(index.groups.tolist())
            if any(o not in exact for runs in w.bins.values() for run in runs for o in run):
                break
        else:
            self.fail("no world has GAI offsets in another form")

    def test_random_batch_sequences_match_the_oracle(self):
        rng = random.Random(20261005)
        used = dict(record=0, walk=0)
        for sequence in range(500):
            world = self.worlds[sequence % len(self.worlds)]
            cap, min_mapq = rng.choice([0, 1, 3, 7, 10 ** 6]), rng.choice([None, 10])
            cache = rng.choice([1, 4096, 64 << 20])
            mode = rng.choice(["off", "require", "memory"])
            old = oracle.IndexedGam(world.path, cache_bytes=cache, min_mapq=min_mapq)
            new = opened(world.path, mode, cache_bytes=cache, min_mapq=min_mapq)
            if mode != "off":
                self.assertIsNotNone(new.cache_stats["record_index"])
            for nodes in batches(rng, world.nodes):
                self.assertEqual(outcome(new, nodes, cap), outcome(old, nodes, cap), (sequence, mode, nodes))
                if rng.random() < 0.3:
                    self.assertEqual(outcome(new, nodes, cap, "fetch"), outcome(old, nodes, cap, "fetch"))
            self.assertEqual(new.cache_stats["fallback_fetches"], 0)
            used["record" if mode != "off" else "walk"] += new.cache_stats["record_fetches"] or 1
        self.assertGreater(min(used.values()), 100)

    def test_every_visit_is_found_by_a_one_node_query(self):
        for world in self.worlds[:4]:
            new = opened(world.path, "require")
            visits = {}
            for raw in world.raws:
                for node in {m.position.node_id for m in decode(raw).path.mapping}:
                    visits[node] = visits.get(node, 0) + 1
            for node, count in visits.items():
                self.assertEqual(len(list(new.fetch([node]))), count, node)

    def test_sorted_files_read_record_windows(self):
        """A sorted index keeps its records on disk and reads windows of them, moving forward."""
        world = next(w for w in self.worlds if w.summary["sorted_by_lo1"] and w.summary["W"] < 4096)
        new = opened(world.path, "require")
        self.assertIsNone(new._records._records)
        for i in range(0, len(world.nodes), 25):
            list(new.fetch(world.nodes[i:i + 25]))
            start, rows = new._records._cache
            self.assertLessEqual(len(rows), new._records.header["records"])
        self.assertGreater(new.cache_stats["preads"], 0)


class RefusalAndErrorTest(unittest.TestCase):
    """GAMs the record path must not take: the build refuses them and the GAI walk raises as the oracle does."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.directory = Path(self.tmp.name)
        environment = patch.dict(os.environ)
        environment.start()
        self.addCleanup(environment.stop)
        os.environ.pop(ENVIRONMENT, None)

    def tearDown(self):
        self.tmp.cleanup()

    def assert_parity(self, world, modes=("off", "auto", "memory"), sequences=12, seed=3):
        """Every mode fetches what the oracle fetches, or fails as it does; returns the oracle's errors."""
        rng, errors = random.Random(seed), set()
        for sequence in range(sequences):
            cap, min_mapq = rng.choice([0, 3, 10 ** 6]), rng.choice([None, 10])
            for mode in modes:
                old = oracle.IndexedGam(world.path, cache_bytes=4096, min_mapq=min_mapq)
                new = opened(world.path, mode, cache_bytes=4096, min_mapq=min_mapq)
                for nodes in batches(random.Random(sequence), world.nodes) + [[PROBE], [PROBE, *world.nodes[:9]]]:
                    expected = outcome(old, nodes, cap)
                    self.assertEqual(outcome(new, nodes, cap), expected, (mode, nodes))
                    if expected[0] == "error":
                        errors.add(expected[1:3])
        return errors

    def world(self, name, **options):
        (self.directory / name).mkdir()
        return World(self.directory / name, "w", zlib.crc32(name.encode()), records=150, **options)

    def test_gai_runs_the_walk_cannot_read_are_refused(self):
        def record_start(world):
            """A record start inside a group (not a group start)."""
            for start, end, kind, messages in world.groups:
                if kind == "GAM" and len(messages) > 1:
                    with pysam.BGZFile(str(world.path), "rb") as stream:
                        stream.seek(start)
                        stream.read(1 + 1 + 3)
                        stream.read(len(encode_varint(len(messages[0]))) + len(messages[0]))
                        return stream.tell(), start, end
            raise AssertionError("no multi-record group")

        def other_end(world, bins):
            """A run of one group (not the last) ending at the other form of its end."""
            start, other = next((s, world.other_form(e)) for s, e, _, _ in world.groups[:-1] if world.other_form(e))
            bins.setdefault(common_bin(PROBE, PROBE), []).append((start, other))

        # Each bad run sits alone in the bin of node PROBE (no record visits it), so fetches of PROBE read it alone.
        boundary = {("ValueError", "GAI run does not end on a GAM group boundary")}
        cases = {
            "run end mid-group": (lambda w, b: b.setdefault(common_bin(PROBE, PROBE), []).append(
                (record_start(w)[1], record_start(w)[0])), boundary),
            "run start mid-group": (lambda w, b: b.setdefault(common_bin(PROBE, PROBE), []).append(
                (record_start(w)[0], record_start(w)[2])), None),
            "run end in another form": (other_end, boundary),
        }
        for name, (edit, errors) in cases.items():
            for attempt in range(20):  # a world with a group end at a block start after a block with data
                world = self.world(f"{name.replace(' ', '_')}_{attempt}")
                if any(world.other_form(e) for _, e, _, _ in world.groups[:-1]):
                    break
            world.write_index(edit=lambda bins: edit(world, bins))
            with self.assertRaises(Refused, msg=name):
                quiet_build(world.path)
            self.assertFalse(gam_record_index.gri_path(world.path).exists())
            memory = opened(world.path, "memory")
            self.assertIsNone(memory.cache_stats["record_index"])
            self.assertIn("not indexable", memory.cache_stats["reason"])
            found = self.assert_parity(world, sequences=6)
            self.assertTrue(found, name)
            if errors is not None:
                self.assertEqual(found, errors, name)

    def test_a_gam_changed_after_opening(self):
        world = self.world("changed")
        quiet_build(world.path)
        readers = [oracle.IndexedGam(world.path)] + [opened(world.path, m) for m in ("off", "auto", "memory")]
        with world.path.open("ab") as stream:
            stream.write(b"changed")
        for reader in readers:
            self.assertEqual(outcome(reader, world.nodes[:5], 3)[:3],
                             ("error", "ValueError", "GAM changed after opening its index"))

    def test_stale_or_truncated_indexes_fall_back_or_raise(self):
        world = self.world("stale")
        quiet_build(world.path)
        index = gam_record_index.gri_path(world.path)
        self.assertEqual(opened(world.path, "auto").cache_stats["record_index"], str(index))
        with self.assertRaisesRegex(ValueError, "exists; move it aside"):
            quiet_build(world.path)
        stat = world.path.stat()
        os.utime(world.path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10 ** 9))
        auto = opened(world.path, "auto")
        self.assertIsNone(auto.cache_stats["record_index"])
        self.assertIn("the GAM's size or mtime differs", auto.cache_stats["reason"])
        with self.assertRaisesRegex(ValueError, "require, but .* is stale"):
            opened(world.path, "require")
        self.assert_parity(world, modes=("auto",), sequences=4)
        os.utime(world.path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        self.assertIsNotNone(opened(world.path, "require").cache_stats["record_index"])
        gai = Path(str(world.path) + ".gai")  # the same GAI bytes, rewritten: a new mtime
        gai.write_bytes(gai.read_bytes())
        self.assertIn("the GAI's size, mtime or sha256 differs", opened(world.path, "auto").cache_stats["reason"])
        world.write_index()
        os.unlink(index)
        self.assertEqual([opened(world.path, "auto").cache_stats[k] for k in ("record_index", "reason")], [None, None])
        with self.assertRaisesRegex(ValueError, "require, but there is no"):
            opened(world.path, "require")
        quiet_build(world.path)
        data = index.read_bytes()
        index.write_bytes(data[:-100])
        self.assertIn("truncated", opened(world.path, "auto").cache_stats["reason"])
        with self.assertRaisesRegex(ValueError, "require, but .*truncated"):
            opened(world.path, "require")
        index.write_bytes(b"junk")
        self.assertIn("unreadable", opened(world.path, "auto").cache_stats["reason"])
        index.write_bytes(data)
        with patch.object(gam_record_index, "vg_pb2_sha256", lambda: "0" * 64):
            self.assertIn("vg_pb2.py differs", opened(world.path, "auto").cache_stats["reason"])
        with patch.dict(os.environ, {ENVIRONMENT: "fast"}), self.assertRaisesRegex(ValueError, "must be one of"):
            IndexedGam(world.path)

    def test_corrupt_blocks_are_refused_and_raise_as_before(self):
        world = self.world("corrupt")
        data = bytearray(world.path.read_bytes())
        address = world.blocks[len(world.blocks) // 2][0]
        while not dict(world.blocks)[address]:
            address = next(a for a, _ in world.blocks if a > address)
        size = struct.unpack_from("<H", data, address + 16)[0] + 1
        data[address + size - 8] ^= 1  # CRC32
        world.path.write_bytes(bytes(data))
        world.write_index()
        with self.assertRaisesRegex(Refused, "Error reading from BGZFile"):
            quiet_build(world.path)
        self.assert_parity(world, sequences=6)
        found = {outcome(oracle.IndexedGam(world.path), world.nodes, 0)[:3]}
        self.assertEqual(found, {("error", "OSError", "Error closing BGZFile object")})

    def test_bad_gai_files_fail_as_before(self):
        world = self.world("gai")
        gai = Path(str(world.path) + ".gai")
        good = gai.read_bytes()
        plain = gzip.decompress(good)
        variants = {"truncated gzip": good[:len(good) // 2], "not gzip": b"GAI!" + plain[4:],
                    "gzip crc": good[:-8] + bytes([good[-8] ^ 1]) + good[-7:], "trailing": good + b"x",
                    "magic": gzip.compress(b"GAX" + plain[3:]), "version": gzip.compress(b"GAI!" + encode_varint(2)),
                    "offset past EOF": gzip.compress(b"GAI!" + encode_varint(1) + encode_varint(1) + encode_varint(0)
                                                     + encode_varint(1) + encode_varint(1) + encode_varint(1 << 40)),
                    "trailing data": gzip.compress(plain + b"\x00"), "missing": None}
        for name, content in variants.items():
            if content is None:
                gai.unlink()
            else:
                gai.write_bytes(content)
            found = []
            for make in (lambda: oracle.IndexedGam(world.path), lambda: opened(world.path, "auto")):
                try:
                    make()
                    found.append(None)
                except Exception as error:  # noqa: BLE001 -- compared
                    found.append((type(error).__name__, str(error)))
            self.assertEqual(found[1], found[0], name)
            self.assertIsNotNone(found[0], name)

    def test_short_read_gams_need_force(self):
        path, _ = tiny_gam(self.directory)
        with self.assertRaisesRegex(Refused, "short-read GAM"):
            gam_record_index.build(path)

    def test_check_finds_a_damaged_index(self):
        world = self.world("check")
        quiet_build(world.path)
        report = gam_record_index.check(world.path, sample=50, groups=5)
        self.assertEqual(report["problems"], [])
        self.assertGreater(report["sampled_records"], 0)
        index = gam_record_index.gri_path(world.path)
        header = gam_record_index.info(world.path)
        records = header["sections"]["records"]
        data = bytearray(index.read_bytes())
        data[records["offset"] + 40] ^= 0xFF  # lo2 of the first record
        index.write_bytes(bytes(data))
        problems = gam_record_index.check(world.path, sample=500, groups=50)["problems"]
        self.assertIn("trailer sha256 differs (the file is damaged)", problems)
        self.assertGreater(len(problems), 1)


class PrepareTest(unittest.TestCase):
    def test_prepare_stamps_the_record_index_and_refuses_a_stale_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "w").mkdir()
            world = World(root / "w", "w", 5, records=120)
            graph = graph_fixture(root / "graph.sqlite", [(n, "ACGT", 5) for n in world.nodes])
            nodes = root / "nodes.txt"
            nodes.write_text("".join(f"{n}\n" for n in world.nodes))

            def prepare(name):
                argv = ["prepare", "--root", str(root / name), "--gam", str(world.path), "--nodes", str(nodes),
                        "--graph-index", str(graph), "--haplotypes", "90", "--tasks", "2", "--processes", "1",
                        "--snv-min-af", ".1", "--indel-min-af", ".1", "--merge-shard-size", "0", "--decoder", "python"]
                with redirect_stdout(io.StringIO()):
                    orchestrate.main(argv)
                return json.loads((root / name / "config.json").read_text())

            self.assertNotIn("record_index", prepare("without")["inputs"])
            quiet_build(world.path)
            config = prepare("with")
            self.assertEqual(config["inputs"]["record_index"]["path"], str(gam_record_index.gri_path(world.path)))
            self.assertTrue((root / "with" / "source" / orchestrate.PACKAGE / "gam_record_index.py").exists())
            os.utime(world.path)
            with self.assertRaisesRegex(ValueError, "is stale .*rebuild it .* or move it aside"):
                prepare("stale")
            self.assertFalse((root / "stale").exists())


if __name__ == "__main__":
    unittest.main()

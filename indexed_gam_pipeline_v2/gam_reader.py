"""Read vg GAI (v0/v1) bins and seek BGZF GAM groups without the vg executable.

Format reference: https://github.com/vgteam/vg/blob/master/src/stream_index.cpp
A query scans every bin intersecting the requested node IDs, merges their
virtual-offset runs, decodes the GAM groups in those runs and finally filters
complete Alignment messages by exact node membership.
"""
import bisect
from collections import OrderedDict, defaultdict
import gzip
import io
from pathlib import Path
import sys

import pysam

from indexed_gam_pipeline_v2 import vg_pb2


# --- protobuf stream framing -------------------------------------------------

def varint(stream, allow_eof=False):
    value = 0
    for shift in range(0, 70, 7):
        b = stream.read(1)
        if not b:
            if shift == 0 and allow_eof:
                return None
            raise ValueError("Truncated varint")
        value |= (b[0] & 127) << shift
        if not b[0] & 128:
            if value >= 1 << 64:
                raise ValueError("Varint exceeds uint64")
            return value
    raise ValueError("Invalid varint")


def encode_varint(value):
    out = bytearray()
    while value > 127:
        out.append((value & 127) | 128)
        value >>= 7
    out.append(value)
    return bytes(out)


def exact(stream, size):
    data = stream.read(size)
    if len(data) != size:
        raise ValueError("Truncated GAM/GAI record")
    return data


def group(stream):
    """Read one tagged group as raw messages; None denotes clean EOF.

    Groups with another tag (e.g. giraffe's PARAMS_JSON) are consumed and returned
    empty, so scans and index runs stay aligned with the stream.
    """
    count = varint(stream, allow_eof=True)
    if count is None:
        return None
    if count == 0:
        return []
    tag = exact(stream, varint(stream))
    messages = [exact(stream, varint(stream)) for _ in range(count - 1)]
    return messages if tag == b"GAM" else []


def decode(raw):
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw)
    return alignment


def scan_gam(path, max_alignments=None):
    """Sequential scan of every record; never consults the GAI."""
    seen = 0
    with gzip.open(path, "rb") as stream:
        while True:
            messages = group(stream)
            if messages is None:
                return
            for raw in messages:
                if max_alignments is not None and seen >= max_alignments:
                    return
                yield decode(raw)
                seen += 1


# --- GAI construction ---------------------------------------------------------

def build_index(gam, output):
    """Create a GAI v1 from complete BGZF GAM groups and their node ID spans."""
    bins = {}
    groups = alignments = 0
    with pysam.BGZFile(str(gam), "rb") as stream:
        while True:
            start = stream.tell()
            messages = group(stream)
            if messages is None:
                break
            end = stream.tell()
            node_ids = [m.position.node_id for raw in messages
                        for m in decode(raw).path.mapping if m.position.node_id > 0]
            alignments += len(messages)
            groups += 1
            if not node_ids:
                continue
            lo, hi = min(node_ids), max(node_ids)
            shift = max(1, (lo ^ hi).bit_length())
            if shift > 64:
                raise ValueError("Node ID exceeds uint64")
            number = (lo >> shift) + (1 << (64 - shift)) - 1
            bins.setdefault(number, []).append((start, end))
    with gzip.open(output, "wb") as stream:
        stream.write(b"GAI!" + encode_varint(1) + encode_varint(len(bins)))
        for number, runs in sorted(bins.items()):
            stream.write(encode_varint(number) + encode_varint(len(runs)))
            for start, end in runs:
                stream.write(encode_varint(start) + encode_varint(end))
        stream.write(encode_varint(0))
    return {"groups": groups, "alignments": alignments, "bins": len(bins)}


def equivalent_eof(stream, expected, actual):
    """Both offsets denote the same logical EOF (end of last block vs. start of EOF block)."""
    stream.seek(expected)
    expected_is_eof = not stream.read(1)
    stream.seek(actual)
    actual_is_eof = not stream.read(1)
    stream.seek(actual)
    return expected_is_eof and actual_is_eof


# --- indexed retrieval --------------------------------------------------------

class IndexedGam:
    """Fetch complete alignments touching a node set, with a bounded LRU group cache.

    The cache retains decoded-once groups (raw protobuf bytes plus a node -> record
    posting list), so consecutive node batches that hit the same BGZF groups do not
    re-decode every record. The byte budget covers the retained group data only.
    """

    def __init__(self, gam, index=None, cache_bytes=64 * 1024 * 1024):
        self.gam = Path(gam)
        self.index = Path(index or str(gam) + ".gai")
        if cache_bytes < 0:
            raise ValueError("GAM cache size must be nonnegative")
        self.cache_limit = cache_bytes
        self._groups = OrderedDict()
        self._cache_bytes = 0
        stat = self.gam.stat()
        self._source_stamp = (stat.st_size, stat.st_mtime_ns)
        self.cache_stats = dict(group_hits=0, group_misses=0, indexed_records=0,
                                peak_accounted_bytes=0, limit_bytes=cache_bytes)
        self.bins = []
        self._load_index(stat.st_size)

    def _load_index(self, gam_size):
        with gzip.open(self.index, "rb") as stream:
            data = io.BytesIO(stream.read())
        magic = data.read(4)
        if magic == b"GAI!":
            self.version = varint(data)
        else:
            self.version = 0
            data.seek(0)
        if self.version not in (0, 1):
            raise ValueError(f"Unsupported GAI version: {self.version}")
        for _ in range(varint(data)):
            number = varint(data)
            # Heap-like prefix numbering: (2**n - 1) + n-bit prefix.
            bits = (number + 1).bit_length() - 1
            if bits > 63:
                raise ValueError("Invalid GAI bin")
            prefix = number - ((1 << bits) - 1)
            lo = prefix << (64 - bits)
            hi = ((prefix + 1) << (64 - bits)) - 1
            runs = []
            for _ in range(varint(data)):
                start, end = varint(data), varint(data)
                if start >= end or (end >> 16) > gam_size:
                    raise ValueError("Invalid GAI offsets or mismatched GAM/index")
                runs.append((start, end))
            self.bins.append((lo, hi, runs))
        for _ in range(varint(data)):  # window speedup table, not needed here
            varint(data)
            varint(data)
        if data.read(1):
            raise ValueError("Unexpected trailing GAI data")

    def ranges(self, nodes):
        """Merged virtual-offset runs of every bin intersecting the node set."""
        nodes = sorted(set(nodes))
        if not nodes:
            return []
        runs = []
        for lo, hi, offsets in self.bins:
            i = bisect.bisect_left(nodes, lo)
            if i < len(nodes) and nodes[i] <= hi:
                runs.extend(offsets)
        merged = []
        for start, end in sorted(runs):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        return merged

    def _indexed_group(self, stream, metrics):
        start = stream.tell()
        cached = self._groups.get(start)
        if cached is not None:
            self._groups.move_to_end(start)
            self.cache_stats["group_hits"] += 1
            end, messages, postings, _ = cached
            stream.seek(end)
            return messages, postings
        self.cache_stats["group_misses"] += 1
        messages = group(stream)
        if messages is None:
            raise ValueError("GAI run extends past GAM EOF")
        postings = defaultdict(list)
        for i, raw in enumerate(messages):
            metrics["decoded_alignments"] += 1
            for node in {m.position.node_id for m in decode(raw).path.mapping}:
                postings[node].append(i)
        postings = {node: tuple(indices) for node, indices in postings.items()}
        self.cache_stats["indexed_records"] += len(messages)
        # Containers, raw bytes, keys and posting ints; a fixed allowance for bookkeeping.
        size = (512 + sys.getsizeof(messages) + sum(map(sys.getsizeof, messages))
                + sys.getsizeof(postings)
                + sum(sys.getsizeof(n) + sys.getsizeof(ids) + sum(map(sys.getsizeof, ids))
                      for n, ids in postings.items()))
        if size <= self.cache_limit:
            while self._groups and self._cache_bytes + size > self.cache_limit:
                _, evicted = self._groups.popitem(last=False)
                self._cache_bytes -= evicted[3]
            self._groups[start] = (stream.tell(), messages, postings, size)
            self._cache_bytes += size
            self.cache_stats["peak_accounted_bytes"] = max(
                self.cache_stats["peak_accounted_bytes"], self._cache_bytes)
        return messages, postings

    def fetch(self, nodes, metrics=None):
        """Yield every complete alignment whose path visits any of `nodes`, once per record."""
        stat = self.gam.stat()
        if (stat.st_size, stat.st_mtime_ns) != self._source_stamp:
            raise ValueError("GAM changed after opening its index")
        wanted = set(nodes)
        ranges = self.ranges(wanted)
        if metrics is None:
            metrics = {}
        metrics.update(runs=len(ranges), groups=0, decoded_alignments=0, returned_alignments=0)
        with pysam.BGZFile(str(self.gam), "rb") as stream:
            for start, end in ranges:
                stream.seek(start)
                while stream.tell() < end:
                    metrics["groups"] += 1
                    if self.cache_limit:
                        messages, postings = self._indexed_group(stream, metrics)
                        # Original record order; a record touching several nodes is yielded once.
                        for i in sorted({i for n in wanted for i in postings.get(n, ())}):
                            metrics["decoded_alignments"] += 1
                            metrics["returned_alignments"] += 1
                            yield decode(messages[i])
                        continue
                    messages = group(stream)
                    if messages is None:
                        raise ValueError("GAI run extends past GAM EOF")
                    for raw in messages:
                        metrics["decoded_alignments"] += 1
                        alignment = decode(raw)
                        if any(m.position.node_id in wanted for m in alignment.path.mapping):
                            metrics["returned_alignments"] += 1
                            yield alignment
                actual_end = stream.tell()
                if actual_end != end and not equivalent_eof(stream, end, actual_end):
                    raise ValueError("GAI run does not end on a GAM group boundary")

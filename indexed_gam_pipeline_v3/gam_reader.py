"""Read vg GAI v1 bins and seek BGZF GAM groups without the vg executable.

Format reference: https://github.com/vgteam/vg/blob/master/src/stream_index.cpp
A query scans every bin intersecting the requested node IDs, merges their
virtual-offset runs, decodes the GAM groups in those runs and finally filters
complete Alignment messages by exact node membership.
"""
import bisect
from collections import OrderedDict, defaultdict
import gzip
import hashlib
import heapq
import io
from pathlib import Path
import sys

import numpy as np
import pysam

from . import vg_pb2


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


def equivalent_eof(stream, expected, actual):
    """Both offsets denote the same logical EOF (end of last block vs. start of EOF block)."""
    stream.seek(expected)
    expected_is_eof = not stream.read(1)
    stream.seek(actual)
    actual_is_eof = not stream.read(1)
    stream.seek(actual)
    return expected_is_eof and actual_is_eof


# --- indexed retrieval --------------------------------------------------------

def sample_key(raw):
    """Downsampling rank of a record: its first 8 BLAKE2b bytes (stable across runs and Pythons)."""
    return int.from_bytes(hashlib.blake2b(raw, digest_size=8).digest(), "big")


class IndexedGam:
    """Fetch complete alignments touching a node set, with a bounded LRU group cache.

    The cache retains decoded-once groups (raw protobuf bytes plus a node -> record
    posting list), so consecutive node batches that hit the same BGZF groups do not
    re-decode every record. The byte budget covers the retained group data only.
    With `min_mapq`, records with MAPQ <= min_mapq are left out of the postings (and their
    bytes out of the cache) when a group is first read: fetch never yields or re-parses them.
    """

    def __init__(self, gam, index=None, cache_bytes=64 * 1024 * 1024, min_mapq=None):
        self.gam = Path(gam)
        self.index = Path(index or str(gam) + ".gai")
        if cache_bytes <= 0:
            raise ValueError("GAM cache size must be positive")
        self.cache_limit = cache_bytes
        self.min_mapq = min_mapq
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
        if data.read(4) != b"GAI!":
            raise ValueError("Unsupported GAI: no 'GAI!' magic (GAI v0 is not supported)")
        self.version = varint(data)
        if self.version != 1:
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

    def _bounds(self):
        """Bin bounds as arrays, so a query tests every bin at once (short-read GAIs have millions of
        bins); rebuilt when self.bins is replaced or grown (the arrays keep the list they describe)."""
        described = getattr(self, "_bound_arrays", None)
        if described is None or described[0] is not self.bins or described[1] != len(self.bins):
            self._bound_arrays = (self.bins, len(self.bins), np.array([b[0] for b in self.bins], dtype=np.uint64),
                                  np.array([b[1] for b in self.bins], dtype=np.uint64))
        return self._bound_arrays[2:]

    def ranges(self, nodes):
        """Merged virtual-offset runs of every bin intersecting the node set."""
        nodes = np.unique(np.fromiter(nodes, dtype=np.uint64))
        if not len(nodes):
            return []
        # A bin [lo, hi] intersects the nodes if the first node >= lo is <= hi.
        lo, hi = self._bounds()
        first = np.searchsorted(nodes, lo, side="left")
        inside = first < len(nodes)
        hit = np.zeros(len(self.bins), dtype=bool)
        hit[inside] = nodes[first[inside]] <= hi[inside]
        runs = []
        for b in np.flatnonzero(hit).tolist():
            runs.extend(self.bins[b][2])
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
            alignment = decode(raw)
            if self.min_mapq is not None and alignment.mapping_quality <= self.min_mapq:
                messages[i] = b""  # never yielded; its bytes need not stay cached
                continue
            for node in {m.position.node_id for m in alignment.path.mapping}:
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

    def fetch(self, nodes, metrics=None, sample=None):
        """Yield every complete alignment whose path visits any of `nodes`, once per record, in file order.

        sample=N: yield only the N records with the smallest sample_key (all if there are fewer),
        still in file order; metrics gets sampled_from = the number of records they were drawn from.
        """
        stat = self.gam.stat()
        if (stat.st_size, stat.st_mtime_ns) != self._source_stamp:
            raise ValueError("GAM changed after opening its index")
        wanted = set(nodes)
        ranges = self.ranges(wanted)
        if metrics is None:
            metrics = {}
        metrics.update(runs=len(ranges), groups=0, decoded_alignments=0, returned_alignments=0)
        if sample is not None:
            yield from self._sampled(wanted, ranges, metrics, sample)
            return
        for raw in self._matching(wanted, ranges, metrics):
            metrics["decoded_alignments"] += 1
            metrics["returned_alignments"] += 1
            yield decode(raw)

    def _sampled(self, wanted, ranges, metrics, size):
        """The `size` matching records with the smallest sample_key, in file order."""
        kept, seen = [], 0  # max-heap on the key: (-key, -order, raw)
        for order, raw in enumerate(self._matching(wanted, ranges, metrics)):
            seen += 1
            item = (-sample_key(raw), -order, raw)
            if len(kept) < size:
                heapq.heappush(kept, item)
            elif item > kept[0]:
                heapq.heapreplace(kept, item)
        metrics["sampled_from"] = seen
        for _, _, raw in sorted(kept, key=lambda item: -item[1]):
            metrics["decoded_alignments"] += 1
            metrics["returned_alignments"] += 1
            yield decode(raw)

    def _matching(self, wanted, ranges, metrics):
        """Raw bytes of every record in `ranges` that visits a node of `wanted`, in file order."""
        # Take the cached groups of this fetch first (refreshed, held for the walk), then read the missing
        # ones. A plain LRU walk loses every hit when consecutive fetches cycle through more groups than
        # the cache holds (long reads: neighbouring batches read the same ~50 groups of ~100-270 MiB);
        # this way the groups still cached from the previous fetch are all hits. Records are still
        # yielded in file order, so the fetched records are unchanged.
        pinned = {}
        if self._groups:
            starts = [s for s, _ in ranges]
            for start in list(self._groups):
                k = bisect.bisect_right(starts, start) - 1
                if k >= 0 and start < ranges[k][1]:
                    self._groups.move_to_end(start)
                    pinned[start] = self._groups[start]
                    self.cache_stats["group_hits"] += 1
        with pysam.BGZFile(str(self.gam), "rb") as stream:
            for start, end in ranges:
                stream.seek(start)
                while stream.tell() < end:
                    metrics["groups"] += 1
                    entry = pinned.pop(stream.tell(), None)  # released once used (it may have left the cache)
                    if entry is not None:
                        stream.seek(entry[0])
                        messages, postings = entry[1], entry[2]
                    else:
                        messages, postings = self._indexed_group(stream, metrics)
                    # Original record order; a record touching several nodes is yielded once. The
                    # smaller of (batch nodes, group nodes) is scanned for the intersection.
                    if len(wanted) < len(postings):
                        hits = [n for n in wanted if n in postings]
                    else:
                        hits = [n for n in postings if n in wanted]
                    for i in sorted({i for n in hits for i in postings[n]}):
                        yield messages[i]
                actual_end = stream.tell()
                if actual_end != end and not equivalent_eof(stream, end, actual_end):
                    raise ValueError("GAI run does not end on a GAM group boundary")

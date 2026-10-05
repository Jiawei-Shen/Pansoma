"""Read vg GAI bins ('GAI!' magic, format number 1) and seek BGZF GAM groups without the vg executable.

Format reference: https://github.com/vgteam/vg/blob/master/src/stream_index.cpp
A query scans every bin intersecting the requested node IDs, merges their
virtual-offset runs, decodes the GAM groups in those runs and finally filters
complete Alignment messages by exact node membership.
"""
import bisect
from collections import OrderedDict, defaultdict
import gzip
import hashlib
import io
from pathlib import Path
import sys

import numpy as np
import pysam

from .. import vg_pb2


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


def scan_gam(path):
    """Sequential scan of every record; never consults the GAI."""
    with gzip.open(path, "rb") as stream:
        while True:
            messages = group(stream)
            if messages is None:
                return
            for raw in messages:
                yield decode(raw)


def equivalent_eof(stream, expected, actual):
    """Both offsets denote the same logical EOF (end of last block vs. start of EOF block)."""
    stream.seek(expected)
    expected_is_eof = not stream.read(1)
    stream.seek(actual)
    actual_is_eof = not stream.read(1)
    stream.seek(actual)
    return expected_is_eof and actual_is_eof


# --- indexed retrieval --------------------------------------------------------

def record_key(alignment):
    """Read-cap rank of a record: the first 8 bytes of the builder's record digest (SHA-256 of the
    deterministic serialization, candidates Read.digest), so fetch_capped keeps the records with the
    smallest digest (build.READ_CAP_RULE)."""
    return int.from_bytes(hashlib.sha256(alignment.SerializeToString(deterministic=True)).digest()[:8], "big")


NO_NODES = frozenset()


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
            raise ValueError("Unsupported GAI: no 'GAI!' magic (vg's older GAI format is not supported)")
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

    def _query(self, nodes, metrics):
        stat = self.gam.stat()
        if (stat.st_size, stat.st_mtime_ns) != self._source_stamp:
            raise ValueError("GAM changed after opening its index")
        wanted = set(nodes)
        ranges = self.ranges(wanted)
        metrics.update(runs=len(ranges), groups=0, decoded_alignments=0, returned_alignments=0)
        return wanted, ranges

    def fetch(self, nodes, metrics=None):
        """Yield every complete alignment whose path visits any of `nodes`, once per record, in file order."""
        metrics = {} if metrics is None else metrics
        wanted, ranges = self._query(nodes, metrics)
        for messages, postings, hits in self._group_hits(wanted, ranges, metrics):
            for i in sorted({i for n in hits for i in postings[n]}):  # file order, each record once
                metrics["decoded_alignments"] += 1
                metrics["returned_alignments"] += 1
                yield decode(messages[i])

    def fetch_capped(self, nodes, cap, metrics=None):
        """fetch with the read cap applied while reading: a node visited by more than `cap` records keeps the
        `cap` with the smallest record_key (the builder's record digest; ties in file order), every other node
        all of its records. A node's sample does not depend on the other nodes asked for.

        Returns [(alignment, left_out)] in file order: left_out = the visited nodes of `nodes` whose sample
        leaves the record out (empty for most records); a record that no visited node keeps is not returned.
        metrics gets sampled = {node: records} of the capped nodes. cap=0: every record, nothing left out.
        """
        metrics = {} if metrics is None else metrics
        wanted, ranges = self._query(nodes, metrics)
        raws, groups, counts = [], [], defaultdict(int)
        for messages, postings, hits in self._group_hits(wanted, ranges, metrics):
            local = np.unique(np.fromiter((i for n in hits for i in postings[n]), dtype=np.int64))
            groups.append((len(raws), local, [(n, postings[n]) for n in hits]))
            raws.extend(messages[i] for i in local.tolist())
            for n in hits:
                counts[n] += len(postings[n])
        capped = sorted(n for n, c in counts.items() if cap and c > cap)
        left_out = defaultdict(set)
        visits = None
        if capped:
            members, visits, capped_set = defaultdict(list), np.zeros(len(raws), dtype=np.int64), set(capped)
            for base, local, hits in groups:
                for n, posting in hits:
                    at = base + np.searchsorted(local, np.asarray(posting, dtype=np.int64))
                    visits[at] += 1  # visited nodes of `nodes` per record
                    if n in capped_set:
                        members[n].append(at)
            keys, keyed = np.zeros(len(raws), dtype=np.uint64), np.zeros(len(raws), dtype=bool)
            for n in capped:
                at = np.concatenate(members[n])
                todo = at[~keyed[at]]
                for i in todo.tolist():
                    keys[i] = record_key(decode(raws[i]))
                keyed[todo] = True
                metrics["decoded_alignments"] += len(todo)
                for i in at[np.lexsort((at, keys[at]))[cap:]].tolist():
                    left_out[i].add(n)
        metrics["sampled"] = {n: counts[n] for n in capped}
        fetched = []
        for i, raw in enumerate(raws):
            out = left_out.get(i)
            if out and len(out) == visits[i]:
                continue  # every visited node of the batch leaves it out
            fetched.append((decode(raw), frozenset(out) if out else NO_NODES))
        metrics["decoded_alignments"] += len(fetched)
        metrics["returned_alignments"] += len(fetched)
        return fetched

    def _group_hits(self, wanted, ranges, metrics):
        """(records, postings, hits) of every GAM group in `ranges` that holds a record visiting a node of
        `wanted`, in file order: hits = those nodes."""
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
                    # The smaller of (batch nodes, group nodes) is scanned for the intersection.
                    if len(wanted) < len(postings):
                        hits = [n for n in wanted if n in postings]
                    else:
                        hits = [n for n in postings if n in wanted]
                    if hits:
                        yield messages, postings, hits
                actual_end = stream.tell()
                if actual_end != end and not equivalent_eof(stream, end, actual_end):
                    raise ValueError("GAI run does not end on a GAM group boundary")

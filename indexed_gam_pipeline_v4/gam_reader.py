"""Read vg GAI bins ('GAI!' magic, format number 1) and fetch the GAM records visiting a node set without vg.

Format reference: https://github.com/vgteam/vg/blob/master/src/stream_index.cpp
A query scans every bin intersecting the requested node IDs and merges their virtual-offset runs; it returns
the records of the GAM groups in those runs that visit a requested node (and pass min_mapq), in file order.
The records are found one of two ways, with the same result:

* GAI walk: every group of the runs is read (BgzfReader: htslib's BGZF seek/tell/read over 16 MiB preads) and
  indexed node -> records (`_postings`); a bounded LRU cache keeps the indexed groups.
* record path, when the GAM has a per-record index (<gam>.gri, gam_record_index; built for long reads): only
  the records of those groups whose node intervals hold a requested node and whose MAPQ passes are read
  (coalesced preads of their BGZF blocks) and indexed by the same `_postings`.
"""
import bisect
from collections import OrderedDict, defaultdict
import gzip
import hashlib
import io
import os
from pathlib import Path
import sys
import zlib

import numpy as np

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


# --- BGZF ---------------------------------------------------------------------

MAX_BLOCK = 1 << 16  # BGZF block size limit, compressed and inflated
EXTENT_GAP = 256 << 10  # record path: selected records closer than this share one read
EXTENT_BYTES = 64 << 20  # record path: bytes per read at most (a larger record is read alone)


def _corrupt(reason):
    """pysam's error for a BGZF read htslib fails (the same type and message); `reason` is its cause."""
    error = OSError("Error reading from BGZFile")
    error.__cause__ = ValueError(reason)
    return error


def block_size(header):
    """Total size of the BGZF block whose first 18 bytes are `header`, by htslib's header check (htslib would
    read a gzip member without the BGZF extra field on in gzip mode; such a block is refused here)."""
    if (len(header) != 18 or header[:3] != b"\x1f\x8b\x08" or not header[3] & 4
            or header[10:16] != b"\x06\x00BC\x02\x00"):
        raise _corrupt("invalid or truncated BGZF block header")
    size = int.from_bytes(header[16:18], "little") + 1
    if size < 18:
        raise _corrupt("invalid BGZF block size")
    return size


def inflate(block, size, strict=False):
    """Data of one BGZF block of `size` bytes as htslib inflates it: block[18:] raw-inflated (whatever follows
    the deflate stream ignored) to at most 64 KiB, its CRC32 checked; ISIZE is checked only when `strict`."""
    if len(block) != size:
        raise _corrupt("truncated BGZF block")
    try:
        data = zlib.decompress(memoryview(block)[18:], -15)
    except zlib.error as error:
        raise _corrupt(f"invalid deflate data: {error}") from None
    if len(data) > MAX_BLOCK or zlib.crc32(data) != int.from_bytes(block[size - 8:size - 4], "little"):
        raise _corrupt("BGZF block CRC32 mismatch or over 64 KiB")
    if strict and len(data) != int.from_bytes(block[size - 4:], "little"):
        raise _corrupt("BGZF block ISIZE mismatch")
    return data


def _pread(fd, size, offset):
    """`size` bytes at `offset` (fewer at EOF): os.pread repeated over the short reads BeeGFS returns."""
    parts, got = [], 0
    while got < size:
        part = os.pread(fd, size - got, offset + got)
        if not part:
            break
        parts.append(part)
        got += len(part)
    return parts[0] if len(parts) == 1 else b"".join(parts)


class BgzfReader:
    """The read side of pysam.BGZFile (htslib bgzf_seek / bgzf_tell / bgzf_read) over large os.pread windows.

    htslib reads 32 KiB per request, several times slower on a busy BeeGFS node than a few large reads; this
    reader reads WINDOW bytes at once (bounded by `limit`: the end of the current GAI run plus one block) and
    inflates each block with zlib as htslib does (`inflate`; `strict` also checks ISIZE, which `position`
    relies on). The positions are htslib's: tell() = block address << 16 | offset in the block's data; a read
    that ends a block moves to the next block's start (offset 0) without loading it; a read after a seek to a
    block's end continues in the next block; empty blocks are skipped (a seek's offset is kept across them); at
    EOF read() returns what is left (b"" when nothing) and tell() keeps its value. Corrupt input raises
    OSError("Error reading from BGZFile") as pysam does (the cause names the problem), and close() after it
    OSError("Error closing BGZFile object"), so a `with` block ends with pysam's error too.
    """

    WINDOW = 16 << 20

    def __init__(self, path, strict=False):
        self.strict, self._failed = strict, False
        self.fd = os.open(str(path), os.O_RDONLY)
        self._window, self._window_start, self._read_limit = b"", 0, None
        self._address = self._offset = 0  # htslib block_address, block_offset
        self._data = b""  # the loaded block's data; b"" = none (htslib block_length 0)
        self._next = 0  # file position after the loaded block (htslib htell)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        if self.fd is not None:
            os.close(self.fd)
            self.fd = None
            if self._failed:
                raise OSError("Error closing BGZFile object")

    def limit(self, end):
        """Read ahead no further than the block at virtual offset `end` (None: WINDOW bytes)."""
        self._read_limit = None if end is None else (end >> 16) + MAX_BLOCK

    def _view(self, start, size):
        """Where file byte `start` is in the read-ahead window, read anew from `start` unless the window holds
        [start, start + size) (fewer bytes at EOF)."""
        at = start - self._window_start
        if at < 0 or at + size > len(self._window):
            ahead = self.WINDOW if self._read_limit is None else min(self.WINDOW, self._read_limit - start)
            self._window, self._window_start, at = _pread(self.fd, max(size, ahead), start), start, 0
        return at

    def _bytes(self, start, size):
        """File bytes [start, start + size) (fewer at EOF) through the read-ahead window."""
        at = self._view(start, size)
        return self._window[at:at + size]

    def _load(self):
        """htslib bgzf_read_block: the next non-empty block from the file position, or none at EOF."""
        while True:
            address = self._next
            header = self._bytes(address, 18)
            if not header:
                self._data = b""
                return
            try:
                size = block_size(header)
                data = inflate(self._bytes(address, size), size, self.strict)
            except OSError:
                self._failed = True
                raise
            self._next = address + size
            if data:  # a seek's offset applies to the first block with data
                self._address, self._data = address, data
                return

    def _fill(self):
        """bgzf_read's block step: True when self._data[self._offset] is the next byte, False at EOF."""
        while self._offset >= len(self._data):
            self._load()
            available = len(self._data) - self._offset
            if available > 0:
                return True
            if available < 0:
                self._failed = True
                raise _corrupt("virtual offset past the end of its block")
            if not self._data:
                return False
            self._address, self._offset, self._data = self._next, 0, b""  # a seek to the end of a block
        return True

    def _consume(self, size):
        self._offset += size
        if self._offset == len(self._data):
            self._address, self._offset, self._data = self._next, 0, b""

    def seek(self, offset):
        self._address, self._offset, self._data, self._next = offset >> 16, offset & 0xFFFF, b"", offset >> 16
        return offset

    def tell(self):
        return (self._address << 16) | (self._offset & 0xFFFF)

    def read(self, size):
        parts = []
        while size > 0 and self._fill():
            take = min(size, len(self._data) - self._offset)
            parts.append(self._data[self._offset:self._offset + take])
            self._consume(take)
            size -= take
        return parts[0] if len(parts) == 1 else b"".join(parts)

    def position(self, offset=None, check=False):
        """(block address, offset in its data) of the byte a read from virtual offset `offset` (default: the
        current position) returns first; None at EOF. Offsets that differ only in form (the end of a block and
        the next block's start; an empty block before the data) have the same position. Reads block headers and
        ISIZE only (`check`: also inflates the blocks passed over); tell() and the next read are unchanged."""
        if offset is None:
            if self._offset < len(self._data):
                return self._address, self._offset
            offset = self.tell()
        address, at = offset >> 16, offset & 0xFFFF
        while True:
            header = self._bytes(address, 18)
            if not header:
                if at:
                    raise _corrupt("virtual offset past the end of its block")
                return None
            size = block_size(header)
            start = self._view(address, size)  # the whole block: the next read of it reads nothing more
            tail = self._window[start + size - 4:start + size]
            if len(tail) != 4 or len(self._window) < start + size:
                raise _corrupt("truncated BGZF block")
            isize = int.from_bytes(tail, "little")
            if check and (not isize or at >= isize):
                inflate(self._bytes(address, size), size, self.strict)
            if isize:
                if at < isize:
                    return address, at
                if at > isize:
                    raise _corrupt("virtual offset past the end of its block")
                at = 0
            address += size

    def _varint(self, allow_eof=False):
        """varint(self, allow_eof), decoded from the block data while it lasts."""
        if not self._fill():
            if allow_eof:
                return None
            raise ValueError("Truncated varint")
        data, at = self._data, self._offset
        value = shift = 0
        for i in range(at, min(len(data), at + 10)):
            b = data[i]
            value |= (b & 127) << shift
            if not b & 128:
                if value >= 1 << 64:
                    raise ValueError("Varint exceeds uint64")
                self._consume(i + 1 - at)
                return value
            shift += 7
        if len(data) - at >= 10:
            raise ValueError("Invalid varint")
        return varint(self, allow_eof)  # continues in the next block: byte by byte, nothing consumed yet

    def _exact(self, size):
        """exact(self, size), sliced from the block data when it holds all of it."""
        if 0 < size <= len(self._data) - self._offset:
            data = self._data[self._offset:self._offset + size]
            self._consume(size)
            return data
        if size > 1 << 26:  # pysam allocates the length first: MemoryError for an impossible one (pages untouched)
            try:
                np.empty(size, dtype=np.uint8)
            except MemoryError:
                raise MemoryError from None
        data = self.read(size)
        if len(data) != size:
            raise ValueError("Truncated GAM/GAI record")
        return data

    def read_group(self, positions=None):
        """group(self) without a Python call per byte: the same values (None at a clean EOF, [] for an empty or
        foreign-tagged group) and the same ValueErrors. With `positions` (a list), the (start, end) virtual
        offsets of each returned record (its length varint .. its last byte) are appended."""
        count = self._varint(True)
        if count is None:
            return None
        if count == 0:
            return []
        tag = self._exact(self._varint())
        if positions is None:
            messages = [self._exact(self._varint()) for _ in range(count - 1)]
        else:
            mark, messages = len(positions), []
            for _ in range(count - 1):
                start = self.tell()
                messages.append(self._exact(self._varint()))
                positions.append((start, self.tell()))
            if tag != b"GAM":
                del positions[mark:]
        return messages if tag == b"GAM" else []


def _varint_at(data, at):
    """(value, offset after it) of the varint at data[at], or (None, None) when there is none."""
    value = shift = 0
    for i in range(at, min(len(data), at + 10)):
        value |= (data[i] & 127) << shift
        if not data[i] & 128:
            return (value, i + 1) if value < 1 << 64 else (None, None)
        shift += 7
    return None, None


def _record(buffer, base, start, end, blocks):
    """Message bytes of the record at virtual offsets [start, end) (length varint .. last byte) in `buffer`
    (file bytes from `base`); `blocks` caches inflated blocks {address: (data, size)} while records move on."""
    address, stop, stop_at = start >> 16, end >> 16, end & 0xFFFF
    for old in [a for a in blocks if a < address]:
        del blocks[old]
    parts = []
    while address != stop or stop_at:
        cached = blocks.get(address)
        if cached is None:
            at = address - base
            header = buffer[at:at + 18]
            if len(header) != 18:
                raise ValueError("record index does not match GAM")
            size = block_size(header)
            cached = blocks[address] = (inflate(buffer[at:at + size], size), size)
        data, size = cached
        if address == stop:
            parts.append(data[:stop_at])
            break
        parts.append(data)
        address += size
        if address > stop:
            raise ValueError("record index does not match GAM")
    flat = parts[0] if len(parts) == 1 else b"".join(parts)
    length, at = _varint_at(flat, start & 0xFFFF)
    if length is None or at + length != len(flat):
        raise ValueError("record index does not match GAM")
    return flat[at:]


def read_records(path, starts, ends, stats):
    """Message bytes of the records at virtual offsets [starts[i], ends[i]) (file order; ends as tell() gives
    them): records closer than EXTENT_GAP share one pread of their BGZF blocks (at most EXTENT_BYTES unless one
    record is larger), each block checked and inflated once. ValueError when a record does not end at its end
    offset; stats gets preads and record_bytes_read."""
    starts, ends, messages = starts.tolist(), ends.tolist(), []
    if not starts:
        return messages
    fd = os.open(str(path), os.O_RDONLY)
    try:
        size = os.fstat(fd).st_size

        def block_end(end):
            return end >> 16 if not end & 0xFFFF else min(size, (end >> 16) + MAX_BLOCK)

        i = 0
        while i < len(starts):
            first, last, j = starts[i] >> 16, block_end(ends[i]), i + 1
            while j < len(starts) and (starts[j] >> 16) <= last + EXTENT_GAP and \
                    max(last, block_end(ends[j])) - first <= EXTENT_BYTES:
                last, j = max(last, block_end(ends[j])), j + 1
            buffer = _pread(fd, last - first, first)
            stats["preads"] += 1
            stats["record_bytes_read"] += len(buffer)
            blocks = {}
            messages.extend(_record(buffer, first, starts[k], ends[k], blocks) for k in range(i, j))
            i = j
    finally:
        os.close(fd)
    return messages


# --- GAI ----------------------------------------------------------------------

def read_gai(path, gam_size):
    """(format number, bins [(lo, hi, [(start, end), ...])], sha256 of the file) of a GAI with the 'GAI!' magic,
    format number 1, whose offsets lie in a GAM of `gam_size` bytes. The file is read once (its sha256 stamps a
    record index) and inflated as gzip.open would."""
    raw = Path(path).read_bytes()
    with gzip.GzipFile(fileobj=io.BytesIO(raw), mode="rb") as stream:
        data = io.BytesIO(stream.read())
    if data.read(4) != b"GAI!":
        raise ValueError("Unsupported GAI: no 'GAI!' magic (vg's older GAI format is not supported)")
    version = varint(data)
    if version != 1:
        raise ValueError(f"Unsupported GAI version: {version}")
    bins = []
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
        bins.append((lo, hi, runs))
    for _ in range(varint(data)):  # window speedup table, not needed here
        varint(data)
        varint(data)
    if data.read(1):
        raise ValueError("Unexpected trailing GAI data")
    return version, bins, hashlib.sha256(raw).hexdigest()


# --- indexed retrieval --------------------------------------------------------

def record_key(alignment):
    """Read-cap rank of a record: the first 8 bytes of the builder's record digest (SHA-256 of the
    deterministic serialization, candidates Read.digest), so fetch_capped keeps the records with the
    smallest digest (build.READ_CAP_RULE)."""
    return int.from_bytes(hashlib.sha256(alignment.SerializeToString(deterministic=True)).digest()[:8], "big")


NO_NODES = frozenset()


class IndexedGam:
    """Fetch complete alignments touching a node set: the GAI walk with a bounded LRU group cache, or the
    record path when the GAM has a usable record index (gam_record_index.open_for; PANSOMA_GAM_RECORD_INDEX,
    read when the reader is made and applied at its first fetch: readers that never fetch, such as discovery
    and gam_prep check, neither load nor require an index).

    The cache retains decoded-once groups (raw protobuf bytes plus a node -> record
    posting list), so consecutive node batches that hit the same BGZF groups do not
    re-decode every record. The byte budget covers the retained group data only.
    With `min_mapq`, records with MAPQ <= min_mapq are left out of the postings (and their
    bytes out of the cache) when a group is first read: fetch never yields or re-parses them.
    After the first fetch cache_stats also holds the record index state (mode, record_index, reason) and its
    counters.
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
        from .gam_record_index import environment_mode
        self._records, self._record_mode, self._record_state = None, environment_mode(), None

    def _load_index(self, gam_size):
        self.version, self.bins, sha256 = read_gai(self.index, gam_size)
        stat = self.index.stat()
        self._index_stamp = (stat.st_size, stat.st_mtime_ns, sha256)

    def open_record_index(self):
        """Open the record index by the reader's mode (gam_record_index.open_for), once; the first fetch does."""
        if self._record_state is not None:
            return
        from . import gam_record_index
        self._records, state = gam_record_index.open_for(self.gam, self._source_stamp, self.index,
                                                         self._index_stamp, self.bins, self._record_mode)
        self._record_state = state
        self.cache_stats.update(state, record_fetches=0, fallback_fetches=0, selected_records=0,
                                record_bytes_read=0, preads=0)
        if self._records is not None or state["reason"] or state["mode"] != "auto":  # not: auto, no <gam>.gri
            how = f"record index {state['record_index']}" if self._records is not None else "GAI walk"
            print(f"GAM reader: {how} (PANSOMA_GAM_RECORD_INDEX={state['mode']}"
                  f"{': ' + state['reason'] if state['reason'] else ''})", file=sys.stderr, flush=True)

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

    def _postings(self, messages, metrics):
        """node -> indices of the records of `messages` visiting it (ascending); a record with MAPQ <= min_mapq
        is in none and its bytes are dropped (messages[i] = b"")."""
        postings = defaultdict(list)
        for i, raw in enumerate(messages):
            metrics["decoded_alignments"] += 1
            alignment = decode(raw)
            if self.min_mapq is not None and alignment.mapping_quality <= self.min_mapq:
                messages[i] = b""  # never yielded; its bytes need not stay cached
                continue
            for node in {m.position.node_id for m in alignment.path.mapping}:
                postings[node].append(i)
        return {node: tuple(indices) for node, indices in postings.items()}

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
        messages = stream.read_group()
        if messages is None:
            raise ValueError("GAI run extends past GAM EOF")
        postings = self._postings(messages, metrics)
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
        `wanted`, in file order: hits = those nodes. With a record index whose GAI map knows every run end
        point, the record path yields one group of the selected records (`_record_hits`); otherwise (no
        index, or bins edited after opening) the GAI walk reads the groups (`_walk_hits`)."""
        self.open_record_index()
        runs = self._ordinal_runs(ranges)
        if runs is not None:
            self.cache_stats["record_fetches"] += 1
            return self._record_hits(wanted, runs, metrics)
        if self._records is not None:
            self.cache_stats["fallback_fetches"] += 1
        return self._walk_hits(wanted, ranges, metrics)

    def _ordinal_runs(self, ranges):
        """[(first group, stop group)] of the merged runs, or None unless the record path applies: each run
        starts where the GAI walk reads a group start and ends where it stops cleanly (gam_record_index
        gai_map), after its start and after the previous run."""
        if self._records is None:
            return None
        runs = []
        for start, end in ranges:
            first, stop = self._records.ordinal(start, "start"), self._records.ordinal(end, "end")
            if first < 0 or stop <= first or (runs and first < runs[-1][1]):
                return None
            runs.append((first, stop))
        return runs

    def _record_hits(self, wanted, runs, metrics):
        """The record path: one group of the records in the groups [first, stop) of `runs` that may visit a node
        of `wanted` and pass min_mapq (gam_record_index.RecordIndex.select), in file order, through the same
        `_postings`; so the records and their order equal those of the GAI walk."""
        rows = self._records.select(np.unique(np.fromiter(wanted, dtype=np.uint64)), self.min_mapq, runs)
        metrics["groups"] += len(np.unique(rows["group"]))
        messages = read_records(self.gam, rows["start_vo"], rows["end_vo"], self.cache_stats)
        self.cache_stats["selected_records"] += len(messages)
        postings = self._postings(messages, metrics)
        hits = [n for n in postings if n in wanted]
        if hits:
            yield messages, postings, hits

    def _walk_hits(self, wanted, ranges, metrics):
        """The GAI walk: every group of `ranges`, indexed through the LRU cache."""
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
        with BgzfReader(self.gam) as stream:
            for start, end in ranges:
                stream.limit(end)
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

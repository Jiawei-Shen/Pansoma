"""Per-record index of a sorted long-read GAM (<gam>.gri): the records a node batch needs, without the GAI walk.

    python -m indexed_gam_pipeline_v4.gam_record_index build --gam G [--index I] [--processes 8] [--force]
    python -m indexed_gam_pipeline_v4.gam_record_index check --gam G [--index I] [--sample 2000] [--groups 8]
    python -m indexed_gam_pipeline_v4.gam_record_index info --gam G

A long-read GAM group holds 1,000 records (14-37 MB compressed) and vg's GAI files it under the smallest
node-ID bin holding all of them, so the GAI runs of one batch cover 40-100 groups (1-3 GB) for a few hundred
records. With <gam>.gri, IndexedGam reads only the records a batch can use (the record path of gam_reader):
those with a visited node in the batch, MAPQ above min_mapq and a group in the batch's merged GAI runs, i.e.
the records the GAI walk would return, then filters them with the same code.

File (little-endian; every section 64-byte aligned):
    b"PGRI\\0\\1\\r\\n", u64 header length, JSON header: format, gam (path, size, mtime_ns), index (path, size,
    mtime_ns, sha256), vg_pb2_sha256, protobuf, groups, records, long_records, W, sorted_by_lo1,
    checkpoint_stride, build, sections {name: {offset, count, dtype}}, trailer_offset
    groups       u64[groups + 1]: each group's start as BgzfReader.tell() gives it on a walk from the file
                 start (PARAMS_JSON and empty groups included), then the tell() at EOF
    gai_map      (offset u64, start i64, end i64), by offset: every distinct GAI run start and end; start = the
                 group a GAI walk from there reads first (any form of that group start), end = the group before
                 which a walk to there stops cleanly (exactly its start, or a form of EOF before EOF's tell():
                 groups); -1 where the walk would read a wrong group or raise
    records      56 B per record with at least one mapping, file order: start_vo, end_vo (tell() before its
                 length varint and after its last byte), group (u32), mapq (i32), lo1, hi1, lo2, hi2: its
                 visited node IDs ({m.position.node_id for m in path.mapping}) split at their largest gap
                 when that skips an ID (else lo2 = 1, hi2 = 0)
    checkpoints  u64: lo1 of every checkpoint_stride-th record; when lo1 never decreases (gamsort order)
    order        u64: the records by lo1 (stable); when it does decrease (fixtures, unsorted GAMs)
    long         64 B: every record whose span (largest - smallest node ID) exceeds W, and its row number;
                 W = the smallest power of two >= 2^12 that at most 1 % of the spans exceed (at most 2^26)
    trailer      sha256 of all bytes before it (`check`)

A batch with nodes b0..b1 reads the records with lo1 in [b0 - W, b1] (a window of the sorted table found
with the checkpoints and extended as batches move on) plus the long ones; every record visiting a batch node
is among them. The build walks the GAM once in parallel (segments cut at GAI run starts) and refuses to
write an index (exit 1, nothing written) unless every block inflates with matching CRC32 and ISIZE, every
record parses, the segment walks meet, and every GAI run starts at a group start and ends at a group start
or EOF after it: exactly the GAI runs the GAI walk reads without an error, so the record path selects the
same groups. IndexedGam uses <gam>.gri when its stamps match the GAM (size, mtime_ns), the GAI (size,
mtime_ns, sha256) and vg_pb2.py (sha256); PANSOMA_GAM_RECORD_INDEX: auto (default; else the GAI walk), off,
require (a missing or stale index raises), memory (the arrays scanned when the reader opens; tests).
"""
import argparse
from array import array
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time

from google.protobuf import __version__ as PROTOBUF_VERSION
from google.protobuf.message import DecodeError
import numpy as np

from . import vg_pb2
from .gam_reader import BgzfReader, decode, read_gai, read_records

FORMAT = "pansoma-gam-record-index-1"
MAGIC = b"PGRI\0\1\r\n"
ENVIRONMENT = "PANSOMA_GAM_RECORD_INDEX"
MODES = ("auto", "off", "require", "memory")
RECORD = np.dtype([("start_vo", "<u8"), ("end_vo", "<u8"), ("group", "<u4"), ("mapq", "<i4"),
                   ("lo1", "<u8"), ("hi1", "<u8"), ("lo2", "<u8"), ("hi2", "<u8")])
LONG = np.dtype(RECORD.descr + [("row", "<u8")])
GAI_MAP = np.dtype([("offset", "<u8"), ("start", "<i8"), ("end", "<i8")])
DTYPES = dict(groups=np.dtype("<u8"), gai_map=GAI_MAP, records=RECORD, checkpoints=np.dtype("<u8"),
              order=np.dtype("<u8"), long=LONG)
CHECKPOINT_STRIDE = 4096
W_BITS = (12, 26)
LONG_FRACTION = 0.01
SHORT_READ_BYTES = 2048  # mean compressed bytes per record below which build refuses without --force
SEGMENTS_PER_PROCESS = 4
ALIGN = 64


class Refused(ValueError):
    """The GAM cannot be indexed: the GAI walk is used (it reads and raises as it always did)."""


def gri_path(gam):
    """<gam>.gri next to the GAM file itself (its resolved path), wherever the GAM is opened from."""
    return Path(str(Path(gam).resolve()) + ".gri")


def vg_pb2_sha256():
    return hashlib.sha256(Path(vg_pb2.__file__).read_bytes()).hexdigest()


def _file_stamp(path):
    stat = Path(path).stat()
    return stat.st_size, stat.st_mtime_ns


def _reason(error):
    """An exception with the causes and contexts behind it, on one line."""
    parts = []
    while error is not None and len(parts) < 4:
        parts.append(f"{type(error).__name__}: {error}")
        error = error.__cause__ or error.__context__
    return " <- ".join(parts)


def default_processes():
    return int(os.environ.get("SLURM_CPUS_PER_TASK") or 8)


# --- scan ----------------------------------------------------------------------

def segments(size, bins, count):
    """discovery.segments from the GAI bins: up to `count` contiguous [start, end) segments cut at GAI run
    starts near equal file fractions (the last ends at EOF: None)."""
    starts = sorted({start for _, _, runs in bins for start, _ in runs} - {0})
    cuts, k = [0], 0
    for i in range(1, count):
        target = size * i // count
        while k < len(starts) and (starts[k] >> 16) < target:
            k += 1
        if k < len(starts) and starts[k] > cuts[-1]:
            cuts.append(starts[k])
    return list(zip(cuts, cuts[1:] + [None]))


def intervals(nodes):
    """(lo1, hi1, lo2, hi2) of sorted unique node IDs: split at the largest gap (the first of equal ones) when
    it skips an ID, else one interval (lo2 = 1, hi2 = 0)."""
    if len(nodes) > 1:
        gaps = [b - a for a, b in zip(nodes, nodes[1:])]
        i = max(range(len(gaps)), key=gaps.__getitem__)
        if gaps[i] > 1:
            return nodes[0], nodes[i], nodes[i + 1], nodes[-1]
    return nodes[0], nodes[-1], 1, 0


def _scan_segment(job):
    """Walk one segment of the GAM from `start` until the data position of `end` (the next segment's start;
    None = EOF): {start, starts (each group's tell(); the first as sought), stop (tell() there), rows (one
    RECORD per record with a mapping; group = index in the segment)}. Refused on anything unreadable."""
    gam, start, end = job
    starts, positions = array("Q"), []
    columns = [array(code) for code in ("Q", "Q", "I", "i", "Q", "Q", "Q", "Q")]
    try:
        with BgzfReader(gam, strict=True) as reader:
            target = None if end is None else reader.position(end)
            reader.seek(start)
            while True:
                here = reader.position()
                if here == target:
                    break
                if here is None or (target is not None and here > target):
                    raise Refused(f"the walk from {start} passed the next segment start {end}")
                starts.append(reader.tell())
                positions.clear()
                messages = reader.read_group(positions)
                if messages is None:
                    raise Refused(f"group at {starts[-1]} ends past EOF")
                for raw, (first, last) in zip(messages, positions):
                    alignment = decode(raw)
                    nodes = sorted({m.position.node_id for m in alignment.path.mapping})
                    if not nodes:
                        continue
                    if nodes[0] < 0:
                        raise Refused(f"record at {first} visits a negative node ID")
                    row = (first, last, len(starts) - 1, alignment.mapping_quality, *intervals(nodes))
                    for column, value in zip(columns, row):
                        column.append(value)
            stop = reader.tell()
            if end is None and reader.read_group() is not None:  # inflates (checks) the blocks after the data
                raise Refused("data after EOF")
    except Refused:
        raise
    except (ValueError, OSError, DecodeError) as error:
        raise Refused(f"segment from {start}: {_reason(error)}") from None
    rows = np.zeros(len(columns[0]), dtype=RECORD)
    for name, column in zip(RECORD.names, columns):
        rows[name] = np.frombuffer(column, dtype=column.typecode)
    return dict(start=start, starts=np.frombuffer(starts, dtype=np.uint64).copy(), stop=stop, rows=rows)


def scan(gam, bins, processes=1):
    """One walk over the GAM (`processes` workers over segments cut at GAI run starts): {groups, records,
    gai_map}. Refused unless every GAI run starts at a group start and ends at a group start or EOF after it
    (the runs the GAI walk reads without an error), every block and record reads, and the segments meet."""
    size = Path(gam).stat().st_size
    cuts = segments(size, bins, max(1, processes) * SEGMENTS_PER_PROCESS if processes > 1 else 1)
    jobs = [(str(gam), start, end) for start, end in cuts]
    if processes <= 1:
        results = [_scan_segment(job) for job in jobs]
    else:
        with ProcessPoolExecutor(processes, mp_context=multiprocessing.get_context("fork")) as pool:
            results = list(pool.map(_scan_segment, jobs))
    offsets = sorted({v for _, _, runs in bins for run in runs for v in run})
    with BgzfReader(gam, strict=True) as reader:
        # A segment's first group start is the previous walk's end as a walk from the file start gives it
        # (a GAI run start can be another form of the same position); blocks between the two are checked.
        starts, rows, canonical = [], [], 0
        for result in results:
            if result["start"] != canonical:
                try:
                    same = reader.position(canonical, check=True) == reader.position(result["start"])
                except OSError as error:
                    raise Refused(f"blocks before the segment from {result['start']}: {_reason(error)}") from None
                if not same:
                    raise Refused(f"the segment from {result['start']} does not start where the walk to it ended")
            if len(result["starts"]):
                result["starts"][0] = canonical
                canonical = result["stop"]
            result["rows"]["group"] += sum(len(s) for s in starts)
            starts.append(result["starts"])
            rows.append(result["rows"])
        groups = np.concatenate(starts + [np.array([canonical], dtype=np.uint64)])
        gai_map = map_offsets(reader, groups, offsets)
    lookup = {int(o): (int(s), int(e)) for o, s, e in gai_map.tolist()}
    for _, _, runs in bins:
        for start, end in runs:
            first, stop = lookup[start][0], lookup[end][1]
            if first < 0:
                raise Refused(f"GAI run ({start}, {end}) does not start at a GAM group start")
            if stop < 0:
                raise Refused(f"GAI run ({start}, {end}) does not end at a GAM group start or EOF")
            if stop <= first:
                raise Refused(f"GAI run ({start}, {end}) ends before the group it starts with")
    return dict(groups=groups, records=np.concatenate(rows) if rows else np.zeros(0, RECORD), gai_map=gai_map)


def map_offsets(reader, groups, offsets):
    """gai_map of the GAI offsets: group ordinals of a GAI walk starting / ending there (see the module doc)."""
    n = len(groups) - 1
    values = np.array(offsets, dtype=np.uint64)
    k = np.searchsorted(groups, values)
    exact = (k <= n) & (groups[np.minimum(k, n)] == values)
    result = np.zeros(len(values), dtype=GAI_MAP)
    result["offset"] = values
    result["start"] = np.where(exact & (k < n), k, -1)
    result["end"] = np.where(exact, k, -1)
    for i in np.flatnonzero(~exact).tolist():
        value, j = offsets[i], int(k[i])
        try:
            here = reader.position(value)
        except OSError:
            continue  # the walk would raise there
        if here is None:  # a form of EOF: a walk to it ends after the last group (equivalent_eof)
            if (n == 0 or int(groups[n - 1]) < value) and value < int(groups[n]):
                result["end"][i] = n
        else:  # another form of a group start: a walk reads that group first, but cannot stop there
            for c in (j - 1, j):
                if 0 <= c < n and reader.position(int(groups[c])) == here:
                    result["start"][i] = c
    return result


def layout(records):
    """(W, row numbers of the long records, sorted_by_lo1) of a records table."""
    hi = np.where(records["lo2"] <= records["hi2"], records["hi2"], records["hi1"])
    span = hi - records["lo1"]
    bits = W_BITS[0]
    while bits < W_BITS[1] and np.count_nonzero(span > (1 << bits)) > LONG_FRACTION * len(records):
        bits += 1
    lo1 = records["lo1"]
    return 1 << bits, np.flatnonzero(span > (1 << bits)), bool(np.all(lo1[1:] >= lo1[:-1]))


def long_table(records, rows):
    table = np.zeros(len(rows), dtype=LONG)
    for name in RECORD.names:
        table[name] = records[name][rows]
    table["row"] = rows
    return table


# --- the index as the reader uses it ---------------------------------------------

def _hit(nodes, lo, hi):
    """Per interval [lo, hi]: does it hold a node of `nodes` (sorted unique)? ranges()' test."""
    first = np.searchsorted(nodes, lo, side="left")
    inside = first < len(nodes)
    hit = np.zeros(len(lo), dtype=bool)
    hit[inside] = nodes[first[inside]] <= hi[inside]
    return hit


class RecordIndex:
    """A record index: a .gri (`load`; a sorted file's records read in lo1 windows on demand) or a scan held in
    memory (`from_scan`). select() picks the rows of one batch; ordinal() maps GAI offsets to groups."""

    def __init__(self, header, groups, gai_map, long, records=None, order=None, checkpoints=None, fd=None):
        self.header, self.W, self.groups, self._map = header, header["W"], groups, gai_map
        self._long = np.zeros(len(long), dtype=RECORD)
        for name in RECORD.names:
            self._long[name] = long[name]
        self._long_rows = long["row"].astype(np.int64)
        self._records, self._order, self._checkpoints, self._fd = records, order, checkpoints, fd
        if records is not None:
            self._lo1 = records["lo1"] if order is None else records["lo1"][order]
        self._cache = (0, np.zeros(0, dtype=RECORD))

    @classmethod
    def from_scan(cls, scanned):
        records = scanned["records"]
        width, rows, ordered = layout(records)
        header = dict(format=FORMAT, groups=len(scanned["groups"]) - 1, records=len(records), long_records=len(rows),
                      W=width, sorted_by_lo1=ordered, checkpoint_stride=CHECKPOINT_STRIDE)
        order = None if ordered else np.argsort(records["lo1"], kind="stable").astype(np.uint64)
        return cls(header, scanned["groups"], scanned["gai_map"], long_table(records, rows), records, order)

    def close(self):
        fd, self._fd = getattr(self, "_fd", None), None
        if fd is not None:
            os.close(fd)

    __del__ = close

    def ordinal(self, offset, role):
        """gai_map[role] ('start' or 'end') of a virtual offset; -1 when the GAI has no such offset."""
        offsets = self._map["offset"]
        i = int(np.searchsorted(offsets, offset))
        return int(self._map[role][i]) if i < len(offsets) and int(offsets[i]) == offset else -1

    def _rows(self, start, stop):
        """Records [start, stop) of a sorted file, through a window kept from `start` on (batches move forward)."""
        first, cached = self._cache
        if first <= start and stop <= first + len(cached):
            return cached[start - first:stop - first]
        if first <= start < first + len(cached):
            cached = np.concatenate([cached[start - first:], self._read(first + len(cached), stop)])
        else:
            cached = self._read(start, stop)
        self._cache = (start, cached)
        return cached[:stop - start]

    def _read(self, start, stop):
        spec = self.header["sections"]["records"]
        return _section_bytes(self._fd, spec["offset"] + start * RECORD.itemsize, stop - start, RECORD)

    def window(self, low, high):
        """(rows, row numbers) of the records with lo1 in [low, high], ascending lo1."""
        if self._records is not None:
            i, j = np.searchsorted(self._lo1, low, side="left"), np.searchsorted(self._lo1, high, side="right")
            numbers = np.arange(i, j) if self._order is None else self._order[i:j].astype(np.int64)
            return self._records[numbers], numbers
        stride, count = self.header["checkpoint_stride"], self.header["records"]
        start = max(0, int(np.searchsorted(self._checkpoints, low, side="left")) - 1) * stride
        stop = min(count, int(np.searchsorted(self._checkpoints, high, side="right")) * stride)
        rows = self._rows(start, stop)
        i, j = np.searchsorted(rows["lo1"], low, side="left"), np.searchsorted(rows["lo1"], high, side="right")
        return rows[i:j], np.arange(start + i, start + j)

    def select(self, nodes, min_mapq, runs):
        """Rows (file order) of the records that may visit a node of `nodes` (sorted unique uint64): an
        interval holds one; MAPQ > min_mapq (unless None); group in [first, stop) of a run of `runs`."""
        if not len(nodes) or not runs:
            return np.zeros(0, dtype=RECORD)
        low, high = int(nodes[0]), int(nodes[-1])
        rows, numbers = self.window(max(0, low - self.W), high)
        near = self._long["lo1"] <= high
        rows = np.concatenate([rows, self._long[near]])
        numbers, first = np.unique(np.concatenate([numbers, self._long_rows[near]]), return_index=True)
        rows = rows[first]
        keep = _hit(nodes, rows["lo1"], rows["hi1"]) | _hit(nodes, rows["lo2"], rows["hi2"])
        if min_mapq is not None:
            keep &= rows["mapq"] > min_mapq
        group = rows["group"].astype(np.int64)
        firsts, stops = (np.array(column, dtype=np.int64) for column in zip(*runs))
        at = np.searchsorted(firsts, group, side="right") - 1
        keep &= (at >= 0) & (group < stops[np.maximum(at, 0)])
        return rows[keep]


# --- file ----------------------------------------------------------------------

def _align(offset):
    return -(-offset // ALIGN) * ALIGN


def _section_bytes(fd, offset, count, dtype):
    size = count * dtype.itemsize
    data = b""
    while len(data) < size:
        part = os.pread(fd, size - len(data), offset + len(data))
        if not part:
            raise ValueError("record index truncated")
        data += part
    return np.frombuffer(data, dtype=dtype)


def write(path, header, sections):
    """MAGIC, header length, header JSON (offsets of the 64-byte aligned sections filled in), sections, sha256
    trailer: into <path>.tmp (created exclusively), fsynced, then renamed to `path`."""
    length = 0
    while True:
        offset, table = _align(16 + length), {}
        for name, values in sections:
            table[name] = dict(offset=offset, count=len(values), dtype=np.lib.format.dtype_to_descr(values.dtype))
            offset = _align(offset + values.nbytes)
        header.update(sections=table, trailer_offset=offset)
        text = json.dumps(header, indent=1).encode()
        if len(text) <= length:
            break
        length = len(text) + ALIGN
    temporary = Path(str(path) + ".tmp")
    fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        with os.fdopen(fd, "wb") as stream:
            digest, written = hashlib.sha256(), 0

            def put(data):
                nonlocal written
                stream.write(data)
                digest.update(data)
                written += len(data)

            put(MAGIC + length.to_bytes(8, "little") + text.ljust(length))
            for name, values in sections:
                put(bytes(table[name]["offset"] - written))
                put(np.ascontiguousarray(values).reshape(-1).view(np.uint8))
            put(bytes(header["trailer_offset"] - written))
            stream.write(digest.digest())
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def read_header(fd):
    head = os.pread(fd, 16, 0)
    if len(head) != 16 or head[:8] != MAGIC:
        raise ValueError("not a GAM record index (no PGRI magic)")
    length = int.from_bytes(head[8:], "little")
    text = os.pread(fd, length, 16)
    if len(text) != length:
        raise ValueError("record index truncated")
    header = json.loads(text)
    if header.get("format") != FORMAT:
        raise ValueError(f"unknown record index format {header.get('format')!r}")
    if os.fstat(fd).st_size != header["trailer_offset"] + 32:
        raise ValueError("record index truncated or overlong")
    for name, spec in header["sections"].items():
        if np.lib.format.descr_to_dtype([tuple(f) for f in spec["dtype"]] if isinstance(spec["dtype"], list)
                                        else spec["dtype"]) != DTYPES[name]:
            raise ValueError(f"unexpected dtype of section {name}")
    return header


def section(fd, header, name):
    spec = header["sections"][name]
    return _section_bytes(fd, spec["offset"], spec["count"], DTYPES[name])


def load(path):
    """RecordIndex of a .gri: header, size and dtypes checked (the trailer only by `check`); a sorted file's
    records stay on disk, an unsorted one is read whole with its order."""
    fd = os.open(str(path), os.O_RDONLY)
    try:
        header = read_header(fd)
        parts = {name: section(fd, header, name) for name in header["sections"] if name != "records"}
        if header["sorted_by_lo1"]:
            return RecordIndex(header, parts["groups"], parts["gai_map"], parts["long"],
                               checkpoints=parts["checkpoints"], fd=fd)
        records = section(fd, header, "records")
        os.close(fd)
        return RecordIndex(header, parts["groups"], parts["gai_map"], parts["long"], records, parts["order"])
    except BaseException:
        try:
            os.close(fd)
        except OSError:
            pass
        raise


def stale(header, gam_stamp, index_stamp):
    """Why an index header does not describe this GAM (size, mtime_ns), GAI (size, mtime_ns, sha256) and
    vg_pb2.py; [] when it does."""
    problems = []
    if (header["gam"]["size"], header["gam"]["mtime_ns"]) != tuple(gam_stamp):
        problems.append("the GAM's size or mtime differs")
    if (header["index"]["size"], header["index"]["mtime_ns"], header["index"]["sha256"]) != tuple(index_stamp):
        problems.append("the GAI's size, mtime or sha256 differs")
    if header["vg_pb2_sha256"] != vg_pb2_sha256():
        problems.append("vg_pb2.py differs")
    return problems


def current_stamps(gam, index):
    """(GAM stamp, GAI stamp) as IndexedGam compares them with an index header."""
    raw = Path(index).read_bytes()
    return _file_stamp(gam), (*_file_stamp(index), hashlib.sha256(raw).hexdigest())


def validate(path, gam, index=None):
    """Problems of <gam>.gri `path` against the GAM and GAI as they are now ([] = usable)."""
    gam_stamp, index_stamp = current_stamps(gam, index or str(gam) + ".gai")
    fd = os.open(str(path), os.O_RDONLY)
    try:
        return stale(read_header(fd), gam_stamp, index_stamp)
    except ValueError as error:
        return [str(error)]
    finally:
        os.close(fd)


def open_for(gam, gam_stamp, index, index_stamp, bins):
    """(RecordIndex or None, {mode, record_index, reason}) of IndexedGam, by PANSOMA_GAM_RECORD_INDEX:
    auto: <gam>.gri when it matches the GAM, the GAI and vg_pb2.py, else None (the GAI walk); off: None;
    require: the matching <gam>.gri, else ValueError; memory: the same arrays from a scan now (None when the
    scan refuses the GAM). reason: why an index there is not used (None when none is there)."""
    mode = os.environ.get(ENVIRONMENT) or "auto"
    if mode not in MODES:
        raise ValueError(f"{ENVIRONMENT} must be one of {', '.join(MODES)}, not {mode!r}")
    if mode == "off":
        return None, dict(mode=mode, record_index=None, reason=None)
    if mode == "memory":
        try:
            return RecordIndex.from_scan(scan(gam, bins)), dict(mode=mode, record_index="memory", reason=None)
        except Refused as error:
            return None, dict(mode=mode, record_index=None, reason=f"not indexable: {error}")
    path, records, problems = gri_path(gam), None, []
    if not path.exists():
        if mode == "require":
            raise ValueError(f"{ENVIRONMENT}=require, but there is no {path}; build it with "
                             f"python -m {__package__}.gam_record_index build --gam {gam}")
        return None, dict(mode=mode, record_index=None, reason=None)
    try:
        records = load(path)
        problems = stale(records.header, gam_stamp, index_stamp)
    except (OSError, ValueError, KeyError) as error:
        problems = [f"unreadable ({error})"]
    if not problems:
        return records, dict(mode=mode, record_index=str(path), reason=None)
    if records is not None:
        records.close()
    reason = f"{path.name} is stale: {'; '.join(problems)}"
    if mode == "require":
        raise ValueError(f"{ENVIRONMENT}=require, but {reason} ({path}); rebuild it with "
                         f"python -m {__package__}.gam_record_index build --gam {gam}, or move it aside")
    return None, dict(mode=mode, record_index=None, reason=reason)


# --- commands ------------------------------------------------------------------

def mean_record_bytes(gam, bins, points=4):
    """Mean compressed bytes per record of the first GAM group with records after `points` GAI run starts spread
    over the file (block-granular, so short reads come out near 0); None when no group has records."""
    total = count = 0
    with BgzfReader(gam) as reader:
        for start, _ in segments(Path(gam).stat().st_size, bins, points):
            reader.seek(start)
            positions = []
            while not positions and reader.read_group(positions) is not None:
                pass
            if positions:
                total += (positions[-1][1] >> 16) - (positions[0][0] >> 16)
                count += len(positions)
    return total / count if count else None


def build(gam, index=None, processes=8, force=False):
    """Write <gam>.gri (never over an existing one); Refused (nothing written) as described in the module doc,
    and for a short-read GAM unless `force`. Returns a summary."""
    started = time.perf_counter()
    target, index = gri_path(gam), Path(index or str(gam) + ".gai")
    gam, index = Path(gam).resolve(), index.resolve()
    if target.exists():
        raise ValueError(f"{target} exists; move it aside to rebuild it")
    gam_stamp = _file_stamp(gam)
    _, bins, sha256 = read_gai(index, gam_stamp[0])
    index_stamp = (*_file_stamp(index), sha256)
    mean = mean_record_bytes(gam, bins)
    if not force and mean is not None and mean < SHORT_READ_BYTES:
        raise Refused(f"{mean:.0f} compressed bytes per record: a short-read GAM, whose GAI walk reads little "
                      f"more than its records; --force indexes it anyway")
    scanned = scan(gam, bins, processes)
    if _file_stamp(gam) != gam_stamp or (*_file_stamp(index), sha256) != index_stamp:
        raise Refused("the GAM or GAI changed while it was indexed")
    records, groups = scanned["records"], scanned["groups"]
    width, rows, ordered = layout(records)
    header = dict(format=FORMAT, created=time.strftime("%Y-%m-%dT%H:%M:%S"),
                  gam=dict(path=str(gam), size=gam_stamp[0], mtime_ns=gam_stamp[1]),
                  index=dict(path=str(index), size=index_stamp[0], mtime_ns=index_stamp[1], sha256=sha256),
                  vg_pb2_sha256=vg_pb2_sha256(), protobuf=PROTOBUF_VERSION, groups=len(groups) - 1,
                  records=len(records), long_records=len(rows), W=width, sorted_by_lo1=ordered,
                  checkpoint_stride=CHECKPOINT_STRIDE,
                  build=dict(processes=processes, mean_record_bytes=mean,
                             seconds=round(time.perf_counter() - started, 1)))
    sections = [("groups", groups), ("gai_map", scanned["gai_map"]), ("records", records)]
    if ordered:
        sections.append(("checkpoints", np.ascontiguousarray(records["lo1"][::CHECKPOINT_STRIDE])))
    else:
        sections.append(("order", np.argsort(records["lo1"], kind="stable").astype(np.uint64)))
    sections.append(("long", long_table(records, rows)))
    write(target, header, sections)
    return dict(path=str(target), bytes=target.stat().st_size, groups=header["groups"], records=len(records),
                long_records=len(rows), W=width, sorted_by_lo1=ordered, seconds=round(time.perf_counter() - started, 1))


def check(gam, index=None, sample=2000, groups=8, seed=0):
    """Report of <gam>.gri with its `problems` ([] = good): the trailer sha256, the stamps, the structure, the
    GAI map against the GAI, `sample` random records re-read from the GAM (MAPQ, intervals and group bounds as
    the scan computes them) and `groups` random groups walked again (every record's offsets and fields)."""
    started = time.perf_counter()
    gam, index = Path(gam), Path(index or str(gam) + ".gai")
    path = gri_path(gam)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        size = os.fstat(stream.fileno()).st_size
        left = size - 32
        while left > 0:
            chunk = stream.read(min(left, 64 << 20))
            if not chunk:
                break
            digest.update(chunk)
            left -= len(chunk)
        trailer = stream.read(32)
    problems = [] if trailer == digest.digest() else ["trailer sha256 differs (the file is damaged)"]
    fd = os.open(str(path), os.O_RDONLY)
    try:
        header = read_header(fd)
        problems += stale(header, *current_stamps(gam, index))
        parts = {name: section(fd, header, name) for name in header["sections"]}
    finally:
        os.close(fd)
    records, starts, gai_map, n = parts["records"], parts["groups"], parts["gai_map"], header["groups"]

    def need(condition, problem):
        if not condition:
            problems.append(problem)

    need(len(starts) == n + 1 and bool(np.all(starts[1:] > starts[:-1])), "group starts not increasing")
    need(len(records) == header["records"], "record count differs from the header")
    if len(records):
        group = records["group"].astype(np.int64)
        need(bool(np.all(np.diff(group) >= 0)) and int(group[-1]) < n, "record groups out of order or range")
        need(bool(np.all(records["start_vo"] < records["end_vo"]))
             and bool(np.all(records["end_vo"][:-1] <= records["start_vo"][1:])), "record offsets not increasing")
        need(bool(np.all(records["start_vo"] > starts[group])) and bool(np.all(records["end_vo"] <= starts[group + 1])),
             "a record outside its group")
        single = (records["lo2"] == 1) & (records["hi2"] == 0)
        need(bool(np.all(records["lo1"] <= records["hi1"]))
             and bool(np.all(single | ((records["hi1"] + 1 < records["lo2"]) & (records["lo2"] <= records["hi2"])))),
             "invalid node intervals")
    width, rows, ordered = layout(records)
    need((width, ordered, len(rows)) == (header["W"], header["sorted_by_lo1"], header["long_records"]),
         "W, sorted_by_lo1 or long_records differ from the records")
    need(np.array_equal(parts["long"], long_table(records, rows)), "long table differs from the records")
    if ordered:
        need(np.array_equal(parts["checkpoints"], records["lo1"][::header["checkpoint_stride"]]), "checkpoints differ")
    else:
        need(np.array_equal(parts["order"], np.argsort(records["lo1"], kind="stable")), "order differs")
    _, bins, _ = read_gai(index, Path(gam).stat().st_size)
    offsets = sorted({v for _, _, runs in bins for run in runs for v in run})
    need(gai_map["offset"].tolist() == offsets, "gai_map offsets differ from the GAI's")
    lookup = {int(o): (int(s), int(e)) for o, s, e in gai_map.tolist()}
    need(all(0 <= lookup.get(s, (-1, -1))[0] < lookup.get(e, (-1, -1))[1] for _, _, runs in bins for s, e in runs),
         "a GAI run the record path cannot take")
    rng = np.random.default_rng(seed)
    picked = (np.sort(rng.choice(len(records), min(sample, len(records)), replace=False)) if len(records)
              else np.zeros(0, dtype=np.int64))
    messages = read_records(gam, records["start_vo"][picked], records["end_vo"][picked],
                            dict(preads=0, record_bytes_read=0))
    for row, raw in zip(records[picked], messages):
        alignment = decode(raw)
        nodes = sorted({m.position.node_id for m in alignment.path.mapping})
        fields = tuple(int(row[name]) for name in ("mapq", "lo1", "hi1", "lo2", "hi2"))
        need(nodes and (alignment.mapping_quality, *intervals(nodes)) == fields,
             f"record at {int(row['start_vo'])} differs from the GAM")
    walked = np.sort(rng.choice(n, min(groups, n), replace=False)) if n else np.zeros(0, dtype=np.int64)
    with BgzfReader(gam, strict=True) as reader:
        for g in walked.tolist():
            reader.seek(int(starts[g]))
            positions = []
            messages = reader.read_group(positions)
            found = []
            for raw, (first, last) in zip(messages or [], positions):
                alignment = decode(raw)
                nodes = sorted({m.position.node_id for m in alignment.path.mapping})
                if nodes:
                    found.append((first, last, g, alignment.mapping_quality, *intervals(nodes)))
            need(reader.tell() == int(starts[g + 1]), f"group {g} does not end at the next group start")
            need(found == [tuple(int(v) for v in r) for r in records[records["group"] == g].tolist()],
                 f"group {g}'s records differ from the GAM")
    return dict(path=str(path), bytes=size, problems=problems, sampled_records=len(picked), walked_groups=len(walked),
                seconds=round(time.perf_counter() - started, 1))


def info(gam):
    fd = os.open(str(gri_path(gam)), os.O_RDONLY)
    try:
        return dict(read_header(fd), bytes=os.fstat(fd).st_size)
    finally:
        os.close(fd)


def make_parser():
    parser = argparse.ArgumentParser(prog=f"python -m {__package__}.gam_record_index", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("build", help="write <gam>.gri (refused for a short-read GAM unless --force)")
    p.add_argument("--gam", required=True, help="sorted BGZF GAM with its GAI")
    p.add_argument("--index", help="default: GAM path + .gai")
    p.add_argument("--processes", type=int, default=default_processes(),
                   help="parallel walkers (default: $SLURM_CPUS_PER_TASK, else 8)")
    p.add_argument("--force", action="store_true", help="index a short-read GAM too")
    p = commands.add_parser("check", help="verify <gam>.gri against the GAM and GAI; exit 1 on a problem")
    p.add_argument("--gam", required=True)
    p.add_argument("--index", help="default: GAM path + .gai")
    p.add_argument("--sample", type=int, default=2000, help="random records re-read from the GAM")
    p.add_argument("--groups", type=int, default=8, help="random groups walked again")
    p = commands.add_parser("info", help="print the header of <gam>.gri")
    p.add_argument("--gam", required=True)
    return parser


def main(argv=None):
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "build":
            result = build(args.gam, args.index, args.processes, args.force)
        elif args.command == "check":
            result = check(args.gam, args.index, args.sample, args.groups)
        else:
            result = info(args.gam)
    except (ValueError, OSError) as error:
        parser.exit(1, f"Error: {_reason(error)}\n")
    print(json.dumps(result, indent=2))
    if result.get("problems"):
        sys.exit(1)


if __name__ == "__main__":
    main()

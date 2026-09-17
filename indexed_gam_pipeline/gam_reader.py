"""Read VG GAI v0/v1 bins and seek BGZF GAM groups without vg or an NPU store.

Format reference: https://github.com/vgteam/vg/blob/master/src/stream_index.cpp
Queries conservatively scan intersecting bins (the window lower-bound speedup
is intentionally unnecessary here), merge their virtual-offset runs, then
filter complete Alignment messages by exact node membership.
"""
import bisect
import gzip
import io
from pathlib import Path

import pysam
import vg_pb2


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
    """One tagged GAM group; None denotes clean EOF."""
    count = varint(stream, allow_eof=True)
    if count is None:
        return None
    if count == 0:
        return []
    tag = exact(stream, varint(stream))
    if tag != b"GAM":
        raise ValueError(f"Expected tagged GAM group, got {tag[:30]!r}")
    return [exact(stream, varint(stream)) for _ in range(count - 1)]


def decode(raw):
    alignment = vg_pb2.Alignment()
    alignment.ParseFromString(raw)
    return alignment


def scan_gam(path, max_alignments=None):
    """Independent sequential scan; does not consult the GAI."""
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


class IndexedGam:
    def __init__(self, gam, index=None):
        self.gam = Path(gam)
        self.index = Path(index or str(gam) + ".gai")
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
        self.bins = []
        gam_size = self.gam.stat().st_size
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
        self.window_count = varint(data)
        for _ in range(self.window_count):
            varint(data)
            varint(data)
        if data.read(1):
            raise ValueError("Unexpected trailing GAI data")

    def ranges(self, nodes):
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

    def fetch(self, nodes, metrics=None):
        wanted = set(nodes)
        ranges = self.ranges(wanted)
        if metrics is None:
            metrics = {}
        metrics.update(runs=len(ranges), groups=0, decoded_alignments=0,
                       returned_alignments=0)
        with pysam.BGZFile(str(self.gam), "rb") as stream:
            for start, end in ranges:
                stream.seek(start)
                while stream.tell() < end:
                    messages = group(stream)
                    if messages is None:
                        raise ValueError("GAI run extends past GAM EOF")
                    metrics["groups"] += 1
                    for raw in messages:
                        metrics["decoded_alignments"] += 1
                        alignment = decode(raw)
                        if any(m.position.node_id in wanted for m in alignment.path.mapping):
                            metrics["returned_alignments"] += 1
                            yield alignment
                if stream.tell() != end:
                    raise ValueError("GAI run does not end on a GAM group boundary")

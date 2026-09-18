"""Persistent node -> distinct GFA W-record count, built by a single GFA scan."""
import gzip
import hashlib
import json
import uuid
from pathlib import Path
import re
import sqlite3
import time

import numpy as np

CACHE_VERSION = 'gfa-distinct-w-records-v1'
MAX_COUNT = np.iinfo(np.int32).max


def walk_nodes(walk):
    """Numeric vg node IDs; strand and repeated visits do not multiply counts."""
    if not walk or walk[:1] not in (b'>', b'<') or re.search(rb'[^0-9<>]', walk):
        raise ValueError('W walk must use >/< followed by positive numeric vg node IDs')
    nodes = np.fromstring(walk.replace(b'>', b' ').replace(b'<', b' '), sep=' ', dtype=np.int64)
    if len(nodes) != walk.count(b'>') + walk.count(b'<') or np.any(nodes <= 0):
        raise ValueError('Invalid or out-of-range node ID in W walk')
    # NumPy saturates decimal integers above int64 max instead of raising.
    # Validate those rare boundary cases before accepting the parsed IDs.
    if np.any(nodes == np.iinfo(np.int64).max):
        if any(int(token.group()) > np.iinfo(np.int64).max for token in re.finditer(rb'\d+', walk)):
            raise ValueError('Node ID exceeds signed int64 range')
    return np.unique(nodes)


class _Accumulator:
    """Vectorized disk-backed counts for normal vg IDs; SQLite for sparse huge IDs."""
    DENSE_LIMIT = 100_000_000

    def __init__(self, db, path):
        self.db, self.path = db, path
        self.array = None
        self.capacity = 0

    def add(self, ids):
        cut = np.searchsorted(ids, self.DENSE_LIMIT, side='right')
        low = ids[:cut]
        if len(low):
            maximum = int(low[-1])
            if maximum >= self.capacity:
                if self.array is not None:
                    self.array.flush()
                    self.array._mmap.close()
                self.capacity = min(self.DENSE_LIMIT+1, ((maximum//1_000_000)+1)*1_000_000)
                with open(self.path, 'ab') as stream:
                    stream.truncate(self.capacity*4)
                self.array = np.memmap(self.path, dtype=np.uint32, mode='r+', shape=(self.capacity,))
            self.array[low] += 1
        self.db.executemany('INSERT INTO counts VALUES (?, 1) ON CONFLICT(node_id) '
                            'DO UPDATE SET walk_count=walk_count+1', ((int(n),) for n in ids[cut:]))

    def finish(self):
        if self.array is None:
            return
        for start in range(0, self.capacity, 1_000_000):
            values = self.array[start:start+1_000_000]
            present = np.flatnonzero(values)
            self.db.executemany('INSERT INTO counts VALUES (?, ?)',
                                ((start+int(i), int(values[i])) for i in present))
            self.db.commit()
            if start % 10_000_000 == 0:
                print(f'Finalizing walk lookup: node IDs through {min(start+1_000_000, self.capacity):,}', flush=True)

    def close(self):
        if self.array is not None:
            self.array._mmap.close()
            self.array = None
        if self.path.exists():
            self.path.unlink()


def build_walk_counts(gfa, output):
    """Atomic reusable SQLite cache. Duplicate first-seven-field W records count once.

    Optional tags are not part of W identity. Distinct coordinate intervals,
    samples, haplotypes, sequences or walk strings are distinct records. Only W
    records contribute; S, L and P records do not. Nodes absent from W have zero.
    """
    source = Path(gfa).resolve()
    output = Path(output)
    if output.exists():
        raise ValueError(f'Walk-count cache already exists: {output}')
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f'.building-{uuid.uuid4().hex}')
    if temporary.exists():
        raise ValueError(f'Incomplete cache already exists: {temporary}')
    before = source.stat()
    start = time.perf_counter()
    digest = hashlib.sha256()
    total_lines = w_lines = unique_w = duplicate_w = empty_w = 0
    opener = gzip.open if str(source).endswith('.gz') else open
    accumulator = None
    try:
        with sqlite3.connect(temporary) as db:
            db.execute('PRAGMA journal_mode=OFF')
            db.execute('PRAGMA synchronous=OFF')
            db.execute('PRAGMA cache_size=-131072')
            db.execute('CREATE TABLE counts (node_id INTEGER PRIMARY KEY, walk_count INTEGER NOT NULL)')
            db.execute('CREATE TABLE seen_records (digest BLOB PRIMARY KEY) WITHOUT ROWID')
            db.execute('CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)')
            accumulator = _Accumulator(db, temporary.with_name(temporary.name+'.counts'))
            with opener(source, 'rb') as stream:
                for line in stream:
                    total_lines += 1
                    digest.update(line)
                    if total_lines % 1000000 == 0:
                        print(f'GFA scan: {total_lines:,} lines, {unique_w:,} distinct W records, '
                              f'{time.perf_counter()-start:.1f}s', flush=True)
                    if not line.startswith(b'W\t'):
                        continue
                    w_lines += 1
                    fields = line.rstrip(b'\r\n').split(b'\t', 7)
                    if len(fields) < 7:
                        raise ValueError(f'Malformed W record at GFA line {total_lines}')
                    identity = hashlib.sha256(b'\t'.join(fields[:7])).digest()
                    inserted = db.execute('INSERT OR IGNORE INTO seen_records VALUES (?)', (identity,)).rowcount
                    if not inserted:
                        duplicate_w += 1
                        continue
                    ids = walk_nodes(fields[6])
                    unique_w += 1
                    empty_w += int(len(ids) == 0)
                    if unique_w > MAX_COUNT:
                        raise ValueError('Too many distinct W records for int32 counts')
                    accumulator.add(ids)
                    db.commit()
                    if unique_w <= 5 or unique_w % 100 == 0:
                        print(f'W record {unique_w:,}: {len(ids):,} distinct nodes; '
                              f'{time.perf_counter()-start:.1f}s elapsed', flush=True)
            after = source.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError('GFA changed while the cache was being built')
            accumulator.finish()
            node_count, max_count = db.execute('SELECT count(*), coalesce(max(walk_count),0) FROM counts').fetchone()
            if max_count > MAX_COUNT:
                raise ValueError('Walk counts exceed the candidate-v3 int32 range')
            if not w_lines:
                raise ValueError('GFA contains no W records; cannot build a walk-count channel')
            report = dict(cache_format_version=CACHE_VERSION, status='complete', source_gfa=str(source),
                source_size_bytes=before.st_size, source_mtime_ns=before.st_mtime_ns,
                uncompressed_gfa_sha256=digest.hexdigest(), lines_scanned=total_lines,
                w_records=w_lines, distinct_w_records=unique_w, duplicate_w_records=duplicate_w,
                empty_w_records=empty_w, nodes_in_walks=node_count, max_node_walk_count=max_count,
                elapsed_seconds=time.perf_counter()-start,
                identity='W first seven tab-separated fields; optional tags ignored',
                counting='one contribution per distinct W record per node; orientation and repeated visits ignored')
            db.execute('INSERT INTO metadata VALUES (?, ?)', ('report', json.dumps(report)))
            db.commit()
        temporary.rename(output)
        return report
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise
    finally:
        if accumulator is not None:
            accumulator.close()


class WalkCounts:
    def __init__(self, path):
        self.path = Path(path).resolve()
        if not self.path.is_file():
            raise ValueError(f'Walk-count cache does not exist: {self.path}')
        self.db = sqlite3.connect(self.path.as_uri() + '?mode=ro', uri=True)
        try:
            row = self.db.execute("SELECT value FROM metadata WHERE key='report'").fetchone()
            self.metadata = json.loads(row[0]) if row else {}
            if self.metadata.get('cache_format_version') != CACHE_VERSION or self.metadata.get('status') != 'complete':
                raise ValueError('Incomplete or unsupported walk-count cache')
            source = Path(self.metadata['source_gfa'])
            if source.exists():
                stat = source.stat()
                if (stat.st_size, stat.st_mtime_ns) != (self.metadata['source_size_bytes'], self.metadata['source_mtime_ns']):
                    raise ValueError('Walk-count cache is stale: source GFA size/mtime changed')
        except BaseException:
            self.db.close()
            raise

    def get_counts(self, nodes):
        ids = sorted(set(int(n) for n in nodes))
        result = dict.fromkeys(ids, 0)
        for offset in range(0, len(ids), 900):
            chunk = ids[offset:offset+900]
            sql = 'SELECT node_id, walk_count FROM counts WHERE node_id IN (' + ','.join('?' for _ in chunk) + ')'
            for node, count in self.db.execute(sql, chunk):
                if not 0 <= count <= MAX_COUNT:
                    raise ValueError(f'Invalid walk count at node {node}: {count}')
                result[node] = count
        return result

    def close(self):
        self.db.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gfa', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(build_walk_counts(args.gfa, args.output), indent=2))

"""Persistent, graph-identified GBWT path counts with node-sequence validation."""
import hashlib
import json
from pathlib import Path
import sqlite3
import subprocess

METRIC = 'distinct-gbwt-paths-v1'


def fingerprint(path):
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


class GBZCounts:
    def __init__(self, gbz, helper, cache):
        self.source = Path(gbz).resolve()
        self.helper = str(Path(helper).resolve())
        self.path = Path(cache).resolve()
        self.process = None
        self.source_stat = fingerprint(self.source)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.path)
        try:
            tables = {r[0] for r in self.db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if tables:
                if not {'gbz_metadata', 'gbz_counts'} <= tables:
                    raise ValueError('Not a GBZ occurrence cache; GFA W caches cannot be reused')
                self.metadata = json.loads(self.db.execute('SELECT value FROM gbz_metadata').fetchone()[0])
                if self.metadata.get('metric') != METRIC or self.metadata.get('source') != self.source_stat:
                    raise ValueError('GBZ occurrence cache source/metric mismatch; use a new cache')
            else:
                digest = hashlib.sha256()
                with self.source.open('rb') as stream:
                    for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                        digest.update(block)
                if fingerprint(self.source) != self.source_stat:
                    raise ValueError('GBZ changed while hashing')
                self.metadata = dict(metric=METRIC, source=self.source_stat,
                    source_sha256=digest.hexdigest(),
                    counting='distinct GBWT path IDs visiting either node orientation; repeated visits and reverse-complement copies counted once; includes reference, haplotype and generic paths')
                self.db.execute('CREATE TABLE gbz_metadata (value TEXT NOT NULL)')
                self.db.execute('CREATE TABLE gbz_counts (node_id INTEGER PRIMARY KEY, count INTEGER NOT NULL, sequence TEXT NOT NULL)')
                self.db.execute('INSERT INTO gbz_metadata VALUES (?)', (json.dumps(self.metadata),))
                self.db.commit()
        except BaseException:
            self.db.close()
            raise

    def _start(self):
        if self.process is not None:
            return
        self.process = subprocess.Popen([self.helper, str(self.source)], stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, text=True, bufsize=1)
        line = self.process.stdout.readline()
        if not line:
            raise ValueError('GBZ query helper failed to load the graph; see stderr')
        header = json.loads(line)
        if header.get('protocol_version') != 1 or header.get('metric') != METRIC:
            raise ValueError('Unsupported GBZ query helper protocol/metric')
        self.metadata['index'] = header
        self.db.execute('UPDATE gbz_metadata SET value=?', (json.dumps(self.metadata),))
        self.db.commit()

    def get_counts(self, nodes, sequences=None):
        if fingerprint(self.source) != self.source_stat:
            raise ValueError('GBZ changed during occurrence queries')
        result = {}
        for node in sorted(set(map(int, nodes))):
            row = self.db.execute('SELECT count, sequence FROM gbz_counts WHERE node_id=?', (node,)).fetchone()
            if row is None:
                self._start()
                self.process.stdin.write(json.dumps({'node_id': node}) + '\n')
                self.process.stdin.flush()
                line = self.process.stdout.readline()
                if not line:
                    raise ValueError(f'GBZ query failed for node {node}; see stderr')
                response = json.loads(line)
                if response.get('node_id') != node:
                    raise ValueError('GBZ helper returned the wrong node ID')
                row = response['count'], response['sequence'].upper()
                if not isinstance(row[0], int) or not 0 <= row[0] <= 2147483647:
                    raise ValueError(f'Invalid GBZ path count at node {node}')
                self.db.execute('INSERT INTO gbz_counts VALUES (?, ?, ?)', (node, *row))
            count, sequence = row
            if not isinstance(count, int) or not 0 <= count <= 2147483647:
                raise ValueError(f'Invalid cached GBZ path count at node {node}')
            if sequences is not None and sequence != sequences[node].upper():
                raise ValueError(f'GBZ/graph sequence mismatch at node {node}; node IDs belong to different graphs')
            result[node] = count
        self.db.commit()
        return result

    def close(self):
        if self.process is not None:
            self.process.stdin.close()
            self.process.stdout.close()
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.process.terminate()
                self.process.wait()
        self.db.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

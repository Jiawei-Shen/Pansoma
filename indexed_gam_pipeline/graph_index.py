"""Read-only combined GBZ sequence / distinct logical path index."""
import json
import sqlite3
import time
from pathlib import Path

SCHEMA = 'gbz-all-nodes-v1'
METRIC = 'distinct-gbwt-paths-v1'


class GraphIndex:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.db = sqlite3.connect(self.path.as_uri() + '?mode=ro', uri=True)
        try:
            rows = self.db.execute('SELECT value FROM graph_metadata').fetchall()
            if len(rows) != 1:
                raise ValueError('Missing/ambiguous graph index metadata')
            self.metadata = json.loads(rows[0][0])
            if (self.metadata.get('schema') != SCHEMA or self.metadata.get('metric') != METRIC
                    or self.metadata.get('status') != 'complete'):
                raise ValueError('Incomplete or incompatible graph index')
            self.db.execute('SELECT node_id,seq,distinct_path_count FROM nodes LIMIT 0')
            self.db.execute('BEGIN')
        except BaseException:
            self.db.close()
            raise
        self.performance = dict(queries=0, nodes_read=0, query_seconds=0.)

    def get_nodes(self, nodes):
        start = time.perf_counter()
        ids = sorted(set(map(int, nodes)))
        result = {}
        for offset in range(0, len(ids), 900):
            chunk = ids[offset:offset+900]
            query = 'SELECT node_id,seq,distinct_path_count FROM nodes WHERE node_id IN (' + ','.join('?' for _ in chunk) + ')'
            for node, seq, count in self.db.execute(query, chunk):
                if not seq or not isinstance(count, int) or not 0 <= count <= 2147483647:
                    raise ValueError(f'Invalid graph index record: {node}')
                result[node] = dict(sequence=seq, distinct_path_count=count)
            self.performance['queries'] += 1
        missing = set(ids) - result.keys()
        if missing:
            raise ValueError(f'Nodes missing from GBZ graph index: {sorted(missing)[:10]}')
        self.performance['nodes_read'] += len(result)
        self.performance['query_seconds'] += time.perf_counter() - start
        return result

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.db.close()

"""Build once from GBZ; atomically publish a complete unified SQLite index."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import tempfile
import time
from indexed_gam_pipeline.graph_index import GraphIndex, SCHEMA, METRIC


def fingerprint(path):
    s = path.stat()
    return dict(path=str(path), size=s.st_size, mtime_ns=s.st_mtime_ns)


def build_index(gbz, output, builder):
    start = time.perf_counter()
    gbz, output = Path(gbz).resolve(), Path(output).absolute()
    if output.exists():
        raise FileExistsError(output)
    source = fingerprint(gbz)
    digest = hashlib.sha256()
    with gbz.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(chunk)
    hash_seconds = time.perf_counter() - start
    metadata = dict(schema=SCHEMA, metric=METRIC, source=dict(source, sha256=digest.hexdigest()),
                    count_definition='distinct logical GBWT paths, either orientation, including reference paths; revisits count once')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=output.name+'.building-', dir=output.parent) as temp:
        temp = Path(temp)
        meta = temp/'metadata.json'
        meta.write_text(json.dumps(metadata))
        dbpath = temp/'index.sqlite'
        subprocess.run([str(Path(builder).resolve()), str(gbz), str(dbpath), str(meta)], check=True)
        if fingerprint(gbz) != source:
            raise ValueError('GBZ changed during index build')
        t = time.perf_counter()
        with GraphIndex(dbpath) as index:
            if index.db.execute('SELECT count(*) FROM nodes').fetchone()[0] != index.metadata['nodes']:
                raise ValueError('Incomplete graph index node table')
            if index.db.execute('PRAGMA quick_check').fetchall() != [('ok',)]:
                raise ValueError('SQLite integrity check failed')
            result = index.metadata
        result['timing'].update(source_hash_seconds=hash_seconds, validation_seconds=time.perf_counter()-t,
                               total_build_seconds=time.perf_counter()-start)
        with sqlite3.connect(dbpath) as db:
            db.execute('UPDATE graph_metadata SET value=?', (json.dumps(result),))
        # Same-filesystem link publishes atomically and refuses overwrite, even in a race.
        os.link(dbpath, output)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gbz', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--builder', required=True)
    a = p.parse_args()
    print(json.dumps(build_index(a.gbz, a.output, a.builder), indent=2))


if __name__ == '__main__':
    main()

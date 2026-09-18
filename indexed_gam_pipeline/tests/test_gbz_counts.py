"""GBZ occurrence regressions. Set GBZ_QUERY and GBZ_TOOL for native integration."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from indexed_gam_pipeline.gbz_counts import GBZCounts


class CacheTests(unittest.TestCase):
    def test_reuse_sequence_validation_and_source_identity(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            source = root / 'graph.gbz'
            source.write_bytes(b'fixture')
            cache = root / 'counts.sqlite'
            with GBZCounts(source, root / 'unused', cache) as lookup:
                lookup.db.execute('INSERT INTO gbz_counts VALUES (1, 86, "AC")')
                lookup.db.commit()
            with GBZCounts(source, root / 'unused', cache) as lookup:
                with patch.object(lookup, '_start', side_effect=AssertionError('cache miss')):
                    self.assertEqual(lookup.get_counts([1, 1], {1:'AC'}), {1:86})
                with self.assertRaisesRegex(ValueError, 'sequence mismatch'):
                    lookup.get_counts([1], {1:'GT'})
                source.write_bytes(b'changed graph')
                with self.assertRaisesRegex(ValueError, 'changed'):
                    lookup.get_counts([1])
            with self.assertRaisesRegex(ValueError, 'source/metric mismatch'):
                GBZCounts(source, root / 'unused', cache)

    def test_v4_build_contract(self):
        import argparse
        import io
        from contextlib import redirect_stdout
        import numpy as np
        from test_pipeline import fixture
        from indexed_gam_pipeline.run import build
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            gam, _ = fixture(d)
            nodes = [10,20,30,1000]
            (root/'nodes.txt').write_text('\n'.join(map(str,nodes)))
            (root/'nodes.json').write_text(json.dumps([dict(node_id=n,sequence='AAAAAA') for n in nodes]))
            gbz = root/'graph.gbz'
            gbz.write_bytes(b'cache fixture')
            cache = root/'counts.sqlite'
            with GBZCounts(gbz, root/'unused', cache) as lookup:
                lookup.db.executemany('INSERT INTO gbz_counts VALUES (?, ?, ?)', [(n,86,'AAAAAA') for n in nodes])
                lookup.db.commit()
            args = argparse.Namespace(command='build',format='candidate-v4',gam=str(gam),index=None,
                nodes=str(root/'nodes.txt'),node_json=str(root/'nodes.json'),node_sqlite=None,gfa=None,
                gbz=str(gbz),gbz_query=str(root/'unused'),occurrence_cache=str(cache),
                output=str(root/'out'),batch_nodes=2,max_node_span=100,max_batch_segments=100,
                shard_size=2,min_mapq=10,min_af=.05,min_variants=1,min_allele_bq=10.,
                variant_type='all',max_indel_len=50,rows=200,width=100,debug_rows=True)
            with redirect_stdout(io.StringIO()):
                build(args)
            manifest = json.loads((root/'out/manifest.json').read_text())
            self.assertEqual(manifest['schema_version'],4)
            self.assertEqual(manifest['channels'][-1],'node_distinct_gbwt_path_count')
            self.assertNotIn('walk_count_lookup',manifest)
            for line in (root/'out/variant_summary.ndjson').read_text().splitlines():
                meta = json.loads(line)
                self.assertEqual(meta['tensor_format_version'],'indexed-gam-candidate-v4')
                x = np.load(root/f"out/shard_{meta['shard_index']:05d}_data.npy")[meta['index_within_shard']]
                for ri,row in enumerate(meta['rows']):
                    for ci,col in enumerate(row['columns']):
                        self.assertEqual(x[6,ri,ci],86 if col else 0)

    def test_reject_old_gfa_cache(self):
        import sqlite3
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            source = root / 'graph.gbz'
            source.write_bytes(b'fixture')
            cache = root / 'counts.sqlite'
            with sqlite3.connect(cache) as db:
                db.execute('CREATE TABLE counts (node_id INTEGER, count INTEGER)')
            with self.assertRaisesRegex(ValueError, 'GFA W caches cannot be reused'):
                GBZCounts(source, root / 'unused', cache)


@unittest.skipUnless(os.environ.get('GBZ_QUERY') and os.environ.get('GBZ_TOOL'), 'native GBZ tools not configured')
class NativeTests(unittest.TestCase):
    def test_repeated_visits_and_both_orientations_count_paths_once(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            gfa = root / 'graph.gfa'
            gbz = root / 'graph.gbz'
            gfa.write_text('H\tVN:Z:1.1\nS\t1\tA\nS\t2\tC\n'
                'L\t1\t+\t2\t+\t0M\nL\t2\t+\t1\t+\t0M\nL\t1\t+\t1\t-\t0M\n'
                'W\ts\t1\tchr1\t0\t4\t>1>2>1<1\nW\tt\t1\tchr1\t0\t2\t<2<1\n')
            subprocess.run([os.environ['GBZ_TOOL'], 'convert', str(gfa), str(gbz)], check=True, capture_output=True)
            with GBZCounts(gbz, os.environ['GBZ_QUERY'], root / 'counts.sqlite') as lookup:
                self.assertEqual(lookup.get_counts([1,2], {1:'A',2:'C'}), {1:2,2:2})
                self.assertEqual(lookup.metadata['index']['indexed_paths'], 2)
                with self.assertRaisesRegex(ValueError, 'query failed'):
                    lookup.get_counts([99])


if __name__ == '__main__':
    unittest.main()

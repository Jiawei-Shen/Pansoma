import json
import os
from pathlib import Path
import sqlite3
import subprocess
import tempfile
import unittest
from indexed_gam_pipeline.graph_index import GraphIndex, SCHEMA, METRIC
from indexed_gam_pipeline.build_graph_index import build_index


def fixture_index(path, nodes):
    with sqlite3.connect(path) as db:
        db.execute('CREATE TABLE graph_metadata(value TEXT)')
        db.execute('INSERT INTO graph_metadata VALUES(?)', (json.dumps(dict(schema=SCHEMA, metric=METRIC, status='complete')),))
        db.execute('CREATE TABLE nodes(node_id INTEGER PRIMARY KEY, seq TEXT, distinct_path_count INTEGER)')
        db.executemany('INSERT INTO nodes VALUES(?,?,?)', nodes)


class IndexTests(unittest.TestCase):
    def test_batch_read_missing_zero_and_readonly(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'index.sqlite'
            fixture_index(p, [(i,'AC',i%91) for i in range(1,2001)])
            with GraphIndex(p) as idx:
                self.assertEqual(len(idx.get_nodes(range(1,2001))), 2000)
                self.assertEqual(idx.get_nodes([91,91])[91]['distinct_path_count'], 0)
                self.assertEqual(idx.get_nodes([]), {})
                with self.assertRaisesRegex(ValueError, 'missing'):
                    idx.get_nodes([3000])
                with self.assertRaises(sqlite3.OperationalError):
                    idx.db.execute('DELETE FROM nodes')
            with sqlite3.connect(p) as db:
                db.execute('UPDATE graph_metadata SET value=?', ('{}',))
            with self.assertRaisesRegex(ValueError, 'Incomplete'):
                GraphIndex(p)

    def test_tensor_equivalence(self):
        import argparse
        import io
        from contextlib import redirect_stdout
        import numpy as np
        from test_pipeline import fixture
        from indexed_gam_pipeline.run import build
        from indexed_gam_pipeline.gbz_counts import GBZCounts
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            gam, _ = fixture(d)
            nodes = [10,20,30,1000]
            (root/'nodes.txt').write_text('\n'.join(map(str,nodes)))
            (root/'nodes.json').write_text(json.dumps([dict(node_id=n,sequence='AAAAAA') for n in nodes]))
            gbz = root/'graph.gbz'; gbz.write_bytes(b'fixture')
            with GBZCounts(gbz, root/'unused', root/'counts.sqlite') as lookup:
                lookup.db.executemany('INSERT INTO gbz_counts VALUES(?,?,?)', [(n,86,'AAAAAA') for n in nodes])
                lookup.db.commit()
            args = argparse.Namespace(command='build',format='candidate-v4',gam=str(gam),index=None,
                nodes=str(root/'nodes.txt'),node_json=str(root/'nodes.json'),node_sqlite=None,gfa=None,
                gbz=str(gbz),gbz_query=str(root/'unused'),occurrence_cache=str(root/'counts.sqlite'),
                output=str(root/'old'),batch_nodes=2,max_node_span=100,max_batch_segments=100,
                shard_size=2,min_mapq=10,min_af=.05,min_variants=1,min_allele_bq=10.,
                variant_type='all',max_indel_len=50,rows=200,width=101,debug_rows=True, early_alt_filter=True)
            with redirect_stdout(io.StringIO()): build(args)
            fixture_index(root/'unified.sqlite', [(n,'AAAAAA',86) for n in nodes])
            args.node_json=args.gbz=args.gbz_query=args.occurrence_cache=None
            args.graph_index=str(root/'unified.sqlite');args.output=str(root/'new')
            with redirect_stdout(io.StringIO()): build(args)
            oldfiles=list((root/'old').glob('*.npy'))
            self.assertTrue(oldfiles)
            for old in oldfiles:
                np.testing.assert_array_equal(np.load(old), np.load(root/'new'/old.name))
            self.assertEqual((root/'old/variant_summary.ndjson').read_text(), (root/'new/variant_summary.ndjson').read_text())
            manifest=json.loads((root/'new/manifest.json').read_text())
            self.assertIn('graph_index_seconds',manifest['timing'])
            from indexed_gam_pipeline.validate_examples import validate
            self.assertTrue(validate(root/'new', str(gam), None,
                                     graph_index=str(root/'unified.sqlite'))['passed'])
            args.node_json=str(root/'nodes.json')
            with self.assertRaisesRegex(ValueError, 'replaces legacy'):
                build(args)


@unittest.skipUnless(all(os.environ.get(k) for k in ('GBZ_TOOL','GBZ_QUERY','GBZ_GRAPH_INDEX')), 'native tools not configured')
class NativeTests(unittest.TestCase):
    def test_counts_match_query_and_atomic_publication(self):
        from indexed_gam_pipeline.gbz_counts import GBZCounts
        for second in (2,2001):
            with self.subTest(second=second), tempfile.TemporaryDirectory() as d:
                root=Path(d);gfa=root/'graph.gfa';gbz=root/'graph.gbz';output=root/'index.sqlite'
                gfa.write_text(f'H\tVN:Z:1.1\nS\t1\tA\nS\t{second}\tC\n'
                    f'L\t1\t+\t{second}\t+\t0M\nL\t{second}\t+\t1\t+\t0M\nL\t1\t+\t1\t-\t0M\n'
                    f'W\ts\t1\tchr1\t0\t4\t>1>{second}>1<1\nW\tt\t1\tchr1\t0\t2\t<{second}<1\n')
                subprocess.run([os.environ['GBZ_TOOL'],'convert',str(gfa),str(gbz)],check=True,capture_output=True)
                result=build_index(gbz,output,os.environ['GBZ_GRAPH_INDEX'])
                self.assertEqual(result['nodes'],2)
                with GraphIndex(output) as idx, GBZCounts(gbz,os.environ['GBZ_QUERY'],root/'old.sqlite') as old:
                    records=idx.get_nodes([1,second]);sequences={n:r['sequence'] for n,r in records.items()}
                    self.assertEqual(sequences,{1:'A',second:'C'})
                    counts={n:r['distinct_path_count'] for n,r in records.items()}
                    self.assertEqual(counts,old.get_counts([1,second],sequences))
                    self.assertEqual(counts,{1:2,second:2})
                with self.assertRaises(FileExistsError):build_index(gbz,output,os.environ['GBZ_GRAPH_INDEX'])
                gbz.write_bytes(b'invalid')
                with self.assertRaises(subprocess.CalledProcessError):build_index(gbz,root/'bad.sqlite',os.environ['GBZ_GRAPH_INDEX'])
                self.assertFalse((root/'bad.sqlite').exists())

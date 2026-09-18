"""Single-pass W cache and candidate-v3 channel regressions."""
import argparse
from contextlib import redirect_stdout
import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from indexed_gam_pipeline.walk_counts import build_walk_counts, WalkCounts, walk_nodes
from indexed_gam_pipeline.candidates import Candidate, make_tensor
from indexed_gam_pipeline.run import build
from test_candidates import read, eligible
from test_pipeline import fixture

GFA = (b'H\tVN:Z:1.1\nS\t1\tA\nS\t2\tT\nP\tp\t1+,2+\t*\n'
       b'W\ts\t0\tchr1\t0\t3\t>1>2>1<1\n'
       b'W\ts\t0\tchr1\t0\t3\t>1>2>1<1\tXX:Z:duplicate\n'
       b'W\ts\t1\tchr1\t0\t3\t<1>3\n'
       b'W\ts\t0\tchr1\t3\t6\t>1>4\n'
       b'W\tt\t1\tchr2\t0\t1\t>100000000001\n')


class WalkCountTest(unittest.TestCase):
    def cache(self, folder, data=GFA, compressed=False):
        source=folder/('graph.gfa.gz' if compressed else 'graph.gfa')
        if compressed:
            with gzip.open(source,'wb') as stream:
                stream.write(data)
        else:
            source.write_bytes(data)
        target=folder/'walks.sqlite'
        with redirect_stdout(io.StringIO()):
            report=build_walk_counts(source,target)
        return source,target,report

    def test_distinct_records_visits_strands_sparse_ids_and_reuse(self):
        with tempfile.TemporaryDirectory() as d:
            source,cache,report=self.cache(Path(d))
            self.assertEqual((report['w_records'],report['distinct_w_records'],report['duplicate_w_records']),(5,4,1))
            # Lookup uses only cached counts; GFA content is never reopened.
            with patch('indexed_gam_pipeline.walk_counts.open',side_effect=AssertionError('GFA rescanned'),create=True):
                with WalkCounts(cache) as lookup:
                    self.assertEqual(lookup.get_counts([1,2,3,4,5,100000000001]),
                                     {1:3,2:1,3:1,4:1,5:0,100000000001:1})
                    self.assertEqual(len(lookup.get_counts(range(1,1900))),1899)
            with self.assertRaisesRegex(ValueError,'already exists'):
                build_walk_counts(source,cache)
            source.write_bytes(GFA+b'S\t999\tT\n')
            with self.assertRaisesRegex(ValueError,'stale'):
                WalkCounts(cache)

    def test_compressed_input_and_explicit_invalid_input_failure(self):
        with tempfile.TemporaryDirectory() as d:
            _,cache,report=self.cache(Path(d),compressed=True)
            with WalkCounts(cache) as lookup:
                self.assertEqual(lookup.get_counts([1]),{1:3})
        for walk in (b'>abc',b'>>1',b'>0',b'>-1',b'>1?2',b'*',b'1>2',b'>9223372036854775808',b'>18446744073709551617'):
            with self.assertRaises(ValueError):
                walk_nodes(walk)
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);source=folder/'no_walks.gfa';source.write_text('S\t1\tA\n')
            with self.assertRaisesRegex(ValueError,'no W records'):
                build_walk_counts(source,folder/'failed.sqlite')
            self.assertFalse(list(folder.glob('*.sqlite*')))

    def test_seventh_channel_exact_counts_gaps_padding_and_grouping(self):
        seq={1:'ACGT',2:'TT'}
        alt=read([(2,0,False,[(2,2,'')]),(1,0,False,[(2,2,''),(0,2,'TA'),(2,2,'')])],seq)
        ref=read([(1,0,False,[(4,4,'')])],seq)
        c=alt.observations[0].candidate
        e=eligible(c,[alt,ref])
        six,old=make_tensor(c,e,rows=3,width=10,debug=True)
        seven,new=make_tensor(c,e,rows=3,width=10,debug=True,node_walk_counts={1:40000,2:7})
        self.assertEqual(seven.shape,(7,3,10))
        self.assertEqual(seven.dtype,np.int32)
        np.testing.assert_array_equal(seven[:6],six)
        self.assertEqual(new['row_groups'],old['row_groups'])
        self.assertEqual(new['coverage'],old['coverage'])
        for ri,row in enumerate(new['rows']):
            for ci,col in enumerate(row['columns']):
                self.assertEqual(seven[6,ri,ci],{1:40000,2:7}[col['node_id']] if col else 0)
        self.assertFalse(seven[:,2].any())
        with self.assertRaisesRegex(ValueError,'int32'):
            make_tensor(c,e,node_walk_counts={1:2**31,2:1})
        # Deletions retain the node count even though the read has a gap.
        deletion=read([(1,0,False,[(1,1,''),(1,0,''),(2,2,'')])],seq)
        dc=deletion.observations[0].candidate
        x,m=make_tensor(dc,eligible(dc,[deletion]),node_walk_counts={1:0})
        pos=m['candidate_columns'][0]
        self.assertEqual(x[0,0,pos],6)
        self.assertEqual(x[6,0,pos],0)

    def test_v3_indexed_shard_contract_and_cached_channel(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);gam,_=fixture(d)
            _,cache,_=self.cache(folder,b'W\ts\t0\tchr1\t0\t24\t>10>20>30>1000>10\nW\tt\t0\tchr1\t0\t6\t<10\n')
            (folder/'nodes.txt').write_text('10\n20\n30\n1000\n')
            (folder/'nodes.json').write_text(json.dumps([dict(node_id=n,sequence='AAAAAA') for n in (10,20,30,1000)]))
            args=argparse.Namespace(command='build',format='candidate-v3',gam=str(gam),index=None,
                nodes=str(folder/'nodes.txt'),node_json=str(folder/'nodes.json'),node_sqlite=None,gfa=None,
                walk_counts=str(cache),output=str(folder/'out'),batch_nodes=2,max_node_span=100,
                max_batch_segments=100,shard_size=2,min_mapq=10,min_af=.05,min_variants=1,
                min_allele_bq=10.,variant_type='all',max_indel_len=50,rows=200,width=100,debug_rows=True)
            with redirect_stdout(io.StringIO()):
                build(args)
            manifest=json.loads((folder/'out/manifest.json').read_text())
            self.assertEqual(manifest['schema_version'],3)
            self.assertEqual(manifest['shape'],[7,200,100])
            self.assertEqual(manifest['dtype'],'int32')
            summary=[json.loads(s) for s in (folder/'out/variant_summary.ndjson').read_text().splitlines()]
            self.assertEqual(len(summary),4)
            for m in summary:
                self.assertEqual(m['tensor_format_version'],'indexed-gam-candidate-v3')
                x=np.load(folder/f"out/shard_{m['shard_index']:05d}_data.npy")[m['index_within_shard']]
                for ri,row in enumerate(m['rows']):
                    for ci,col in enumerate(row['columns']):
                        self.assertEqual(x[6,ri,ci],(2 if col['node_id']==10 else 1) if col else 0)
            import sqlite3
            from indexed_gam_pipeline.validate_examples import validate
            graph=folder/'graph.sqlite'
            with sqlite3.connect(graph) as connection:
                connection.execute('CREATE TABLE nodes (node_id TEXT PRIMARY KEY, seq TEXT)')
                connection.executemany('INSERT INTO nodes VALUES (?, ?)',[(str(n),'AAAAAA') for n in (10,20,30,1000)])
            audit=validate(folder/'out',str(gam),None,str(graph),str(cache))
            self.assertTrue(audit['walk_counts_checked'])
            shard=folder/'out/shard_00000_data.npy'
            x=np.load(shard);x[0,6,0,49]+=1;np.save(shard,x)
            with self.assertRaisesRegex(ValueError,'walk count'):
                validate(folder/'out',str(gam),None,str(graph),str(cache))


if __name__=='__main__':
    unittest.main()

import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
import numpy as np
from indexed_gam_pipeline.full_run import graph_index, validate_shards


class FullRunTests(unittest.TestCase):
    def test_complete_graph_index_includes_alternate_nodes(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);gfa=p/'graph.gfa';db=p/'nodes.sqlite'
            gfa.write_text('H\tVN:Z:1.1\nS\t1\tac\nL\t1\t+\t9\t+\t0M\nS\t9\tG\tLN:i:1\n')
            self.assertEqual(graph_index(gfa,db),2)
            with sqlite3.connect(db) as c:
                self.assertEqual(c.execute('SELECT node_id,seq FROM nodes ORDER BY node_id').fetchall(),[(1,'AC'),(9,'G')])
            with self.assertRaisesRegex(ValueError,'already exists'): graph_index(gfa,db)

    def test_shard_validation_checks_partial_last_shard_and_metadata(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            np.save(p/'shard_00000_data.npy',np.zeros((2,7,200,100),dtype=np.int32))
            np.save(p/'shard_00001_data.npy',np.zeros((1,7,200,100),dtype=np.int32))
            version='indexed-gam-candidate-v4'
            (p/'manifest.json').write_text(json.dumps(dict(status='complete',tensor_format_version=version,shards=2,tensors=3)))
            rows=[dict(shard_index=s,index_within_shard=i,coverage=3,alt_count=3,ref_count=0,other_count=0,af=1.,selected_alignments=3,tensor_format_version=version,event_type='SNP') for s,n in enumerate((2,1)) for i in range(n)]
            summary=p/'variant_summary.ndjson';summary.write_text(''.join(json.dumps(r)+'\n' for r in rows))
            self.assertEqual(validate_shards(p,2)['tensors'],3)
            rows[-1]['index_within_shard']=1
            summary.write_text(''.join(json.dumps(r)+'\n' for r in rows))
            with self.assertRaisesRegex(ValueError,'ordering'):validate_shards(p,2)

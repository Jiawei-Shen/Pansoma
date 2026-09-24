import argparse
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from indexed_gam_pipeline.candidates import decode_alignment, make_tensor, overlap
from indexed_gam_pipeline.run import build
from test_candidates import alignment
from test_pipeline import fixture


class TargetObservationsTest(unittest.TestCase):
    def test_context_paths_unsupported_and_target_tensors_unchanged(self):
        seq={1:'ACGTAC',2:'ACGTAC',3:'A'*60}
        for reverse in (False,True):
            a=alignment([(1,0,reverse,[(1,1,'T'),(1,1,''),(0,2,'GG'),(2,0,''),(2,2,'')]),
                         (2,0,reverse,[(1,1,'T'),(1,1,''),(0,2,'GG'),(2,0,''),(2,2,'')]),
                         (1,0,reverse,[(6,6,'')]),
                         (3,0,reverse,[(0,51,'C'*51),(60,60,'')])],seq)
            old,rejected=decode_alignment(a,seq)
            for targets in (set(),{1},{2},{1,2,3}):
                new,new_rejected=decode_alignment(a,seq,target_nodes=targets)
                self.assertEqual(old.columns,new.columns)
                self.assertEqual(old.visits,new.visits)
                self.assertEqual(old.digest,new.digest)
                self.assertEqual(rejected,new_rejected)
                self.assertEqual(new.observations,[o for o in old.observations if o.candidate.node in targets])
                for o in new.observations:
                    c=o.candidate
                    oe=[(old,*overlap(old,c,10))];ne=[(new,*overlap(new,c,10))]
                    x,m=make_tensor(c,oe,rows=2,width=101,debug=True)
                    y,n=make_tensor(c,ne,rows=2,width=101,debug=True)
                    np.testing.assert_array_equal(x,y);self.assertEqual(m,n)

    def test_multibatch_build_matches_unrestricted_observations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);gam,_=fixture(root)
            (root/'nodes.txt').write_text('10\n20\n30\n1000\n')
            (root/'nodes.json').write_text(json.dumps([dict(node_id=n,sequence='AAAAAA') for n in (10,20,30,1000)]))
            args=argparse.Namespace(command='build',format='candidate-v2',gam=str(gam),index=None,
                nodes=str(root/'nodes.txt'),node_json=str(root/'nodes.json'),node_sqlite=None,gfa=None,
                batch_nodes=1,max_node_span=10000,max_batch_segments=100,shard_size=2,
                min_mapq=10,min_af=.05,min_variants=1,min_allele_bq=10.,variant_type='all',
                max_indel_len=50,rows=4,width=101,debug_rows=True,workers=1,early_alt_filter=True,
                decode_mode='full',gam_cache_mb=1)
            def unrestricted(*a,**kw):
                kw['target_nodes']=None
                return decode_alignment(*a,**kw)
            args.output=str(root/'baseline')
            with patch('indexed_gam_pipeline.build_v2.decode_alignment',unrestricted), redirect_stdout(io.StringIO()):
                build(args)
            args.output=str(root/'optimized')
            with redirect_stdout(io.StringIO()):build(args)
            for name in ['variant_summary.ndjson','filtered_candidates.ndjson','unsupported_events.ndjson']:
                self.assertEqual((root/'baseline'/name).read_bytes(),(root/'optimized'/name).read_bytes(),name)
            for p in (root/'baseline').glob('shard_*_data.npy'):
                self.assertEqual(p.read_bytes(),(root/'optimized'/p.name).read_bytes(),p.name)

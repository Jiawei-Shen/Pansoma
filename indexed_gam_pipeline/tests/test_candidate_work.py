import argparse
from collections import Counter
import unittest
import numpy as np

from indexed_gam_pipeline.candidate_work import NodeIndexCache,alt_support_bounds,candidate_results
from indexed_gam_pipeline.candidates import Candidate,overlap,decode_alignment
from test_candidates import alignment


class CandidateWorkTests(unittest.TestCase):
    def fixture(self):
        seq={1:'ACGT',2:'TT'}
        a=alignment([(1,0,False,[(1,1,''),(1,1,'T'),(2,2,'')]),
                     (2,0,False,[(2,2,'')]),(1,0,False,[(1,1,''),(1,1,'T'),(2,2,'')])],seq)
        b=alignment([(1,0,True,[(2,2,''),(1,1,'A'),(1,1,'')])],seq,name='reverse')
        b.quality=bytes([2]*4)
        c=alignment([(1,0,False,[(1,1,''),(2,0,''),(1,1,'')])],seq,name='del')
        d=alignment([(1,0,False,[(2,2,''),(0,2,'TT'),(2,2,'')])],seq,name='ins')
        return [decode_alignment(x,seq)[0] for x in (a,a,b,c,d)]

    def test_upper_bound_deduplicates_visits_not_records(self):
        reads=self.fixture();candidates={o.candidate for r in reads for o in r.observations}
        counts=alt_support_bounds(reads,candidates,10)
        c=Candidate(1,1,'C','T','SNP')
        self.assertEqual(counts[c],2) # duplicate records count twice; repeat visits once
        self.assertEqual(alt_support_bounds(reads,candidates,0)[c],3)
        for candidate in candidates:
            actual=sum(overlap(r,candidate,10)[0]=='alt' for r in reads if overlap(r,candidate,10))
            self.assertLessEqual(actual,counts[candidate])

    def test_index_matches_reference_for_repeats_reverse_indels_and_quality(self):
        reads=self.fixture();cache=NodeIndexCache();idx=cache.get(1,reads)
        for r in reads:
            for candidate in {o.candidate for x in reads for o in x.observations}:
                for bq in (0,10,50):
                    self.assertEqual(overlap(r,candidate,bq),overlap(r,candidate,bq,idx[id(r)]))

    def test_fifo_count_and_estimated_memory_limits(self):
        reads=self.fixture();cache=NodeIndexCache(2,100000)
        cache.get(1,reads);cache.get(2,reads);cache.get(1,reads);cache.get(3,reads)
        self.assertEqual(list(cache.entries),[2,3]) # reads do not turn FIFO into LRU
        self.assertEqual(cache.stats['peak_nodes'],2)
        self.assertEqual(cache.stats['evictions'],1)
        self.assertIsNone(NodeIndexCache(1000,1).get(1,reads))
        self.assertIsNone(NodeIndexCache(0).get(1,reads))
        with self.assertRaises(ValueError):NodeIndexCache(1001)
        cache=NodeIndexCache(1000)
        for n in range(1002):cache.get(n,[])
        self.assertEqual(len(cache.entries),1000)
        self.assertEqual(next(iter(cache.entries)),2)
        self.assertLessEqual(cache.bytes,cache.max_bytes)

    def test_parallel_candidates_preserve_order_tensors_and_metadata(self):
        reads=self.fixture();candidates=sorted({o.candidate for r in reads for o in r.observations})
        args=argparse.Namespace(workers=1,node_index_cache_nodes=0,node_index_cache_mb=1,
            min_allele_bq=10,min_variants=1,min_af=.05,variant_type='all',rows=200,width=101,debug_rows=True)
        tasks=[(1,[c]) for c in candidates]
        def collect():
            with candidate_results(tasks,{1:reads},args,{1:90,2:3}) as results:
                return [result for chunk,_,_ in results for result in chunk]
        baseline=collect()
        args.node_index_cache_nodes=1000
        for workers in (1,2):
            args.workers=workers;optimized=collect()
            self.assertEqual(len(optimized),len(baseline))
            for (x,m),(y,n) in zip(baseline,optimized):
                np.testing.assert_array_equal(x,y);self.assertEqual(m,n)

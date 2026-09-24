"""Window decoder must preserve full decoder tensors and boundary semantics."""
import random
import argparse
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import test_candidates as legacy
from indexed_gam_pipeline.candidates import decode_alignment, Candidate, anchor_window, overlap, make_tensor
from indexed_gam_pipeline.edit_columns import select_decode_mode
from indexed_gam_pipeline.tests import test_pipeline


class ExistingWindowSemantics(legacy.CandidatesTest):
    def setUp(self):
        super().setUp()
        def window(*args, **kwargs):
            return decode_alignment(*args, **kwargs, mode='window')
        self.patch = patch.object(legacy, 'decode_alignment', window)
        self.patch.start()
        self.addCleanup(self.patch.stop)


class EditColumnsTest(unittest.TestCase):
    def test_auto_build_short_and_long_with_parallel_candidates(self):
        import pysam
        from indexed_gam_pipeline.gam_reader import build_index
        from indexed_gam_pipeline.run import build
        for length, selected in ((999, 'full'), (1000, 'window'), (10000, 'window')):
            with tempfile.TemporaryDirectory() as temp:
                folder = Path(temp)
                seq = {1:'A'*length}
                records = [legacy.alignment([(1,0,bool(i%2),[(length//2,length//2,''),
                    (1,1,'T'),(length-length//2-1,length-length//2-1,'')])], seq, name=str(i)) for i in range(12)]
                gam = folder/'reads.gam'
                with pysam.BGZFile(str(gam),'wb') as stream:
                    for a in records:
                        raw = a.SerializeToString()
                        stream.write(test_pipeline.vi(2)+test_pipeline.vi(3)+b'GAM'+test_pipeline.vi(len(raw))+raw)
                        stream.flush()
                build_index(gam, str(gam)+'.gai')
                (folder/'nodes.txt').write_text('1\n')
                (folder/'nodes.json').write_text(json.dumps([dict(node_id=1,sequence=seq[1])]))
                args = argparse.Namespace(command='build',format='candidate-v2',gam=str(gam),index=None,
                    nodes=str(folder/'nodes.txt'),node_json=str(folder/'nodes.json'),node_sqlite=None,gfa=None,
                    batch_nodes=512,max_node_span=10000,max_batch_segments=100,shard_size=2,
                    min_mapq=10,min_af=.05,min_variants=1,min_allele_bq=10.,variant_type='all',
                    max_indel_len=50,rows=4,width=101,debug_rows=True,workers=2,early_alt_filter=True)
                expected = None
                for mode in ('full','window','auto'):
                    args.decode_mode=mode; args.output=str(folder/mode)
                    with redirect_stdout(io.StringIO()): build(args)
                    out = folder/mode
                    manifest = json.loads((out/'manifest.json').read_text())
                    if mode == 'auto':
                        self.assertEqual(manifest['decode_policy']['selected'],selected)
                        self.assertEqual(len(manifest['decode_policy']['sample_lengths']),10)
                    actual={f.name:f.read_bytes() for f in out.glob('shard_*_data.npy')}
                    actual['summary']=(out/'variant_summary.ndjson').read_bytes()
                    if expected is None: expected=actual
                    else: self.assertEqual(expected,actual)

    def test_auto_threshold_and_first_ten(self):
        for lengths, expected in [([999]*11, 'full'), ([1000]*10, 'window'),
                                  ([1001]*10, 'window'), ([100]*10+[10000], 'full'),
                                  ([10000]+[100]*9, 'window'), ([100], 'full'),
                                  ([], 'full')]:
            records = [type('A', (), {'sequence': 'A'*n})() for n in lengths]
            with patch('indexed_gam_pipeline.gam_reader.scan_gam', return_value=iter_generator(records)):
                policy = select_decode_mode('unused')
            self.assertEqual(policy['selected'], expected)
            self.assertEqual(policy['sample_lengths'], lengths[:10])
        with patch('indexed_gam_pipeline.gam_reader.scan_gam', side_effect=AssertionError):
            self.assertEqual(select_decode_mode('unused', 'full')['selected'], 'full')
            self.assertEqual(select_decode_mode('unused', 'window')['selected'], 'window')

    def test_long_run_materialization_is_bounded(self):
        seq = {1: 'A'*100000}
        a = legacy.alignment([(1, 0, False, [(50000,50000,''), (1,1,'T'), (49999,49999,'')])], seq)
        full, _ = decode_alignment(a, seq)
        window, _ = decode_alignment(a, seq, mode='window')
        self.assertEqual(window.columns.materialized, 0)
        candidate = full.observations[0].candidate
        self.assertEqual(overlap(full, candidate, 10), overlap(window, candidate, 10))
        self.assertEqual(anchor_window(full, full.visits[0], candidate, 101),
                         anchor_window(window, window.visits[0], candidate, 101))
        self.assertLess(window.columns.materialized, 250)

    def test_long_background_insertion_is_not_expanded_for_support_or_overflow(self):
        seq = {1:'A'*200}
        for reverse in (False, True):
            a = legacy.alignment([(1,0,reverse,[(65,65,''),(0,10000,'G'*10000),(135,135,'')])],seq)
            full,_ = decode_alignment(a,seq)
            window,_ = decode_alignment(a,seq,mode='window')
            pos = 130 if reverse else 60
            candidate = Candidate(1,pos,'A'*20,'','DEL')
            self.assertEqual(overlap(full,candidate,10),overlap(window,candidate,10))
            self.assertEqual(window.columns.materialized,0)
            self.assertEqual(anchor_window(full,full.visits[0],candidate,101),
                             anchor_window(window,window.visits[0],candidate,101))
            self.assertLess(window.columns.materialized,300)

    def test_randomized_edits_paths_and_windows(self):
        rng = random.Random(1926)
        for case in range(100):
            seq = {n: ''.join(rng.choices('ACGT', k=180)) for n in (1, 2, 3)}
            specs = []
            for node in (1, 2, 1, 3):
                edits = []
                for _ in range(12):
                    op = rng.choice(('M', 'X', 'I', 'D', 'C'))
                    f = rng.randint(1, 10) if op != 'I' else 0
                    t = f if op in ('M', 'X') else 0 if op == 'D' else rng.randint(1, 12)
                    edits.append((f, t, '' if op in ('M', 'D') else ''.join(rng.choices('ACGT', k=t))))
                specs.append((node, rng.randint(0, 10), bool(rng.getrandbits(1)), edits))
            a = legacy.alignment(specs, seq, name=str(case))
            a.quality = b'' if case%10 == 0 else bytes(rng.randrange(0,41) for _ in a.sequence)
            full, rejected = decode_alignment(a, seq)
            window, other = decode_alignment(a, seq, mode='window')
            self.assertEqual(rejected, other)
            self.assertEqual(full.visits, window.visits)
            self.assertEqual(full.observations, window.observations)
            self.assertEqual(full.columns, window.columns[:])
            candidates = [o.candidate for o in full.observations[::max(1,len(full.observations)//10)]]
            candidates += [Candidate(1, p, '', 'TT', 'INS') for p in (0, 10, 50, 90)]
            candidates += [Candidate(1, p, 'A'*20, '', 'DEL') for p in (0, 10, 50)]
            for candidate in candidates:
                self.assertEqual(overlap(full, candidate, 10), overlap(window, candidate, 10))
                for visit in full.visits:
                    if visit.node != candidate.node:
                        continue
                    for width in (1, 2, 11, 101, 201):
                        self.assertEqual(anchor_window(full, visit, candidate, width),
                                         anchor_window(window, visit, candidate, width), (case, candidate, visit, width))
            for candidate in candidates[:2]:
                hit = overlap(full, candidate, 10)
                if hit:
                    x, m = make_tensor(candidate, [(full, *hit)], rows=2, debug=True)
                    y, n = make_tensor(candidate, [(window, *hit)], rows=2, debug=True)
                    np.testing.assert_array_equal(x, y)
                    self.assertEqual(m, n)


def iter_generator(records):
    yield from records


if __name__ == '__main__':
    unittest.main()

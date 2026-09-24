"""--early-af-filter: exact coverage and unchanged accepted outputs; --max-node-reads cap."""
import argparse
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

from indexed_gam_pipeline.candidate_work import exact_coverage, process_node
from indexed_gam_pipeline import candidate_work
from indexed_gam_pipeline.candidates import Candidate, decode_alignment, overlap
from indexed_gam_pipeline.run import build
from test_candidates import alignment
from test_candidate_work import CandidateWorkTests
from test_pipeline import vi
import gzip
import pysam


SEQ = {10: 'ACGTACGTAC', 20: 'TTGCAAGGCT'}


def write_gam(path, rows):
    """Same BGZF/GAI layout as test_pipeline.fixture, with caller-chosen records."""
    bins = {}
    with pysam.BGZFile(str(path), "wb") as stream:
        for a in rows:
            start = stream.tell()
            raw = a.SerializeToString()
            stream.write(vi(2) + vi(3) + b"GAM" + vi(len(raw)) + raw)
            stream.flush()
            end = stream.tell()
            ids = [m.position.node_id for m in a.path.mapping]
            shift = max(1, (min(ids) ^ max(ids)).bit_length())
            bins.setdefault((min(ids) >> shift) + ((1 << (64 - shift)) - 1), []).append((start, end))
    payload = b"GAI!" + vi(1) + vi(len(bins))
    for number, runs in bins.items():
        payload += vi(number) + vi(len(runs))
        for start, end in runs:
            payload += vi(start) + vi(end)
    payload += vi(0)
    with gzip.open(str(path) + ".gai", "wb") as stream:
        stream.write(payload)
    return path


def mixed_af_rows():
    """Node 10: SNV 4:A>C in 4/8 reads, SNV 2:G>T in 1/8, 1-bp DEL at 6 in 1/8, INS at 8 in 1/8.
    Node 20: SNV 3:C>A in 2/5, SNV 7:G>C in 1/5, one all-match record on the reverse strand."""
    rows = []
    def add(node, edits, reverse=False):
        rows.append(alignment([(node, 0, reverse, edits)], SEQ, name=f'r{len(rows)}'))
    for i in range(8):
        e = [(4, 4, '')]
        e += [(1, 1, 'C')] if i < 4 else [(1, 1, '')]
        e += [(5, 5, '')]
        if i == 4: e = [(2, 2, ''), (1, 1, 'T'), (7, 7, '')]
        if i == 5: e = [(4, 4, ''), (1, 1, 'C'), (1, 1, ''), (1, 0, ''), (3, 3, '')]
        if i == 6: e = [(8, 8, ''), (0, 2, 'GG'), (2, 2, '')]
        add(10, e)
    for i in range(5):
        e = [(3, 3, ''), (1, 1, 'A'), (6, 6, '')] if i < 2 else [(10, 10, '')]
        if i == 2: e = [(7, 7, ''), (1, 1, 'C'), (2, 2, '')]
        add(20, e, reverse=(i == 4))  # all-match record on the reverse strand
    return rows


class ExactCoverageTest(unittest.TestCase):
    def test_matches_overlap_hits_for_every_position_and_event(self):
        reads = CandidateWorkTests().fixture()
        seq = {1: 'ACGT', 2: 'TT'}
        by_node = {n: [r for r in reads if any(v.node == n for v in r.visits)] for n in seq}
        probes = {o.candidate for r in reads for o in r.observations}
        for node, s in seq.items():
            for pos in range(len(s) + 1):
                probes.add(Candidate(node, pos, '', 'G', 'INS'))
                for length in (1, 2, 3):
                    if pos + length <= len(s):
                        probes.add(Candidate(node, pos, s[pos:pos+length], '', 'DEL'))
                if pos < len(s):
                    probes.add(Candidate(node, pos, s[pos], 'N' if s[pos] != 'N' else 'A', 'SNP'))
        coverage = exact_coverage(probes, by_node)
        for c in probes:
            expected = sum(1 for r in by_node.get(c.node, []) if overlap(r, c, 10) is not None)
            self.assertEqual(coverage[c], expected, c)

    def test_uncovered_node_has_zero_coverage(self):
        self.assertEqual(exact_coverage([Candidate(99, 0, 'A', 'C', 'SNP')], {}), {Candidate(99, 0, 'A', 'C', 'SNP'): 0})


class MaxNodeReadsTest(unittest.TestCase):
    def test_cap_is_deterministic_and_shared_by_candidates(self):
        seq = {1: 'ACGT'}
        reads = [decode_alignment(alignment([(1, 0, False, [(1, 1, ''), (1, 1, 'T'), (2, 2, '')])], seq, name=f'r{i}'), seq)[0]
                 for i in range(6)]
        args = argparse.Namespace(max_node_reads=3, min_allele_bq=10, min_variants=1, snv_min_af=0.0,
                                  indel_min_af=0.0, min_af=0.0, variant_type='all', rows=4, width=11, debug_rows=False)
        c = Candidate(1, 1, 'C', 'T', 'SNP')
        candidate_work._STATE = ({1: reads}, args, None)
        candidate_work._CACHE = candidate_work.NodeIndexCache(0, 0)
        first = process_node((1, [c]))[0][0][1]
        candidate_work._STATE = ({1: list(reversed(reads))}, args, None)
        second = process_node((1, [c]))[0][0][1]
        self.assertEqual(first['coverage'], 3)
        self.assertEqual(first, second)


class EarlyAfBuildTest(unittest.TestCase):
    def run_build(self, root, name, **extra):
        gam = write_gam(root / 'mixed.gam', mixed_af_rows())
        (root / 'nodes.txt').write_text('10\n20\n')
        (root / 'nodes.json').write_text(json.dumps([dict(node_id=n, sequence=s) for n, s in SEQ.items()]))
        options = dict(command='build', format='candidate-v2', gam=str(gam), index=None,
            nodes=str(root / 'nodes.txt'), node_json=str(root / 'nodes.json'), node_sqlite=None, gfa=None,
            batch_nodes=2, max_node_span=10000, max_batch_segments=100, shard_size=2,
            min_mapq=10, min_af=.3, min_variants=1, min_allele_bq=10., variant_type='all',
            max_indel_len=50, rows=4, width=101, debug_rows=True, workers=1, early_alt_filter=True,
            decode_mode='full', gam_cache_mb=1, output=str(root / name))
        options.update(extra)
        args = argparse.Namespace(**options)
        with redirect_stdout(io.StringIO()):
            build(args)
        return root / name

    def test_accepted_outputs_identical_and_same_candidates_filtered(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = self.run_build(root, 'base', early_af_filter=False)
            fast = self.run_build(root, 'fast', early_af_filter=True)
            manifest = json.loads((fast / 'manifest.json').read_text())
            self.assertGreater(manifest['candidate_optimization']['early_af_rejected'], 0)
            self.assertEqual((base / 'variant_summary.ndjson').read_bytes(), (fast / 'variant_summary.ndjson').read_bytes())
            shards = sorted(p.name for p in base.glob('shard_*_data.npy'))
            self.assertEqual(shards, sorted(p.name for p in fast.glob('shard_*_data.npy')))
            for name in shards:
                self.assertEqual((base / name).read_bytes(), (fast / name).read_bytes(), name)
            ids = lambda p: sorted(json.loads(l)['candidate_id'] for l in (p / 'filtered_candidates.ndjson').open())
            self.assertEqual(ids(base), ids(fast))
            for line in (fast / 'filtered_candidates.ndjson').open():
                meta = json.loads(line)
                if meta.get('support_not_evaluated'):
                    self.assertEqual(meta['reasons'], ['min_af'])
                    self.assertLess(meta['af_upper_bound'], .3)

    def test_af_filter_is_default_and_independent_of_alt_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = self.run_build(root, 'base', early_alt_filter=False, early_af_filter=False)
            fast = self.run_build(root, 'fast', early_alt_filter=False)  # AF filter on by default
            opt = json.loads((fast / 'manifest.json').read_text())['candidate_optimization']
            self.assertTrue(opt['early_af_filter'])
            self.assertEqual(opt['early_rejected'], 0)
            self.assertGreater(opt['early_af_rejected'], 0)
            self.assertEqual((base / 'variant_summary.ndjson').read_bytes(), (fast / 'variant_summary.ndjson').read_bytes())
            for p in base.glob('shard_*_data.npy'):
                self.assertEqual(p.read_bytes(), (fast / p.name).read_bytes(), p.name)

"""Site-level candidates and the default per-node read cap."""
import argparse
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from indexed_gam_pipeline.full_run import validate_shards
from indexed_gam_pipeline.run import build
from test_candidates import alignment
from test_early_af_filter import write_gam
from test_graph_index import fixture_index

SEQ = {10: 'ACGTACGTAC', 20: 'TTGCAAGGCT'}


def site_rows():
    """Node 10: position 4 A>C x4, A>G x3, A>T x1 (fails min_variants); position 7
    DEL x3 and INS GG x2. Node 20: position 3 C>A x3."""
    snv = ['C'] * 4 + ['G'] * 3 + ['T'] + [None] * 3
    indel = {0: 'D', 4: 'D', 8: 'D', 1: 'I', 9: 'I'}
    rows = []
    for i, base in enumerate(snv):
        edits = [(4, 4, ''), (1, 1, base)] if base else [(4, 4, ''), (1, 1, '')]
        edits += {'D': [(2, 2, ''), (1, 0, ''), (2, 2, '')],
                  'I': [(2, 2, ''), (0, 2, 'GG'), (3, 3, '')]}.get(indel.get(i), [(5, 5, '')])
        rows.append(alignment([(10, 0, False, edits)], SEQ, name=f'n10_{i}'))
    for i in range(5):
        edits = [(3, 3, ''), (1, 1, 'A'), (6, 6, '')] if i < 3 else [(10, 10, '')]
        rows.append(alignment([(20, 0, False, edits)], SEQ, name=f'n20_{i}'))
    return rows


def run(root, name, nodes=(10, 20), **extra):
    gam = root / 'sites.gam'
    if not gam.exists():
        write_gam(gam, site_rows())
        fixture_index(root / 'graph.sqlite', [(n, s, 7) for n, s in SEQ.items()])
    (root / f'{name}.nodes').write_text(''.join(f'{n}\n' for n in nodes))
    options = dict(command='build', format='candidate-v4', gam=str(gam), index=None,
        nodes=str(root / f'{name}.nodes'), graph_index=str(root / 'graph.sqlite'), node_json=None,
        node_sqlite=None, gfa=None, gbz=None, gbz_query=None, occurrence_cache=None,
        output=str(root / name / 'shared'), snv_output=str(root / name / 'SNV'),
        indel_output=str(root / name / 'INDEL'), snv_min_af=.2, indel_min_af=.1,
        batch_nodes=512, max_node_span=10000, max_batch_segments=100, shard_size=2048,
        min_mapq=10, min_af=.05, min_variants=2, min_allele_bq=10., variant_type='all',
        max_indel_len=50, rows=200, width=101, debug_rows=False, workers=1, early_alt_filter=True)
    options.update(extra)
    with redirect_stdout(io.StringIO()):
        build(argparse.Namespace(**options))
    return root / name


def records(folder):
    metas = [json.loads(line) for line in (folder / 'variant_summary.ndjson').open()]
    shards = sorted(folder.glob('shard_*_data.npy'))
    tensors = np.concatenate([np.load(p) for p in shards]) if shards else np.zeros((0,))
    return metas, tensors


class SiteTest(unittest.TestCase):
    def test_one_tensor_per_site_equal_to_representative_allele_tensor(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            site, allele = run(root, 'site'), run(root, 'allele', candidate_unit='allele')
            for kind, expected in (('SNV', {'10:4:SNV': ['A>C', 'A>G'], '20:3:SNV': ['C>A']}),
                                   ('INDEL', {'10:7:INDEL': ['DEL', 'INS']})):
                smeta, stensor = records(site / kind)
                ameta, atensor = records(allele / kind)
                got = {m['site_id']: [a['ref'] + '>' + a['alt'] if kind == 'SNV' else a['event_type']
                                      for a in m['alleles']] for m in smeta}
                self.assertEqual(got, expected)
                self.assertEqual(len(ameta), sum(map(len, expected.values())))
                by_id = {m['candidate_id']: atensor[i] for i, m in enumerate(ameta)}
                for i, m in enumerate(smeta):
                    self.assertEqual(m['sample_unit'], 'site-v1')
                    self.assertEqual(m['candidate_id'], m['alleles'][0]['candidate_id'])
                    np.testing.assert_array_equal(stensor[i], by_id[m['candidate_id']])
                    self.assertEqual(m['second_allele_af'], m['alleles'][1]['af'] if len(m['alleles']) > 1 else 0)
                # Every candidate is accounted for exactly once in both modes.
                filtered = lambda f: [json.loads(l)['candidate_id'] for l in (f / kind / 'filtered_candidates.ndjson').open()]
                site_ids = [a['candidate_id'] for m in smeta for a in m['alleles']] + filtered(site)
                allele_ids = [m['candidate_id'] for m in ameta] + filtered(allele)
                self.assertEqual(len(site_ids), len(set(site_ids)))
                self.assertEqual(sorted(site_ids), sorted(allele_ids))
                report = validate_shards(site / kind, 2048)
                self.assertEqual(report['sites'], len(expected))
            snv = records(site / 'SNV')[0]
            first = next(m for m in snv if m['site_id'] == '10:4:SNV')
            self.assertEqual([a['alt_count'] for a in first['alleles']], [4, 3])
            self.assertEqual(first['alt'], 'C')

    def test_duplicate_site_is_rejected_by_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = run(Path(tmp), 'site') / 'SNV'
            lines = (folder / 'variant_summary.ndjson').read_text().splitlines()
            second = json.loads(lines[1]); first = json.loads(lines[0])
            second.update(site_id=first['site_id'], node_id=first['node_id'], start=first['start'],
                          alleles=first['alleles'], allele_count=first['allele_count'],
                          candidate_id=first['candidate_id'], event_type=first['event_type'])
            lines[1] = json.dumps(second)
            (folder / 'variant_summary.ndjson').write_text('\n'.join(lines) + '\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate site'):
                validate_shards(folder, 2048)

    def test_read_cap_defaults_to_800(self):
        from indexed_gam_pipeline.candidate_work import node_reads
        reads = [type('R', (), {'digest': f'{i:04d}'})() for i in range(900)]
        kept = node_reads({1: reads}, 1, argparse.Namespace())
        self.assertEqual([r.digest for r in kept], [f'{i:04d}' for i in range(800)])
        self.assertEqual(len(node_reads({1: reads}, 1, argparse.Namespace(max_node_reads=0))), 900)


if __name__ == '__main__':
    unittest.main()

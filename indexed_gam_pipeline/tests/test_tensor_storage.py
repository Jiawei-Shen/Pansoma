"""Storage saturation must not wrap bytes or change raw candidate evidence."""
import unittest
import numpy as np
from indexed_gam_pipeline.candidates import decode_alignment, make_tensor
from indexed_gam_pipeline.tensor_storage import STORAGE_VERSION, manifest_dtype
from test_candidates import alignment, eligible


class TensorStorageTest(unittest.TestCase):
    def test_path_count_boundaries_preserve_evidence_and_padding(self):
        seq = {1: 'ACA'}
        a = alignment([(1, 0, False, [(1, 1, ''), (1, 1, 'T'), (1, 1, '')])], seq)
        read, _ = decode_alignment(a, seq)
        c = read.observations[0].candidate
        baseline = None
        for count in (0, 1, 3, 4, 7, 86, 127, 128, 254, 255, 256, 331, 507, 508, 511, 512, 40000):
            with self.subTest(count=count):
                x, m = make_tensor(c, eligible(c, [read]), rows=2, width=5,
                                   node_walk_counts={1: count})
                self.assertEqual(x.dtype, np.int8)
                self.assertEqual(x.nbytes, 7 * 2 * 5)
                np.testing.assert_array_equal(x[6, 0], [0, min(count // 4,127), min(count // 4,127), min(count // 4,127), 0])
                self.assertFalse(x[:, 1].any())
                evidence = (m['coverage'], m['alt_count'], m['af'])
                baseline = evidence if baseline is None else baseline
                self.assertEqual(evidence, baseline)

    def test_quality_clips_without_modifying_raw_mapq(self):
        seq = {1: 'ACA'}
        for quality, expected in ((b'', -1), (bytes([0]*3), 0), (bytes([127]*3), 127), (bytes([255]*3), 127)):
            for mapq in (0, 127, 255, 256, 40000):
                with self.subTest(quality=quality, mapq=mapq):
                    a = alignment([(1, 0, False, [(1,1,''),(1,1,'T'),(1,1,'')])], seq)
                    a.quality = quality
                    a.mapping_quality = mapq
                    r, _ = decode_alignment(a, seq)
                    c = r.observations[0].candidate
                    x, _ = make_tensor(c, eligible(c, [r]), rows=1, width=3,
                                       node_walk_counts={1: 331})
                    self.assertTrue(np.all(x[1] == expected))
                    self.assertTrue(np.all(x[3] == min(mapq,127)))
                    self.assertEqual(r.mapq, mapq)

    def test_manifest_distinguishes_new_encoding_and_legacy(self):
        self.assertEqual(manifest_dtype({}), np.dtype('int32'))
        self.assertEqual(manifest_dtype(dict(dtype='int8', tensor_storage_version=STORAGE_VERSION)), np.dtype('int8'))
        for m in (dict(dtype='uint8'), dict(dtype='int8'), dict(dtype='uint8',tensor_storage_version='unknown')):
            with self.assertRaises(ValueError):
                manifest_dtype(m)

import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch

from ..chunks import ChunkLoader, build_cache, cache_name
from ..data import EpochSampler, TensorDataset, load_parts
from .fixtures import make_tensor_set


class ChunkLoaderTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        root = Path(cls.tmp.name) / "s" / "tensors"
        make_tensor_set(root, {"SNV": {"chr1": 50, "chr2": 37}}, shard_size=16, seed=2)
        cls.ds = TensorDataset(load_parts([root], ["SNV"], Path(cls.tmp.name) / "c", labelled=False))
        cls.x = np.stack([cls.ds[k][0].numpy() for k in range(len(cls.ds))])

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def sampler(self, **kw):
        keep = np.arange(len(self.ds)) % 5 == 0
        return EpochSampler(len(self.ds), shuffle=True, seed=3, labels=self.ds.label, non_fraction=0.5, keep=keep, **kw)

    def check(self, batches, want):
        """Every batch holds the dataset's own tensor, blocks, scalars and label for its index."""
        got = []
        for x, blocks, scalars, y in batches:
            self.assertEqual((x.dtype, blocks.dtype, scalars.dtype, y.dtype), (torch.int8, torch.int16, torch.float32, torch.int64))
            for k in range(len(y)):  # find the index by the tensor itself (tensors are distinct)
                i = int(np.flatnonzero((self.x == x[k].numpy()).all(axis=(1, 2, 3)))[0])
                self.assertEqual(blocks[k].tolist(), self.ds.blocks[i].tolist())
                self.assertTrue(np.array_equal(scalars[k].numpy(), self.ds.scalars[i]))
                self.assertEqual(int(y[k]), int(self.ds.label[i]))
                got.append(i)
        self.assertEqual(sorted(got), sorted(want))
        return got

    def test_one_rank_yields_the_epoch_exactly_once(self):
        for chunk_rows, window in ((4, 2), (5, 1), (64, 8)):
            s = self.sampler()
            loader = ChunkLoader(self.ds, s, batch_size=7, threads=3, chunk_rows=chunk_rows, window=window, seed=1)
            for epoch in (0, 1):
                loader.set_epoch(epoch)
                batches = list(loader)
                self.assertEqual(len(batches), len(loader))
                self.assertTrue(all(len(b[3]) == 7 for b in batches[:-1]))
                got = self.check(batches, s.epoch_indices().tolist())
                if epoch == 0:
                    first = got
            self.assertNotEqual(first, got)                                   # a new order each epoch

    def test_ranks_get_equal_shares_of_the_epoch(self):
        s = self.sampler()
        loaders = [ChunkLoader(self.ds, s, batch_size=5, threads=2, chunk_rows=4, window=3, rank=r, world=3, seed=1)
                   for r in range(3)]
        per = [list(l) for l in loaders]
        n = [sum(len(b[3]) for b in p) for p in per]
        self.assertEqual(len(set(n)), 1)
        self.assertEqual(len({len(p) for p in per}), 1)                      # the same number of steps
        want = set(s.epoch_indices().tolist())
        got = [i for p in per for i in self.check(p, [int(np.flatnonzero((self.x == b[0][k].numpy()).all(axis=(1, 2, 3)))[0]) for b in p for k in range(len(b[3]))])]
        self.assertTrue(set(got) <= want)
        self.assertGreaterEqual(len(set(got)), len(want) - 4)                # at most one chunk's worth dropped

    def test_the_zstd_cache_gives_the_same_tensors(self):
        cache = Path(self.tmp.name) / "chunks"
        files = sorted(set(self.ds.files))
        stems = build_cache(files, cache, chunk_rows=5, threads=2, log=lambda m: None)
        self.assertEqual(len(list(cache.glob("*.zst"))), len(files))
        again = build_cache(files, cache, chunk_rows=5, threads=2, log=lambda m: self.fail("rebuilt: " + m))
        self.assertEqual(stems, again)                                       # reused
        self.assertNotEqual(cache_name(files[0], 5), cache_name(files[0], 6))
        s = self.sampler()
        loader = ChunkLoader(self.ds, s, batch_size=6, threads=3, chunk_rows=5, window=2, seed=4, cache=stems)
        self.check(list(loader), s.epoch_indices().tolist())
        with self.assertRaises(ValueError):                                  # a cache of another chunk size
            list(ChunkLoader(self.ds, EpochSampler(len(self.ds), shuffle=False), 4, chunk_rows=16, cache=stems))

    def test_a_read_error_reaches_the_consumer(self):
        s = self.sampler()
        loader = ChunkLoader(self.ds, s, batch_size=4, threads=2, chunk_rows=4, window=2)
        loader.headers = [(10 ** 12, 0)] * len(self.ds.files)                # offsets past the end: short reads
        with self.assertRaises(IOError):
            list(loader)


if __name__ == "__main__":
    unittest.main()

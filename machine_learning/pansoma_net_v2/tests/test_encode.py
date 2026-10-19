import unittest

import numpy as np
import torch

from ..encode import (CATEGORICAL, CONTINUOUS, DATASET_NORM, DATASET_NORM_BEFORE_20261019, N_PLANES, PLANES, WIDTH,
                      TensorEncoder, compute_stats)
from .fixtures import make_tensor

STATS = {"base_quality": {"mean": 30.0, "std": 10.0}, "mapping_quality": {"mean": 40.0, "std": 20.0},
         "path_count": {"mean": 80.0, "std": 15.0}}


def batch(n=6, seed=1):
    rng = np.random.default_rng(seed)
    xs, bs = [], []
    for _ in range(n):
        rows = int(rng.integers(3, 30))
        a1 = int(rng.integers(1, rows))
        blocks = [a1, a1 + 1 if a1 + 1 < rows else a1, rows - 1, rows]
        xs.append(make_tensor(rng, blocks))
        bs.append(blocks)
    return torch.from_numpy(np.stack(xs)), torch.tensor(bs, dtype=torch.int16)


def plane(out, name):
    return out[:, PLANES.index(name)]


class EncoderTest(unittest.TestCase):
    def setUp(self):
        self.x, self.blocks = batch()
        self.enc = TensorEncoder(STATS)
        self.out = self.enc(self.x, self.blocks)
        self.covered = self.x[:, 4] > 0

    def test_shape_and_plane_names(self):
        self.assertEqual(N_PLANES, 36)
        self.assertEqual(tuple(self.out.shape), (6, 36, 200, WIDTH))
        self.assertEqual(len(set(PLANES)), N_PLANES)

    def test_dropped_planes_are_left_out(self):
        enc = TensorEncoder(STATS, drop=["path_count"])
        out = enc(self.x, self.blocks)
        self.assertEqual((enc.n_planes, out.shape[1]), (N_PLANES - 1, N_PLANES - 1))
        kept = [p for p in PLANES if p != "path_count"]
        for k, name in enumerate(kept):
            self.assertTrue(torch.equal(out[:, k], plane(self.out, name)))
        self.x[:, 6] = 1                                    # the path count no longer reaches the output
        self.assertTrue(torch.equal(enc(self.x, self.blocks), out))
        with self.assertRaises(ValueError):
            TensorEncoder(STATS, drop=["no_such_plane"])

    def test_padding_is_zero_in_every_plane(self):
        out = self.out[..., :101]
        self.assertTrue(((~self.covered).sum() > 0).item())
        self.assertEqual(float(out.permute(1, 0, 2, 3)[:, ~self.covered].abs().sum()), 0.0)
        self.assertEqual(float(self.out[..., 101:].abs().sum()), 0.0)  # width padding for the 4 x 4 stem

    def test_one_hot_groups(self):
        out = self.out[..., :101]
        for name, c, k in CATEGORICAL:
            group = torch.stack([plane(out, f"{name}={code}") for code in range(1, k + 1)], 1)
            for code in range(1, k + 1):
                self.assertTrue(torch.equal(group[:, code - 1].bool(), self.x[:, c] == code), (name, code))
            valid = (self.x[:, c] >= 1) & (self.x[:, c] <= k)
            self.assertTrue(torch.equal(group.sum(1), valid.float()), name)  # one plane per cell, none for 0
        self.assertTrue(torch.all(plane(out, "strand=1") + plane(out, "strand=2") == self.covered.float()))

    def test_masked_zscore_keeps_zero_meaning(self):
        out = self.out[..., :101]
        for name, c, low in CONTINUOUS:
            z = plane(out, name)
            valid = self.covered & (self.x[:, c] >= low)
            expected = (self.x[:, c].float() - STATS[name]["mean"]) / STATS[name]["std"]
            self.assertTrue(torch.allclose(z[valid], expected[valid]), name)
            self.assertEqual(float(z[~valid].abs().sum()), 0.0, name)
        bq = self.x[:, 1]
        self.assertTrue(torch.equal(plane(out, "bq_missing").bool(), self.covered & (bq < 0)))
        # MAPQ 0 on a read is a value (z = -2 here), distinct from padding; the coverage plane marks it
        mapq0 = self.covered & (self.x[:, 3] == 0)
        if mapq0.any():
            self.assertTrue(torch.all(plane(out, "mapping_quality")[mapq0] == -2.0))
        self.assertTrue(torch.equal(plane(out, "covered").bool(), self.covered))

    def test_differs(self):
        out = self.out[..., :101]
        self.assertTrue(torch.equal(plane(out, "differs").bool(), self.covered & (self.x[:, 0] != self.x[:, 5])))

    def test_row_blocks_fill_the_covered_cells_of_each_read(self):
        out = self.out[..., :101]
        rows = torch.arange(200)
        for i, (a1, alt, ref, other) in enumerate(self.blocks.tolist()):
            spans = {"row_A1": (0, a1), "row_ALT": (a1, alt), "row_REF": (alt, ref), "row_OTHER": (ref, other)}
            for name, (lo, hi) in spans.items():
                in_rows = ((rows >= lo) & (rows < hi))[:, None]
                self.assertTrue(torch.equal(plane(out, name)[i].bool(), self.covered[i] & in_rows), (i, name))
        blocks_sum = sum(plane(out, f"row_{b}") for b in ("A1", "ALT", "REF", "OTHER"))
        self.assertTrue(torch.equal(blocks_sum, self.covered.float()))  # every read in exactly one block

    def test_statistics_come_from_valid_covered_cells_only(self):
        loader = [(self.x[:3], self.blocks[:3], None), (self.x[3:], self.blocks[3:], None)]
        stats = compute_stats(loader, 100)
        self.assertEqual(stats["tensors"], 6)
        for name, c, low in CONTINUOUS:
            v = self.x[:, c][self.covered & (self.x[:, c] >= low)].double()
            self.assertAlmostEqual(stats[name]["mean"], float(v.mean()), places=6)
            self.assertAlmostEqual(stats[name]["std"], float(v.std(unbiased=False)), places=5)
        padded = torch.cat([self.x, torch.zeros_like(self.x)])  # empty tensors change nothing
        again = compute_stats([(padded, None, None)], 100)
        self.assertEqual({k: v for k, v in again.items() if k != "tensors"},
                         {k: v for k, v in stats.items() if k != "tensors"})
        self.assertEqual(compute_stats(loader, 4)["tensors"], 4)

    def test_statistics_round_trip_and_unfitted_encoder(self):
        with self.assertRaises(RuntimeError):
            TensorEncoder()(self.x, self.blocks)
        enc = TensorEncoder()
        enc.load_state_dict(self.enc.state_dict())
        self.assertEqual(enc.stats(), self.enc.stats())
        self.assertTrue(torch.equal(enc(self.x, self.blocks), self.out))

    def test_dataset_statistics_remove_a_dataset_offset(self):
        """norm dataset: a dataset whose DATASET_NORM channels (BQ) are shifted and scaled encodes like the original when
        each is z-scored with its own statistics; path count keeps the encoder's statistics; without a table part is
        ignored."""
        x2 = self.x.clone().to(torch.int16)
        for name, c, low in CONTINUOUS:
            if name in DATASET_NORM:  # a "run" with every valid value doubled and shifted by 3
                valid = (x2[:, 4] > 0) & (x2[:, c] >= low)
                x2[:, c] = torch.where(valid, x2[:, c] * 2 + 3, x2[:, c])
        own = [compute_stats([(self.x, self.blocks, None)], 100), compute_stats([(x2, self.blocks, None)], 100)]
        enc = TensorEncoder(STATS)
        enc.set_part_stats(own)
        zeros, ones = torch.zeros(len(self.x), dtype=torch.long), torch.ones(len(self.x), dtype=torch.long)
        a, b = enc(self.x, self.blocks, part=zeros), enc(x2, self.blocks, part=ones)
        for name in DATASET_NORM:
            self.assertTrue(torch.allclose(plane(a, name), plane(b, name), atol=1e-4), name)
        self.assertTrue(torch.equal(plane(a, "path_count"), plane(self.out, "path_count")))  # the encoder's statistics
        self.assertFalse(torch.allclose(plane(self.enc(x2, self.blocks), "base_quality"),
                                        plane(self.out, "base_quality"), atol=1e-2))  # one global z-score keeps the shift
        mixed = enc(torch.cat([self.x, x2]), torch.cat([self.blocks, self.blocks]), part=torch.cat([zeros, ones]))
        self.assertTrue(torch.equal(mixed[:len(self.x)], a) and torch.equal(mixed[len(self.x):], b))  # per sample
        self.assertEqual([p["base_quality"]["mean"] for p in enc.part_stats()],
                         [float(torch.tensor(o["base_quality"]["mean"], dtype=torch.float32)) for o in own])
        self.assertNotIn("part_mean", enc.state_dict())  # the table belongs to the data, not to the checkpoint
        enc.set_part_stats([])
        self.assertTrue(torch.equal(enc(self.x, self.blocks, part=ones), self.out))  # no table: the buffers
        with self.assertRaises(RuntimeError):
            TensorEncoder().set_part_stats(own)  # the encoder's own statistics come first

    def test_dataset_norm_channels(self):
        """dataset_norm picks the channels z-scored per dataset: BQ for new models (MAPQ keeps the encoder's statistics,
        like path count), BQ and MAPQ for the norm dataset checkpoints of 2026-10-13 to 10-19; unknown names fail."""
        self.assertEqual(DATASET_NORM, ("base_quality",))
        x2 = self.x.clone().to(torch.int16)
        for name, c, low in CONTINUOUS:
            if name in DATASET_NORM_BEFORE_20261019:  # BQ and MAPQ of a "run" doubled and shifted by 3
                valid = (x2[:, 4] > 0) & (x2[:, c] >= low)
                x2[:, c] = torch.where(valid, x2[:, c] * 2 + 3, x2[:, c])
        zeros = torch.zeros(len(self.x), dtype=torch.long)
        own, own2 = [compute_stats([(self.x, self.blocks, None)], 100)], [compute_stats([(x2, self.blocks, None)], 100)]
        ref = TensorEncoder(STATS, dataset_norm=DATASET_NORM_BEFORE_20261019)
        ref.set_part_stats(own)
        b = ref(self.x, self.blocks, part=zeros)
        old = TensorEncoder(STATS, dataset_norm=DATASET_NORM_BEFORE_20261019)
        old.set_part_stats(own2)
        new = TensorEncoder(STATS)
        new.set_part_stats(own2)
        a_old, a_new = old(x2, self.blocks, part=zeros), new(x2, self.blocks, part=zeros)
        for name in DATASET_NORM_BEFORE_20261019:
            self.assertTrue(torch.allclose(plane(a_old, name), plane(b, name), atol=1e-4), name)
        self.assertTrue(torch.allclose(plane(a_new, "base_quality"), plane(b, "base_quality"), atol=1e-4))
        self.assertTrue(torch.equal(plane(a_new, "mapping_quality"), plane(self.enc(x2, self.blocks), "mapping_quality")))
        self.assertEqual(new.dataset_norm, ("base_quality",))
        with self.assertRaises(ValueError):
            TensorEncoder(STATS, dataset_norm=("depth",))

    def test_bfloat16_output(self):
        out = self.enc(self.x, self.blocks, dtype=torch.bfloat16)
        self.assertEqual(out.dtype, torch.bfloat16)
        self.assertTrue(torch.equal(out.float() != 0, self.out != 0))


if __name__ == "__main__":
    unittest.main()

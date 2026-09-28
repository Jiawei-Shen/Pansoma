import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from .. import data
from ..data import (SCALARS, EpochSampler, KindIndex, TensorDataset, block_ends, block_split, load_parts, site_scalars,
                    somatic_truth, validation_truth)
from .fixtures import make_tensor_set

SPEC = {"SNV": {"chr1": 9, "chr2": 6}, "INDEL": {"chr1": 5, "chr2": 3}}


class IndexTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name) / "sample" / "v3_tensors"
        self.truth = make_tensor_set(self.root, SPEC)
        self.cache = Path(self.tmp.name) / "cache"

    def tearDown(self):
        self.tmp.cleanup()

    def test_index_matches_the_directory(self):
        for kind, rows in self.truth.items():
            index = KindIndex(self.root / kind, self.cache)
            self.assertEqual(len(index), len(rows))
            a = index.arrays
            for k, t in enumerate(rows):
                self.assertEqual(index.meta["chroms"][a["chrom"][k]], t["chrom"])
                self.assertEqual(index.meta["shards"][a["shard"][k]], t["file"])
                self.assertEqual(int(a["row"][k]), t["row"])
                self.assertEqual(int(a["label"][k]), t["label"])
                self.assertEqual(a["blocks"][k].tolist(), t["blocks"])
                self.assertEqual(int(a["node"][k]), t["node"])
                self.assertEqual(index.reason(k), t["reason"])
                self.assertEqual(bool(a["off_reference"][k]), t["off_reference"])
                self.assertEqual(int(a["eval_label"][k]), t["eval_label"])
                self.assertEqual(bool(a["in_region"][k]), t["in_region"])
                self.assertEqual(int(a["pos0"][k]), t["node"])  # GRCh38 position (anchor middle off the reference)
                self.assertEqual(index.truth_of(k).tolist(), [t["truth_id"]] if t["label"] == 1 else [])
                # the record's own numbers, not those of an entry in alleles[]
                self.assertTrue(np.allclose(a["scalars"][k], site_scalars(t["numbers"]), atol=1e-5))
                self.assertAlmostEqual(float(a["scalars"][k][SCALARS.index("af")]), t["af"], places=4)
            self.assertEqual(index.candidates(), [t["candidate_id"] for t in rows])

    def test_cache_is_reused_and_rebuilt_when_labels_change(self):
        KindIndex(self.root / "SNV", self.cache)
        with mock.patch.object(data, "build_index", side_effect=AssertionError("rebuilt")):
            cached = KindIndex(self.root / "SNV", self.cache)
        self.assertEqual(cached.candidates(), [t["candidate_id"] for t in self.truth["SNV"]])
        m = self.root / "SNV" / "labels.manifest.json"
        meta = json.loads(m.read_text())
        meta["created"] = "2026-09-29T00:00:00+00:00"
        m.write_text(json.dumps(meta))
        with mock.patch.object(data, "build_index", wraps=data.build_index) as build:
            KindIndex(self.root / "SNV", self.cache)
        self.assertEqual(build.call_count, 1)

    def test_selection_by_chromosome_leaves_out_ignored_tensors(self):
        index = KindIndex(self.root / "SNV", self.cache)
        labels = [t["label"] for t in self.truth["SNV"]]
        chroms = [t["chrom"] for t in self.truth["SNV"]]
        chr2 = index.select(["chr2"])
        self.assertEqual(chr2.tolist(), [k for k in range(len(labels)) if chroms[k] == "chr2" and labels[k] >= 0])
        self.assertEqual(len(index.select(["chr2"], labelled=False)), SPEC["SNV"]["chr2"])
        self.assertIn(-1, labels)  # the fixture has ignored tensors, and they are never selected for training
        self.assertTrue(all(index.arrays["label"][index.select()] >= 0))
        # evaluation adds exactly the off-reference tensors without a truth match, as non
        rows = self.truth["SNV"]
        evaluation = index.select(evaluation=True)
        self.assertEqual(evaluation.tolist(), [k for k, t in enumerate(rows) if t["eval_label"] >= 0])
        added = [k for k in evaluation if rows[k]["label"] == -1]
        self.assertTrue(added and all(rows[k]["reason"] == "off_reference_no_truth_match" for k in added))
        self.assertTrue(all(index.arrays["eval_label"][added] == 0))
        others = [k for k, t in enumerate(rows) if t["reason"] in ("outside_confident_region", "below_snv_min_af")]
        self.assertFalse(set(others) & set(evaluation.tolist()))
        # labels 1 and 2 outside the confident region are trained on but not evaluated (a BED filter drops them)
        outside = [k for k, t in enumerate(rows) if t["label"] in (1, 2) and not t["in_region"]]
        self.assertTrue(outside)
        self.assertFalse(set(outside) & set(evaluation.tolist()))
        self.assertTrue(set(outside) <= set(index.select().tolist()))

    def test_dataset_returns_the_stored_tensor_blocks_and_label(self):
        parts = load_parts([self.root], ["SNV", "INDEL"], self.cache, labelled=False)
        ds = TensorDataset(parts)
        self.assertEqual(len(ds), sum(len(v) for v in self.truth.values()))
        flat = self.truth["SNV"] + self.truth["INDEL"]
        for k in (0, 5, len(ds) - 1):
            x, blocks, scalars, label = ds[k]
            self.assertTrue(np.array_equal(x.numpy(), flat[k]["x"]))
            self.assertEqual(blocks.tolist(), flat[k]["blocks"])
            self.assertTrue(np.allclose(scalars.numpy(), site_scalars(flat[k]["numbers"]), atol=1e-5))
            self.assertEqual(label, flat[k]["label"])
        counts = ds.class_counts()
        self.assertEqual(counts.tolist(), [sum(t["label"] == c for t in flat) for c in range(3)])

    def test_two_samples_with_the_same_set_name_share_a_cache(self):
        other = Path(self.tmp.name) / "sample2" / "v3_tensors"
        truth2 = make_tensor_set(other, {"SNV": {"chr1": 3}}, seed=5)
        a, b = KindIndex(self.root / "SNV", self.cache), KindIndex(other / "SNV", self.cache)
        with mock.patch.object(data, "build_index", side_effect=AssertionError("rebuilt")):
            a2, b2 = KindIndex(self.root / "SNV", self.cache), KindIndex(other / "SNV", self.cache)
        self.assertEqual(len(a2), len(self.truth["SNV"]))
        self.assertEqual(b2.candidates(), [t["candidate_id"] for t in truth2["SNV"]])
        self.assertEqual((len(a), len(b)), (len(a2), len(b2)))

    def test_block_split_keeps_whole_blocks(self):
        index = KindIndex(self.root / "SNV", self.cache)
        pos = index.select(labelled=False)
        train, val = block_split(index, pos, 0.5, block_nodes=14, seed=1)
        self.assertEqual(sorted(np.r_[train, val].tolist()), pos.tolist())
        key = lambda p: (int(index.arrays["chrom"][p]), int(index.arrays["node"][p]) // 14)  # noqa: E731
        self.assertFalse({key(p) for p in train} & {key(p) for p in val})
        again = block_split(index, pos, 0.5, block_nodes=14, seed=1)
        self.assertEqual(val.tolist(), again[1].tolist())                   # deterministic
        self.assertEqual(len(block_split(index, pos, 0.0, 14)[1]), 0)
        self.assertEqual(len(block_split(index, pos, 1.0, 14)[0]), 0)

    def test_truth_of_the_validation_blocks(self):
        root = Path(self.tmp.name) / "big" / "v3_tensors"
        rows = make_tensor_set(root, {"SNV": {"chr1": 60, "chr2": 60}}, shard_size=16, seed=3)["SNV"]
        index = KindIndex(root / "SNV", self.cache)
        truth = somatic_truth(index)                                          # PASS, in the BED, SNP rows
        other = somatic_truth(index, other_kind=True)                         # DEL and INS rows
        self.assertFalse(set(truth) & set(other))
        in_region = [t["truth_id"] for t in rows if t["label"] == 1 and t["in_region"]]
        self.assertTrue(all(tid in truth or tid in other for tid in in_region))
        self.assertTrue(any(tid in other for tid in in_region))              # SNV tensors of an INS truth
        self.assertTrue(any(tid not in {t["truth_id"] for t in rows} for tid in truth))  # truth without a tensor
        pos = index.select(labelled=False)
        train, val = block_split(index, pos, 0.5, block_nodes=21, seed=2)
        vt = validation_truth(index, pos, truth, 0.5, block_nodes=21, seed=2)
        val_ids = {rows[k]["truth_id"] for k in val.tolist() if rows[k]["label"] == 1} & set(truth)
        train_ids = {rows[k]["truth_id"] for k in train.tolist() if rows[k]["label"] == 1} & set(truth)
        self.assertTrue(val_ids and val_ids <= vt)                            # a validation tensor's truth is validation
        self.assertFalse((train_ids - val_ids) & vt)                          # a training tensor's truth is not
        self.assertEqual(vt, validation_truth(index, pos, truth, 0.5, block_nodes=21, seed=2))  # deterministic

    def test_block_ends(self):
        groups = [dict(start_row=0, end_row=4, allele="A1"), dict(start_row=4, end_row=6, allele="A2"),
                  dict(start_row=6, end_row=7, allele="A3"), dict(start_row=7, end_row=9, allele="OTHER")]
        self.assertEqual(block_ends(groups), [4, 7, 7, 9])  # no REF rows: the REF block is empty
        self.assertEqual(block_ends([dict(start_row=0, end_row=3, allele="A1")]), [3, 3, 3, 3])
        with self.assertRaises(ValueError):
            block_ends([dict(start_row=1, end_row=3, allele="A1")])
        with self.assertRaises(ValueError):
            block_ends([dict(start_row=0, end_row=2, allele="REF"), dict(start_row=2, end_row=3, allele="A1")])


class SamplerTest(unittest.TestCase):
    def test_kept_non_every_epoch_and_the_rest_sampled(self):
        labels = np.array([0] * 100 + [0] * 20 + [1] * 5 + [2] * 5)
        keep = np.array([False] * 100 + [True] * 20 + [False] * 10)      # e.g. the high-AF non
        s = EpochSampler(130, seed=1, labels=labels, non_fraction=0.25, keep=keep)
        for epoch in range(3):
            s.set_epoch(epoch)
            got = set(s)
            self.assertTrue(set(range(100, 130)) <= got)                  # kept non, somatic, germline: always
            low = [i for i in got if i < 100]
            self.assertEqual(len(low), 25)                                # 25 % of the other non
            self.assertAlmostEqual(len(low) / 0.25, 100)                  # with weight 1/0.25 they stand for all 100
        self.assertEqual(s.epoch_size(), 55)


    def test_evaluation_slices_are_disjoint_and_complete(self):
        parts = [list(EpochSampler(10, shuffle=False, rank=r, world=3)) for r in range(3)]
        self.assertEqual(sorted(sum(parts, [])), list(range(10)))
        self.assertEqual([len(EpochSampler(10, shuffle=False, rank=r, world=3)) for r in range(3)], [4, 3, 3])

    def test_training_order_is_padded_per_rank_and_changes_by_epoch(self):
        samplers = [EpochSampler(10, num_samples=7, rank=r, world=2, seed=3) for r in range(2)]
        orders = [list(s) for s in samplers]
        self.assertEqual([len(o) for o in orders], [4, 4])
        self.assertEqual([len(s) for s in samplers], [4, 4])
        self.assertEqual(len(set(orders[0] + orders[1])), 7)
        first = list(samplers[0])
        samplers[0].set_epoch(1)
        self.assertNotEqual(first, list(samplers[0]))


if __name__ == "__main__":
    unittest.main()

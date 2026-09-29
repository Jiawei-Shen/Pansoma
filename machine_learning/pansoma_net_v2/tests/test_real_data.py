"""Checks on a real merged tensor set (skipped when it is absent): the index against the manifests and
summaries, the row blocks against the tensors' own channels, and the encoding of real tensors."""
import json
import os
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import numpy as np
import torch

from ..data import CLASSES, KindIndex, TensorDataset, site_scalars
from ..encode import PLANES, TensorEncoder, compute_stats

REAL = Path(os.environ.get("PANSOMA_TEST_TENSORS", "/scratch/jshen/data/pansoma_v2_tensors/COLO829T_Illumina/tensors"))
CODES = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5}


@unittest.skipUnless((REAL / "SNV" / "manifest.json").exists(), f"no merged tensor set at {REAL}")
class RealDataTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.index = {k: KindIndex(REAL / k, cls.tmp.name) for k in ("SNV", "INDEL")}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_index_counts_match_the_manifests(self):
        for kind, index in self.index.items():
            manifest = json.loads((REAL / kind / "manifest.json").read_text())
            labels = json.loads((REAL / kind / "labels.manifest.json").read_text())
            self.assertEqual(len(index), manifest["tensors"])
            counts = Counter(index.arrays["label"].tolist())
            names = {v: k for k, v in labels["labels"].items()}
            self.assertEqual({names[v]: n for v, n in counts.items()}, labels["totals"])
            reasons = Counter(index.meta["reasons"][r] for r in index.arrays["reason"].tolist())
            self.assertEqual(dict(reasons), labels["reasons"])  # every labels.ndjson line read in order
            a = index.arrays
            to_non = (a["eval_label"] == 0) & (a["label"] == -1)
            self.assertEqual(int(to_non.sum()), labels["reasons"].get("off_reference_no_truth_match", 0))
            self.assertTrue(a["off_reference"][to_non].all())
            # the region test agrees with the labeller: its outside_confident_region tensors are outside, its 0s and
            # its off-reference negatives inside; evaluation drops only labelled tensors outside the region
            reason = np.array(index.meta["reasons"])[a["reason"]]
            self.assertFalse(a["in_region"][reason == "outside_confident_region"].any())
            # (except truth_matches_non_representative_allele: the labeller decides it before the BED test, so it
            # can lie outside; a BED filter drops it at test time, and evaluation leaves it out)
            self.assertTrue(a["in_region"][(a["label"] == 0) & (reason != "truth_matches_non_representative_allele")].all())
            self.assertTrue(a["in_region"][reason == "off_reference_no_truth_match"].all())
            dropped = (a["eval_label"] == -1) & (a["label"] >= 0)
            self.assertFalse(a["in_region"][dropped].any())
            # truth ids: every labelled-1 tensor stands for at least one truth allele of the recall table
            with open(REAL / "somatic.recall.tsv") as f:
                known = {int(line.split("\t")[0]) for line in list(f)[1:]}
            n = np.diff(a["truth_ptr"])
            self.assertTrue((n[a["label"] == 1] >= 1).all())
            self.assertFalse(n[a["label"] != 1].any())
            self.assertTrue(set(a["truth_ids"].tolist()) <= known)

    def test_candidates_follow_the_summary(self):
        index = self.index["SNV"]
        chrom = index.meta["chroms"].index("chr21")
        pos = np.flatnonzero(index.arrays["chrom"] == chrom)
        with open(REAL / "SNV" / "chr21_variant_summary.ndjson") as f:
            records = [json.loads(line) for line in f]
        self.assertEqual([index.candidates()[p] for p in pos], [r["candidate_id"] for r in records])
        for p, r in list(zip(pos, records))[:500]:  # the scalars are the record's own numbers
            self.assertTrue(np.allclose(index.arrays["scalars"][p], site_scalars(r), atol=1e-4), r["candidate_id"])
            self.assertEqual(int(index.arrays["node"][p]), r["node_id"])

    def test_row_blocks_agree_with_the_tensors(self):
        index = self.index["SNV"]
        rng = np.random.default_rng(0)
        picks = rng.choice(len(index), 300, replace=False)
        ds = TensorDataset([(index, picks)])
        checked = Counter()
        for k, p in enumerate(picks):
            x, blocks, _, _ = ds[k]
            x, (a1, alt, ref, other) = x.numpy(), blocks.tolist()
            has_read = x[4].max(1) > 0
            self.assertTrue(has_read[:other].all() and not has_read[other:].any(), index.candidates()[p])
            a1_base = CODES[index.candidates()[p].split(">")[1]]
            site = x[2, :, 50]
            self.assertTrue((site[:a1] == a1_base).all())                        # A1 rows carry the A1 base
            self.assertTrue(((site[a1:alt] != 0) & (site[a1:alt] != a1_base)).all())  # other ALT alleles
            self.assertTrue((site[alt:ref] == x[5, alt:ref, 50]).all())          # REF rows: the graph base
            self.assertFalse(x[2, ref:other].any())                              # OTHER rows: no site allele
            checked.update(dict(A1=a1, ALT=alt - a1, REF=ref - alt, OTHER=other - ref))
        self.assertTrue(all(checked[b] > 0 for b in ("A1", "ALT", "REF", "OTHER")), checked)

    def test_encoding_of_real_tensors(self):
        index = self.index["INDEL"]
        picks = np.random.default_rng(1).choice(len(index), 64, replace=False)
        ds = TensorDataset([(index, picks)])
        x = torch.stack([ds[k][0] for k in range(len(ds))])
        blocks = torch.stack([ds[k][1] for k in range(len(ds))])
        stats = compute_stats([(x, blocks, None)], 64)
        out = TensorEncoder(stats)(x, blocks)
        covered = x[:, 4] > 0
        self.assertTrue(torch.isfinite(out).all())
        self.assertEqual(float(out[..., :101].permute(1, 0, 2, 3)[:, ~covered].abs().sum()), 0.0)
        for name in ("base_quality", "mapping_quality", "path_count"):
            z = out[:, PLANES.index(name), :, :101][covered]
            self.assertLess(abs(float(z.mean())), 0.2, name)  # centred on its own statistics
        self.assertEqual(len(CLASSES), 3)


if __name__ == "__main__":
    unittest.main()

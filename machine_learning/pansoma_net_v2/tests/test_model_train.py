import gzip
import json
import tempfile
import unittest
from pathlib import Path

import torch
import torch.nn as nn

from .. import predict, train
from ..data import KindIndex
from ..model import PansomaNetV2
from .fixtures import make_tensor_set
from .test_encode import STATS, batch

SMALL = dict(depths=(1, 1, 1, 1), dims=(16, 32, 64, 128), front=(8,))  # CBAM needs >= 16 channels
SMALL_ARGS = ["--depths", "1", "1", "1", "1", "--dims", "16", "32", "64", "128", "--front", "8",
              "--batch-size", "4", "--num-workers", "0", "--amp", "off", "--val-chroms", "chr2", "--stats-samples", "10"]


class ModelTest(unittest.TestCase):
    def test_forward_and_checkpoint_round_trip(self):
        x, blocks = batch(4)
        model = PansomaNetV2(3, stats=STATS, **SMALL).eval()
        with torch.no_grad():
            logits = model(x, blocks)
        self.assertEqual(tuple(logits.shape), (4, 3))
        payload = dict(format="pansoma_net_v2", config=model.config, model_state_dict=model.state_dict())
        again = PansomaNetV2.from_checkpoint(payload).eval()
        self.assertEqual(again.encoder.stats(), model.encoder.stats())
        with torch.no_grad():
            self.assertTrue(torch.equal(again(x, blocks), logits))

    def test_ignored_labels_do_not_enter_the_loss(self):
        logits = torch.randn(5, 3)
        y = torch.tensor([0, 1, -1, 2, -1])
        criterion = nn.CrossEntropyLoss(weight=torch.tensor([1.0, 5.0, 2.0]), ignore_index=-1)
        keep = y >= 0
        self.assertTrue(torch.allclose(criterion(logits, y), criterion(logits[keep], y[keep])))


class LoaderTest(unittest.TestCase):
    def test_forkserver_workers_read_the_same_tensors(self):
        from ..data import EpochSampler, TensorDataset, load_parts
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "s" / "t"
            make_tensor_set(root, {"SNV": {"chr1": 10}}, shard_size=3)
            ds = TensorDataset(load_parts([root], ["SNV"], Path(tmp) / "c", labelled=False))
            loader = train.make_loader(ds, EpochSampler(len(ds), shuffle=False), 2, 4, persistent=False)
            got = torch.cat([x for x, _, _ in loader])
            self.assertTrue(torch.equal(got, torch.stack([ds[k][0] for k in range(len(ds))])))


class TrainPredictTest(unittest.TestCase):
    """One small end-to-end run on CPU: train 2 epochs, resume to 3, predict."""

    def test_train_resume_predict(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "v3_tensors"
            truth = make_tensor_set(root, {"SNV": {"chr1": 24, "chr2": 10}, "INDEL": {"chr1": 12, "chr2": 6}},
                                    shard_size=8)
            out = tmp / "run"
            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "2"] + SMALL_ARGS)
            for f in ("best.pth", "last.pth", "stats.json", "metrics.jsonl", "train.log"):
                self.assertTrue((out / f).exists(), f)
            rows = [json.loads(line) for line in (out / "metrics.jsonl").read_text().splitlines()]
            self.assertEqual([r["epoch"] for r in rows], [1, 2])
            ckpt = torch.load(out / "last.pth", weights_only=False)
            stats = json.loads((out / "stats.json").read_text())
            for name, s in ckpt["stats"].items():
                self.assertAlmostEqual(s["mean"], stats[name]["mean"], places=4)
            labelled_chr1 = sum(1 for kind in truth.values() for t in kind if t["chrom"] == "chr1" and t["label"] >= 0)
            self.assertEqual(sum(d["train"] for d in ckpt["data"]), labelled_chr1)  # -1 never trained on
            self.assertEqual(ckpt["chroms"]["val"], ["chr2"])
            self.assertEqual(sum(sum(r) for r in rows[-1]["val"]["confusion"]),
                             sum(1 for kind in truth.values() for t in kind if t["chrom"] == "chr2" and t["label"] >= 0))

            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "3", "--resume", str(out / "last.pth")]
                       + SMALL_ARGS)
            resumed = torch.load(out / "last.pth", weights_only=False)
            self.assertEqual(resumed["epoch"], 3)
            self.assertEqual(resumed["stats"], ckpt["stats"])  # statistics are not refitted

            pred = tmp / "pred"
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(root), "--output", str(pred),
                          "--num-workers", "0", "--amp", "off", "--batch-size", "5"])
            with gzip.open(pred / "sample.v3_tensors.SNV.predictions.tsv.gz", "rt") as f:
                lines = f.read().splitlines()
            self.assertEqual(len(lines) - 1, len(truth["SNV"]))  # every tensor, -1 included
            self.assertEqual([line.split("\t")[1] for line in lines[1:]], [t["candidate_id"] for t in truth["SNV"]])
            probs = [sum(map(float, line.split("\t")[3:6])) for line in lines[1:]]
            self.assertTrue(all(abs(p - 1) < 1e-3 for p in probs))
            metrics = json.loads((pred / "sample.v3_tensors.SNV.metrics.json").read_text())
            self.assertEqual(metrics["labelled"], sum(t["label"] >= 0 for t in truth["SNV"]))
            self.assertEqual(metrics["ignored"], sum(t["label"] < 0 for t in truth["SNV"]))
            self.assertEqual(len(KindIndex(root / "SNV", pred / "index_cache")), len(truth["SNV"]))


if __name__ == "__main__":
    unittest.main()

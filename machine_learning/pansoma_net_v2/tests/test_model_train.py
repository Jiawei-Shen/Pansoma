import gzip
import json
import tempfile
import unittest
from pathlib import Path

import torch
import torch.nn as nn

from .. import predict, train
from ..data import SCALARS, EpochSampler, KindIndex, TensorDataset, load_parts
from ..model import PansomaNetV2
from .fixtures import make_tensor_set
from .test_encode import STATS, batch

SMALL = dict(depths=(1, 1, 1, 1), dims=(16, 32, 64, 128), front=(8,))  # CBAM needs >= 16 channels
SMALL_ARGS = ["--depths", "1", "1", "1", "1", "--dims", "16", "32", "64", "128", "--front", "8",
              "--batch-size", "4", "--num-workers", "0", "--amp", "off", "--stats-samples", "10",
              "--val-fraction", "0.5", "--val-block-nodes", "21"]
SPEC = {"SNV": {"chr1": 10, "chr2": 24, "chr3": 16}, "INDEL": {"chr1": 6, "chr2": 12, "chr3": 8}}


class ModelTest(unittest.TestCase):
    def test_forward_and_checkpoint_round_trip(self):
        x, blocks = batch(4)
        scalars = torch.randn(4, len(SCALARS))
        for n in (0, len(SCALARS)):
            model = PansomaNetV2(3, stats=STATS, scalars=n, **SMALL).eval()
            if n:
                model.set_scalar_stats(torch.zeros(n), torch.full((n,), 2.0))
            with torch.no_grad():
                logits = model(x, blocks, scalars)
            self.assertEqual(tuple(logits.shape), (4, 3))
            payload = dict(format="pansoma_net_v2", config=model.config, model_state_dict=model.state_dict())
            again = PansomaNetV2.from_checkpoint(payload).eval()
            self.assertEqual(again.encoder.stats(), model.encoder.stats())
            with torch.no_grad():
                self.assertTrue(torch.equal(again(x, blocks, scalars), logits))
        with torch.no_grad():  # the scalars reach the logits
            self.assertFalse(torch.equal(model(x, blocks, scalars), model(x, blocks, scalars + 1)))

    def test_ignored_labels_do_not_enter_the_loss(self):
        logits = torch.randn(5, 3)
        y = torch.tensor([0, 1, -1, 2, -1])
        criterion = nn.CrossEntropyLoss(weight=torch.tensor([1.0, 5.0, 2.0]), ignore_index=-1)
        keep = y >= 0
        self.assertTrue(torch.allclose(criterion(logits, y), criterion(logits[keep], y[keep])))


class LoaderTest(unittest.TestCase):
    def test_forkserver_workers_read_the_same_tensors(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "s" / "t"
            make_tensor_set(root, {"SNV": {"chr1": 10}}, shard_size=3)
            ds = TensorDataset(load_parts([root], ["SNV"], Path(tmp) / "c", labelled=False))
            loader = train.make_loader(ds, EpochSampler(len(ds), shuffle=False), 2, 4, persistent=False)
            got = torch.cat([batch[0] for batch in loader])
            self.assertTrue(torch.equal(got, torch.stack([ds[k][0] for k in range(len(ds))])))


class TrainPredictTest(unittest.TestCase):
    """Small end-to-end runs on CPU: train 2 epochs (chr1 left out, validation from node blocks of chr2/chr3),
    resume to 3, predict chr1; and one run with the scalars."""

    def test_train_resume_predict(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "v3_tensors"
            truth = make_tensor_set(root, SPEC, shard_size=8)
            out = tmp / "run"
            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "2"] + SMALL_ARGS)
            for f in ("best.pth", "last.pth", "stats.json", "metrics.jsonl", "train.log"):
                self.assertTrue((out / f).exists(), f)
            rows = [json.loads(line) for line in (out / "metrics.jsonl").read_text().splitlines()]
            self.assertEqual([r["epoch"] for r in rows], [1, 2])
            ckpt = torch.load(out / "last.pth", weights_only=False)
            self.assertEqual(ckpt["chroms"]["test"], ["chr1"])
            self.assertIsNotNone(ckpt["somatic_threshold"])
            labelled_23 = sum(1 for kind in truth.values() for t in kind if t["chrom"] != "chr1" and t["label"] >= 0)
            val_n = rows[-1]["val"]["tensors"]
            self.assertGreater(val_n, 0)
            self.assertEqual(sum(d["train"] for d in ckpt["data"]) + val_n, labelled_23)  # chr1 and -1 never used
            self.assertIn("somatic_ap", rows[-1]["val"])

            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "3", "--resume", str(out / "last.pth")]
                       + SMALL_ARGS)
            resumed = torch.load(out / "last.pth", weights_only=False)
            self.assertEqual(resumed["epoch"], 3)
            self.assertEqual(resumed["stats"], ckpt["stats"])  # statistics are not refitted

            pred = tmp / "pred"
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(root), "--output", str(pred),
                          "--chroms", "chr1", "--num-workers", "0", "--amp", "off", "--batch-size", "5"])
            with gzip.open(pred / "sample.v3_tensors.SNV.predictions.tsv.gz", "rt") as f:
                lines = f.read().splitlines()
            chr1 = [t for t in truth["SNV"] if t["chrom"] == "chr1"]
            self.assertEqual([line.split("\t")[1] for line in lines[1:]], [t["candidate_id"] for t in chr1])
            self.assertTrue(all(abs(sum(map(float, line.split("\t")[3:6])) - 1) < 1e-3 for line in lines[1:]))
            report = json.loads((pred / "sample.v3_tensors.SNV.metrics.json").read_text())
            self.assertEqual(report["tensors"], sum(t["label"] >= 0 for t in chr1))
            self.assertEqual(report["ignored"], sum(t["label"] < 0 for t in chr1))
            self.assertEqual(report["threshold"], torch.load(out / "best.pth", weights_only=False)["somatic_threshold"])
            self.assertEqual(len(KindIndex(root / "SNV", pred / "index_cache")), len(truth["SNV"]))

    def test_scalars_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "v3_tensors"
            make_tensor_set(root, SPEC, shard_size=8)
            out = tmp / "run"
            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "1", "--kinds", "SNV", "--scalars"]
                       + SMALL_ARGS)
            ckpt = torch.load(out / "best.pth", weights_only=False)
            self.assertEqual(ckpt["config"]["scalars"], len(SCALARS))
            self.assertEqual(ckpt["scalars"], list(SCALARS))
            self.assertIn("scalar_mean", ckpt["model_state_dict"])
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(root), "--output", str(tmp / "p"),
                          "--kinds", "SNV", "--num-workers", "0", "--amp", "off"])
            self.assertTrue((tmp / "p" / "sample.v3_tensors.SNV.metrics.json").exists())


if __name__ == "__main__":
    unittest.main()

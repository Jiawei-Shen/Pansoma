import gzip
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

import numpy as np
import torch
import torch.nn as nn

from .. import predict, train
from ..data import SCALARS, EpochSampler, KindIndex, TensorDataset, block_split, load_parts
from ..model import PansomaNetV2, no_decay
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

    def test_dropped_planes_round_trip(self):
        x, blocks = batch(3)
        model = PansomaNetV2(3, stats=STATS, drop_planes=["path_count"], **SMALL).eval()
        self.assertEqual(model.front[0].in_channels, 35)
        payload = dict(format="pansoma_net_v2", config=model.config, model_state_dict=model.state_dict())
        again = PansomaNetV2.from_checkpoint(payload).eval()
        self.assertEqual(again.config["drop_planes"], ["path_count"])
        with torch.no_grad():
            self.assertTrue(torch.equal(again(x, blocks), model(x, blocks)))
            x2 = x.clone(); x2[:, 6] = 1
            self.assertTrue(torch.equal(model(x2, blocks), model(x, blocks)))
        old = dict(payload, config={k: v for k, v in model.config.items() if k != "drop_planes"})
        with self.assertRaises(RuntimeError):   # a config without drop_planes is a 36-plane model
            PansomaNetV2.from_checkpoint(old)

    def test_blocks_v1_and_v2(self):
        x, blocks = batch(2)
        v1, v2 = (PansomaNetV2(3, stats=STATS, block=b, **SMALL).eval() for b in ("v1", "v2"))
        self.assertEqual(v2.backbone.stages[0][0].grn.gamma.shape, (1, 1, 1, 64))    # per channel (4 x 16)
        self.assertEqual(v1.backbone.stages[0][0].grn.gamma.shape, (1,))
        # a v2 block with a zero branch is the identity; a v1 block applies GELU to the sum
        for model, expect in ((v2, lambda h: h), (v1, torch.nn.functional.gelu)):
            blk = model.backbone.stages[0][0]
            nn.init.zeros_(blk.pwconv2.weight), nn.init.zeros_(blk.pwconv2.bias)
            h = torch.randn(2, 16, 5, 7)
            with torch.no_grad():
                self.assertTrue(torch.allclose(blk(h), expect(h), atol=1e-6))
        # GRN v2: per-channel norm over the positions, divided by its mean over the channels
        g = PansomaNetV2(3, stats=STATS, **SMALL).backbone.stages[0][0].grn
        with torch.no_grad():
            g.gamma.fill_(1.0)
            h = torch.randn(2, 5, 7, 64)
            gx = h.pow(2).sum((1, 2), keepdim=True).sqrt()
            self.assertTrue(torch.allclose(g(h), h * gx / (gx.mean(-1, keepdim=True) + 1e-6) + h, atol=1e-5))
        # checkpoints from before "block" load as v1, and give the same logits
        payload = dict(format="pansoma_net_v2", config={k: v for k, v in v1.config.items() if k != "block"},
                       model_state_dict=v1.state_dict())
        again = PansomaNetV2.from_checkpoint(payload).eval()
        self.assertEqual(again.config["block"], "v1")
        with torch.no_grad():
            self.assertTrue(torch.equal(again(x, blocks), v1(x, blocks)))
        with self.assertRaises(RuntimeError):                                   # v1 weights do not fit v2
            PansomaNetV2(3, stats=STATS, block="v2", **SMALL).load_state_dict(v1.state_dict())

    def test_weight_decay_groups(self):
        model = PansomaNetV2(3, stats=STATS, **SMALL)
        params = dict(model.named_parameters())
        names = {n for n, p in params.items() if no_decay(n, p)}
        self.assertTrue(all(params[n].ndim <= 1 or ".grn." in n for n in names))   # biases, norm weights, GRN
        self.assertTrue(any(".grn.gamma" in n for n in names) and any(".norm.weight" in n for n in names))
        self.assertFalse(any(n.endswith(("pwconv1.weight", "dwconv.weight", "head.weight")) for n in names))

    def test_ignored_labels_do_not_enter_the_loss(self):
        logits = torch.randn(5, 3)
        y = torch.tensor([0, 1, -1, 2, -1])
        criterion = nn.CrossEntropyLoss(weight=torch.tensor([1.0, 5.0, 2.0]), ignore_index=-1)
        keep = y >= 0
        self.assertTrue(torch.allclose(criterion(logits, y), criterion(logits[keep], y[keep])))


class ClassWeightsTest(unittest.TestCase):
    def test_power_of_balanced(self):
        counts = np.array([900, 10, 90])
        balanced = train.class_weights("balanced", counts).numpy()
        np.testing.assert_allclose(balanced, 1000 / (3 * counts), rtol=1e-6)
        np.testing.assert_allclose(train.class_weights("sqrt", counts).numpy(), np.sqrt(balanced), rtol=1e-6)
        np.testing.assert_allclose(train.class_weights("pow:0.5", counts).numpy(), np.sqrt(balanced), rtol=1e-6)
        np.testing.assert_allclose(train.class_weights("pow:0.75", counts).numpy(), balanced ** 0.75, rtol=1e-6)
        with self.assertRaisesRegex(SystemExit, "pow:E needs a number"):
            train.class_weights("pow:x", counts)


class MicroBatchTest(unittest.TestCase):
    def test_micro_batches_give_the_whole_batch_update(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8)
            args = ["--tensors", str(root), "--epochs", "1", "--kinds", "SNV", "--seed", "3", "--drop-path", "0"] + SMALL_ARGS
            args[args.index("--batch-size") + 1] = "8"
            train.main(args + ["--output", str(tmp / "one")])
            train.main(args + ["--output", str(tmp / "two"), "--micro-batches", "2"])
            a = torch.load(tmp / "one" / "last.pth", weights_only=False)["model_state_dict"]
            b = torch.load(tmp / "two" / "last.pth", weights_only=False)["model_state_dict"]
            for k in a:
                self.assertTrue(torch.allclose(a[k].float(), b[k].float(), atol=1e-5, rtol=1e-4), k)
            rows = [json.loads(line) for line in (tmp / "one" / "metrics.jsonl").read_text().splitlines()]
            rows2 = [json.loads(line) for line in (tmp / "two" / "metrics.jsonl").read_text().splitlines()]
            self.assertAlmostEqual(rows[0]["train_loss"], rows2[0]["train_loss"], places=5)


class LoaderTest(unittest.TestCase):
    def test_forkserver_workers_read_the_same_tensors(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "s" / "t"
            make_tensor_set(root, {"SNV": {"chr1": 10}}, shard_size=3)
            ds = TensorDataset(load_parts([root], ["SNV"], Path(tmp) / "c", labelled=False))
            loader = train.make_loader(ds, EpochSampler(len(ds), shuffle=False), 2, 4, persistent=False)
            got = torch.cat([batch[0] for batch in loader])
            self.assertTrue(torch.equal(got, torch.stack([ds[k][0] for k in range(len(ds))])))


class RetryTest(unittest.TestCase):
    def test_worker_start_is_retried_then_raised(self):
        calls, logged = [], []

        def flaky():
            calls.append(1)
            if len(calls) < 3:
                raise RuntimeError("DataLoader worker (pid(s) 1) exited unexpectedly")
            return "ok"
        self.assertEqual(train.retry_workers(flaky, logged.append, "test", wait=0), "ok")
        self.assertEqual((len(calls), len(logged)), (3, 2))

        def gone():
            raise FileNotFoundError(2, "gone")
        with self.assertRaises(FileNotFoundError):
            train.retry_workers(gone, logged.append, "t", wait=0)


class ExcludeTest(unittest.TestCase):
    """--exclude-candidates: those tensors are -1 in training, validation and the test, as if never built."""

    def test_exclude_in_place(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "sample" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8)
            index = KindIndex(root / "SNV", Path(tmp) / "cache")
            a = {k: v.copy() for k, v in index.arrays.items()}
            cands = index.candidates()
            picked = np.flatnonzero(a["eval_label"] >= 0)[::3]
            self.assertTrue(len(picked))
            self.assertEqual(train.exclude(index, {cands[k] for k in picked} | {"1:0:SNP:A>C"}), len(picked))
            label, evl = index.arrays["label"], index.arrays["eval_label"]
            self.assertTrue((label[picked] == -1).all() and (evl[picked] == -1).all())
            other = np.setdiff1d(np.arange(len(label)), picked)
            self.assertTrue((label[other] == a["label"][other]).all() and (evl[other] == a["eval_label"][other]).all())

    def test_excluded_candidates_leave_training_validation_and_test(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8)
            index = KindIndex(root / "SNV", tmp / "check_cache")
            cands, a = index.candidates(), index.arrays
            tr, va = block_split(index, index.select(["chr2", "chr3"], labelled=False), 0.5, 21, 0)
            dropped = [int(k) for k in tr if a["label"][k] >= 0][:3]
            val_drop = [int(k) for k in va if a["eval_label"][k] >= 0][:2]
            chr1 = index.meta["chroms"].index("chr1")
            test_drop = [k for k in range(len(cands)) if a["chrom"][k] == chr1 and a["eval_label"][k] >= 0][:2]
            self.assertTrue(dropped and val_drop and test_drop)
            ids = [cands[k] for k in dropped + val_drop + test_drop]
            (tmp / "drop.tsv").write_text("candidate_id\taf\n" + "".join(f"{c}\t0.09\n" for c in ids))
            out = tmp / "run"
            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "1", "--kinds", "SNV",
                        "--exclude-candidates", str(tmp / "drop.tsv")] + SMALL_ARGS)
            ckpt = torch.load(out / "best.pth", weights_only=False)
            self.assertEqual(sum(d["train"] for d in ckpt["data"]), int((a["label"][tr] >= 0).sum()) - len(dropped))
            self.assertEqual(ckpt["val"]["tensors"], int((a["eval_label"][va] >= 0).sum()) - len(val_drop))
            log = (out / "train.log").read_text()
            self.assertIn(f"--exclude-candidates {tmp / 'drop.tsv'}: {len(ids)} candidate ids", log)
            self.assertIn(f"{index.dir}: {len(ids)} tensors excluded", log)
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(root), "--output", str(tmp / "p"),
                          "--kinds", "SNV", "--chroms", "chr1", "--num-workers", "0", "--amp", "off"])
            with gzip.open(tmp / "p" / "sample.tensors.SNV.predictions.ndjson.gz", "rt") as f:
                recs = {r["candidate_id"]: r for r in map(json.loads, f)}
            for k in test_drop:   # the checkpoint's list is predict's default
                self.assertTrue(recs[cands[k]]["excluded"] and not recs[cands[k]]["in_test"])
            self.assertEqual(sum(r["excluded"] for r in recs.values()), len(test_drop))
            rep = json.loads((tmp / "p" / "sample.tensors.SNV.metrics.json").read_text())
            self.assertEqual((rep["excluded"], rep["exclude_candidates"]), (len(test_drop), str(tmp / "drop.tsv")))
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(root), "--output", str(tmp / "p0"),
                          "--kinds", "SNV", "--chroms", "chr1", "--num-workers", "0", "--amp", "off",
                          "--exclude-candidates"])   # without a file: none
            rep0 = json.loads((tmp / "p0" / "sample.tensors.SNV.metrics.json").read_text())
            self.assertEqual((rep0["excluded"], rep0["tensors"]), (0, rep["tensors"] + len(test_drop)))


class SetTestChromsTest(unittest.TestCase):
    """--set-test-chroms: a chromosome left out of one set only; the other sets train and validate on theirs."""

    def test_one_set_leaves_out_its_chr1(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sets = {name: tmp / name / "tensors" for name in ("held", "other")}
            for seed, root in enumerate(sets.values()):
                make_tensor_set(root, SPEC, shard_size=8, seed=seed)
            expected = {}
            for name, chroms in (("held", ["chr2", "chr3"]), ("other", ["chr1", "chr2", "chr3"])):
                index = KindIndex(sets[name] / "SNV", tmp / "check_cache")
                tr, va = block_split(index, index.select(chroms, labelled=False), 0.5, 21, 0)
                expected[str(index.dir)] = (int((index.arrays["label"][tr] >= 0).sum()),
                                            int((index.arrays["eval_label"][va] >= 0).sum()))
            out = tmp / "run"
            train.main(["--tensors", *map(str, sets.values()), "--output", str(out), "--epochs", "1", "--kinds", "SNV",
                        "--test-chroms", "--set-test-chroms", str(sets["held"]), "chr1"] + SMALL_ARGS)
            ckpt = torch.load(out / "best.pth", weights_only=False)
            self.assertEqual({d["directory"]: d["train"] for d in ckpt["data"]},
                             {k: v[0] for k, v in expected.items()})
            self.assertEqual(ckpt["val"]["tensors"], sum(v[1] for v in expected.values()))
            self.assertEqual(ckpt["chroms"]["test"], [])
            self.assertEqual(ckpt["chroms"]["set_test"], {str(sets["held"].resolve()): ["chr1"]})
            with self.assertRaisesRegex(SystemExit, "not one of --tensors"):
                train.main(["--tensors", str(sets["other"]), "--output", str(tmp / "bad"), "--epochs", "1",
                            "--kinds", "SNV", "--set-test-chroms", str(sets["held"]), "chr1"] + SMALL_ARGS)


class PartialTestTest(unittest.TestCase):
    """predict's default: the reasons trained without stay out of the SNV test; INDEL keeps the partial somatic matches."""

    def test_indel_test_keeps_partial_somatic(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8)
            out = tmp / "run"
            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "1",
                        "--ignore-reasons", "residual_partial_somatic_truth", "allele_partial_somatic_truth"] + SMALL_ARGS)
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(root), "--output", str(tmp / "p"),
                          "--num-workers", "0", "--amp", "off"])
            recs, reps = {}, {}
            for kind in ("SNV", "INDEL"):
                with gzip.open(tmp / "p" / f"sample.tensors.{kind}.predictions.ndjson.gz", "rt") as f:
                    recs[kind] = [r for r in map(json.loads, f) if r["reason"] == "residual_partial_somatic_truth"]
                reps[kind] = json.loads((tmp / "p" / f"sample.tensors.{kind}.metrics.json").read_text())
                self.assertTrue(recs[kind])
            self.assertTrue(all(r["ignored"] and not r["in_test"] for r in recs["SNV"]))
            self.assertTrue(all(not r["ignored"] and r["in_test"] and r["test_label"] == 1 for r in recs["INDEL"]))
            self.assertEqual(reps["SNV"]["ignored_reasons"], ["residual_partial_somatic_truth", "allele_partial_somatic_truth"])
            self.assertEqual((reps["INDEL"]["ignored_reasons"], reps["INDEL"]["ignored"]), ([], 0))
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(root), "--output", str(tmp / "p2"),
                          "--kinds", "INDEL", "--num-workers", "0", "--amp", "off",
                          "--ignore-reasons", "residual_partial_somatic_truth"])   # an explicit list still applies
            rep = json.loads((tmp / "p2" / "sample.tensors.INDEL.metrics.json").read_text())
            self.assertEqual(rep["ignored"], len(recs["INDEL"]))


class TrainPredictTest(unittest.TestCase):
    """Small end-to-end runs on CPU: train 2 epochs (chr1 left out, validation from node blocks of chr2/chr3),
    resume to 3, predict chr1; the channel-6 storage checks; runs with the scalars and with --ignore-reasons."""

    def test_train_resume_predict(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "tensors"
            truth = make_tensor_set(root, SPEC, shard_size=8, seed=6)  # seed 6: chr1 has SNV tensors of an INS truth
            out = tmp / "run"
            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "2"] + SMALL_ARGS)
            for f in ("best.pth", "last.pth", "stats.json", "metrics.jsonl", "train.log"):
                self.assertTrue((out / f).exists(), f)
            rows = [json.loads(line) for line in (out / "metrics.jsonl").read_text().splitlines()]
            self.assertEqual([r["epoch"] for r in rows], [1, 2])
            ckpt = torch.load(out / "last.pth", weights_only=False)
            self.assertEqual(ckpt["chroms"]["test"], ["chr1"])
            self.assertIsNotNone(ckpt["somatic_threshold"])
            expected_train = expected_val = 0  # the same block split: training labels / evaluation labels
            for kind in ("SNV", "INDEL"):
                index = KindIndex(root / kind, tmp / "check_cache")
                tr, va = block_split(index, index.select(["chr2", "chr3"], labelled=False), 0.5, 21, 0)
                expected_train += int((index.arrays["label"][tr] >= 0).sum())
                expected_val += int((index.arrays["eval_label"][va] >= 0).sum())
            self.assertEqual(sum(d["train"] for d in ckpt["data"]), expected_train)  # chr1 and -1 never trained on
            self.assertEqual(rows[-1]["val"]["tensors"], expected_val)                # off-reference no-match as non
            self.assertGreater(expected_val, 0)
            self.assertIn("somatic_ap", rows[-1]["val"])
            self.assertGreater(rows[-1]["val"]["truth"]["truth_alleles"], 0)   # truth-level validation each epoch

            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "3", "--resume", str(out / "last.pth")]
                       + SMALL_ARGS)
            resumed = torch.load(out / "last.pth", weights_only=False)
            self.assertEqual(resumed["epoch"], 3)
            self.assertEqual((resumed["path_count_storage"], resumed["path_count_haplotypes"]),
                             ("int8-count-linear100-log2", []))
            self.assertEqual(resumed["stats"], ckpt["stats"])  # statistics are not refitted
            # norm dataset (the default): one BQ statistics row per training dataset, fitted on tensors of all its
            # chromosomes (chr1 too, labels not used), kept on --resume
            self.assertEqual((ckpt["config"]["norm"], ckpt["config"]["dataset_norm"]), ("dataset", ["base_quality"]))
            self.assertEqual([d["directory"] for d in ckpt["dataset_stats"]], [str(root / "SNV"), str(root / "INDEL")])
            self.assertTrue(all(d["stats"]["tensors"] > 0 for d in ckpt["dataset_stats"]))
            self.assertEqual(resumed["dataset_stats"], ckpt["dataset_stats"])
            snv_index = KindIndex(root / "SNV", tmp / "check_cache")
            snv_row = ckpt["dataset_stats"][0]["stats"]
            self.assertEqual(snv_row, train.dataset_stats(snv_index, snv_index.select(None, labelled=False), 10, 12345, 0,
                                                          lambda m: None, "check"))

            pred = tmp / "pred"
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(root), "--output", str(pred),
                          "--chroms", "chr1", "--num-workers", "0", "--amp", "off", "--batch-size", "5"])
            with gzip.open(pred / "sample.tensors.SNV.predictions.ndjson.gz", "rt") as f:
                records = [json.loads(line) for line in f]
            chr1 = [t for t in truth["SNV"] if t["chrom"] == "chr1"]
            self.assertEqual([r["candidate_id"] for r in records], [t["candidate_id"] for t in chr1])
            self.assertTrue(all(abs(r["p_non"] + r["p_somatic"] + r["p_germline"] - 1) < 1e-3 for r in records))
            self.assertEqual([r["in_test"] for r in records], [t["eval_label"] >= 0 for t in chr1])
            self.assertEqual([r["off_reference"] for r in records], [t["off_reference"] for t in chr1])
            self.assertEqual([r["reason"] for r in records], [t["reason"] for t in chr1])
            self.assertTrue(all(r["test_label"] == 0 for r in records if r["reason"] == "off_reference_no_truth_match"))
            report = json.loads((pred / "sample.tensors.SNV.metrics.json").read_text())
            self.assertEqual(report["tensors"], sum(t["eval_label"] >= 0 for t in chr1))
            # a training dataset: predict takes the statistics the model was trained with from the checkpoint
            norm = report["normalization"]
            self.assertEqual((norm["norm"], norm["channels"], norm["source"]), ("dataset", ["base_quality"], "checkpoint"))
            self.assertEqual(norm["dataset_stats"], snv_row)
            self.assertEqual(norm["shift"], {"base_quality": 0.0, "mapping_quality": 0.0, "path_count": 0.0})
            self.assertEqual(norm["checkpoint_stats"], torch.load(out / "best.pth", weights_only=False)["stats"])
            # any other dataset (the same tensors at another path): fitted on all its chromosomes, not on --chroms
            copy = tmp / "copy" / "tensors"
            make_tensor_set(copy, SPEC, shard_size=8, seed=6)
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(copy), "--output", str(tmp / "pc"),
                          "--chroms", "chr1", "--kinds", "SNV", "--num-workers", "0", "--amp", "off", "--batch-size", "5"])
            norm2 = json.loads((tmp / "pc" / "copy.tensors.SNV.metrics.json").read_text())["normalization"]
            self.assertEqual((norm2["source"], norm2["dataset_stats"]["tensors"]), ("fitted", len(truth["SNV"])))
            self.assertEqual(set(norm2["shift"]), {"base_quality", "mapping_quality", "path_count"})
            self.assertEqual(report["left_out"], sum(t["eval_label"] < 0 for t in chr1))
            self.assertEqual(report["off_reference_in_test"], sum(t["off_reference"] and t["eval_label"] >= 0 for t in chr1))
            self.assertEqual(report["threshold"], torch.load(out / "best.pth", weights_only=False)["somatic_threshold"])
            # three decision rules: threshold t (= the report above), argmax, the validation recall-0.9 threshold
            rules = report["rules"]
            self.assertEqual(set(rules), {"threshold_t", "argmax", "val_recall_0.9", "pipeline"})
            self.assertEqual({k: rules["threshold_t"]["per_tensor"][k] for k in ("precision", "recall", "f1", "support")},
                             {k: report["thresholded"]["somatic"][k] for k in ("precision", "recall", "f1", "support")})
            self.assertEqual(rules["argmax"]["per_tensor"]["f1"], report["argmax"]["somatic"]["f1"])
            for k in ("tp", "fp", "precision", "recall", "f1"):
                self.assertEqual(rules["threshold_t"]["truth"][k], report["truth"]["at_threshold"][k])
            best = torch.load(out / "best.pth", weights_only=False)
            self.assertEqual(rules["val_recall_0.9"]["threshold"], best["val"]["somatic_at_recall"]["0.9"]["threshold"])
            argmax_calls = sum(r["in_test"] and max(r["p_non"], r["p_somatic"], r["p_germline"]) == r["p_somatic"]
                               for r in records)
            self.assertEqual(rules["argmax"]["per_tensor"]["calls"], argmax_calls)
            # the SNV off-reference rescue: only off-reference tensors are rescored, and called by argmax
            off = [r for r in records if r["off_reference"]]
            self.assertTrue(off)
            self.assertTrue(all(("offref_call" in r) == r["off_reference"] for r in records))
            self.assertTrue(all(r["offref_call"] == (max(r["p_offref_non"], r["p_offref_somatic"], r["p_offref_germline"])
                                                     == r["p_offref_somatic"]) for r in off))
            o = report["offref_rescue"]
            self.assertEqual((o["tensors"], o["calls"], o["site_halfwidth"], o["floor"]),
                             (len(off), sum(r["offref_call"] for r in off), 10, 90))
            t09 = best["val"]["somatic_at_recall"]["0.9"]["threshold"]
            self.assertEqual(rules["pipeline"]["threshold"], t09)
            self.assertEqual(rules["pipeline"]["per_tensor"]["calls"],
                             sum(r["in_test"] and (r["offref_call"] if r["off_reference"] else r["p_somatic"] >= t09)
                                 for r in records))
            index = KindIndex(root / "SNV", tmp / "check_cache")
            k = next(i for i, r in enumerate(records) if r["off_reference"])
            chr1_pos = index.select(["chr1"], labelled=False)
            plain = TensorDataset([(index, chr1_pos[k:k + 1])])[0][0].numpy()
            edited = predict.OffrefRescue([(index, chr1_pos[k:k + 1])], 10, 90)[0][0].numpy()
            window = np.s_[predict.PATH_COUNT, :, predict.SITE - 10:predict.SITE + 11]
            expected = plain.copy()
            expected[window] = np.where((plain[window] > 0) & (plain[window] < 90), 90, plain[window])
            np.testing.assert_array_equal(edited, expected)  # ch6 raised at the site only, uncovered cells stay 0
            # against the truth table, recomputed from the predictions: every truth allele once
            tr = report["truth"]
            table = [line.split("\t") for line in (root / "somatic.recall.tsv").read_text().splitlines()[1:]]
            wanted = {int(r[0]) for r in table if r[1] == "chr1" and r[3] == "SNP" and r[5] == "True"}
            other = {int(r[0]) for r in table if r[1] == "chr1" and r[3] in ("DEL", "INS") and r[5] == "True"}
            by_id = {t["candidate_id"]: t for t in chr1}
            self.assertTrue(all(r["truth_ids"] == ([by_id[r["candidate_id"]]["truth_id"]] if r["label"] == 1 else [])
                                for r in records))
            called = [r for r in records if r["in_test"] and r["pred"] == "somatic"]
            hit = set().union(*[set(r["truth_ids"]) for r in called if r["test_label"] == 1])
            found, found_other = hit & wanted, hit & other           # an INS truth found by an SNV tensor: a true call
            false_calls = sum(r["test_label"] in (0, 2) for r in called)
            reached = set().union(*[set(r["truth_ids"]) for r in records if r["in_test"] and r["test_label"] == 1])
            at = tr["at_threshold"]
            self.assertEqual((tr["truth_alleles"], at["tp"], at["other_tp"], at["fp"], tr["with_tensor"]),
                             (len(wanted), len(found), len(found_other), false_calls, len(reached & wanted)))
            self.assertGreater(tr["other_kind_truth_with_tensor"], 0)
            self.assertAlmostEqual(at["recall"], len(found) / len(wanted))
            if hit or false_calls:
                self.assertAlmostEqual(at["precision"], len(hit & (wanted | other)) / (len(hit & (wanted | other)) + false_calls))
            self.assertLess(tr["ceiling"], 1.0)                      # truth alleles without a tensor are misses
            with gzip.open(pred / "sample.tensors.INDEL.predictions.ndjson.gz", "rt") as f:
                indel = [json.loads(line) for line in f]
            self.assertFalse(any("offref_call" in r for r in indel))  # the rescue is for SNV only
            self.assertEqual(len(KindIndex(root / "SNV", pred / "index_cache")), len(truth["SNV"]))

    def test_norm_train_and_checkpoints_without_norm(self):
        """--norm train: no per-dataset statistics, predict z-scores with the checkpoint's; a checkpoint from before
        "norm" is a norm train one; predict --norm dataset fits the predicted dataset's statistics for it."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8, seed=6)
            out = tmp / "run"
            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "1", "--kinds", "SNV", "--norm", "train"]
                       + SMALL_ARGS)
            ckpt = torch.load(out / "best.pth", weights_only=False)
            self.assertEqual((ckpt["config"]["norm"], ckpt["dataset_stats"]), ("train", []))
            self.assertFalse((out / "dataset_stats.json").exists())
            old = dict(ckpt, config={k: v for k, v in ckpt["config"].items() if k != "norm"})  # before 2026-10-13
            torch.save(old, tmp / "old.pth")
            self.assertEqual(PansomaNetV2.from_checkpoint(old).config["norm"], "train")
            self.assertEqual(PansomaNetV2.from_checkpoint(ckpt).config["dataset_norm"], [])
            # a norm dataset checkpoint of 2026-10-13 to 10-19 (no dataset_norm): BQ and MAPQ per dataset
            old_ds = dict(ckpt, config={**{k: v for k, v in ckpt["config"].items() if k != "dataset_norm"}, "norm": "dataset"})
            self.assertEqual(PansomaNetV2.from_checkpoint(old_ds).encoder.dataset_norm, ("base_quality", "mapping_quality"))
            for name, checkpoint, extra, want in (("p_train", out / "best.pth", [], None),
                                                  ("p_old", tmp / "old.pth", [], None),
                                                  ("p_dataset", tmp / "old.pth", ["--norm", "dataset"], "dataset")):
                predict.main(["--checkpoint", str(checkpoint), "--tensors", str(root), "--output", str(tmp / name),
                              "--chroms", "chr1", "--num-workers", "0", "--amp", "off", "--kinds", "SNV"] + extra)
                norm = json.loads((tmp / name / "sample.tensors.SNV.metrics.json").read_text())["normalization"]
                self.assertEqual(norm["norm"], want or "train")
                self.assertEqual(norm["dataset_stats"] is not None, want == "dataset")
                if want:  # a norm train checkpoint has no dataset statistics: fitted, BQ (encode.DATASET_NORM)
                    self.assertEqual((norm["channels"], norm["source"], norm["shift"]), (["base_quality"], "fitted", None))
            same = lambda a, b: [json.loads(l)["p_somatic"] for l in gzip.open(tmp / a / "sample.tensors.SNV.predictions.ndjson.gz", "rt")] == \
                [json.loads(l)["p_somatic"] for l in gzip.open(tmp / b / "sample.tensors.SNV.predictions.ndjson.gz", "rt")]  # noqa: E731
            self.assertTrue(same("p_train", "p_old"))  # an old checkpoint predicts as before

    def test_channel_6_storage_is_checked_and_sets_the_rescue_floor(self):
        """Training needs one channel-6 storage and records it; predict refuses another one before predicting any set
        (a checkpoint without the record counts as int8-count-linear100-log2), rescores off-reference SNVs with ch6
        raised to that storage's every-haplotype code (100 for int8-count-haplotypes100, 90 for the old sets) and notes
        another --haplotypes. A model trained without the path-count plane takes sets of any storage."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            new, old, d46 = (tmp / name / "tensors" for name in ("new", "old", "d46"))
            make_tensor_set(new, SPEC, shard_size=8, seed=6, storage="int8-count-haplotypes100", haplotypes=90)
            make_tensor_set(d46, SPEC, shard_size=8, seed=6, storage="int8-count-haplotypes100", haplotypes=464)
            # the deleted v1.1 sets spelled it tensor_storage_version int8-count-linear100-log2-v1, without haplotypes
            make_tensor_set(old, SPEC, shard_size=8, seed=6, storage="int8-count-linear100-log2-v1",
                            storage_key="tensor_storage_version")
            with self.assertRaisesRegex(SystemExit, "code channel 6 differently"):
                train.main(["--tensors", str(new), str(old), "--output", str(tmp / "mixed"), "--epochs", "1"] + SMALL_ARGS)
            out = tmp / "run"
            train.main(["--tensors", str(new), "--output", str(out), "--epochs", "1"] + SMALL_ARGS)
            ckpt = torch.load(out / "best.pth", weights_only=False)
            self.assertEqual((ckpt["path_count_storage"], ckpt["path_count_haplotypes"]), ("int8-count-haplotypes100", [90]))
            args = ["--chroms", "chr1", "--num-workers", "0", "--amp", "off", "--kinds", "SNV"]

            def run_predict(checkpoint, sets, output):
                """predict.main with the floors main() gives OffrefRescue recorded, and its stdout."""
                floors, real = [], predict.OffrefRescue
                rescue = lambda parts, halfwidth, floor: floors.append(floor) or real(parts, halfwidth, floor)  # noqa: E731
                with mock.patch.object(predict, "OffrefRescue", rescue), redirect_stdout(io.StringIO()) as printed:
                    predict.main(["--checkpoint", str(checkpoint), "--tensors", *map(str, sets), "--output", str(output)]
                                 + args)
                return floors, printed.getvalue()

            with self.assertRaisesRegex(SystemExit, "old/tensors/SNV: channel 6 is int8-count-linear100-log2, the "
                                                    "checkpoint was trained on int8-count-haplotypes100"):
                run_predict(out / "best.pth", [new, old], tmp / "p0")  # the new set comes first: still none predicted
            self.assertEqual(list((tmp / "p0").glob("*.predictions.ndjson.gz")), [])
            floors, printed = run_predict(out / "best.pth", [new], tmp / "p1")
            self.assertEqual(floors, [100])  # the rescue rescored with the every-haplotype code of the set's storage
            report = json.loads((tmp / "p1" / "new.tensors.SNV.metrics.json").read_text())
            self.assertEqual(report["offref_rescue"]["floor"], 100)
            self.assertEqual(report["path_count"], dict(storage="int8-count-haplotypes100", haplotypes=90,
                                                        checkpoint_haplotypes=[90]))
            self.assertNotIn("note:", printed)
            floors, printed = run_predict(out / "best.pth", [d46], tmp / "p2")  # same storage, another H
            self.assertEqual(floors, [100])
            self.assertIn("was built with --haplotypes 464, the checkpoint's training sets with [90]", printed)
            self.assertEqual(json.loads((tmp / "p2" / "d46.tensors.SNV.metrics.json").read_text())["path_count"],
                             dict(storage="int8-count-haplotypes100", haplotypes=464, checkpoint_haplotypes=[90]))
            index = KindIndex(new / "SNV", tmp / "cache")
            k = int(np.flatnonzero(index.arrays["off_reference"])[0])
            plain = TensorDataset([(index, np.array([k]))])[0][0].numpy()
            edited = predict.OffrefRescue([(index, np.array([k]))], 10, 100)[0][0].numpy()
            window = np.s_[predict.PATH_COUNT, :, predict.SITE - 10:predict.SITE + 11]
            np.testing.assert_array_equal(edited[window], np.where((plain[window] > 0) & (plain[window] < 100), 100,
                                                                   plain[window]))
            legacy = {k: v for k, v in ckpt.items() if not k.startswith("path_count")}  # trained before the record
            torch.save(legacy, tmp / "legacy.pth")
            with self.assertRaisesRegex(SystemExit, "trained on int8-count-linear100-log2"):
                run_predict(tmp / "legacy.pth", [new], tmp / "p3")
            floors, _ = run_predict(tmp / "legacy.pth", [old], tmp / "p4")
            self.assertEqual(floors, [90])
            self.assertEqual(json.loads((tmp / "p4" / "old.tensors.SNV.metrics.json").read_text())["path_count"],
                             dict(storage="int8-count-linear100-log2", haplotypes=None, checkpoint_haplotypes=[]))
            with self.assertRaisesRegex(SystemExit, "trained on channel 6 int8-count-haplotypes100"):
                train.main(["--tensors", str(old), "--output", str(tmp / "resume"), "--epochs", "2", "--resume",
                            str(out / "last.pth")] + SMALL_ARGS)
            # without the path-count plane the storage does not matter: mixed sets train, any set predicts
            nopc = tmp / "nopc"
            train.main(["--tensors", str(new), str(old), "--output", str(nopc), "--epochs", "1", "--drop-planes",
                        "path_count"] + SMALL_ARGS)
            ckpt = torch.load(nopc / "best.pth", weights_only=False)
            self.assertEqual((ckpt["path_count_storage"], ckpt["path_count_haplotypes"]), (None, [90]))
            floors, _ = run_predict(nopc / "best.pth", [new, old], tmp / "p5")
            self.assertEqual(floors, [100, 90])

    def test_scalars_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8)
            out = tmp / "run"
            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "1", "--kinds", "SNV", "--scalars"]
                       + SMALL_ARGS)
            ckpt = torch.load(out / "best.pth", weights_only=False)
            self.assertEqual(ckpt["config"]["scalars"], len(SCALARS))
            train.main(["--tensors", str(root), "--output", str(tmp / "run_keep"), "--epochs", "1", "--kinds", "SNV",
                        "--non-fraction", "0.5", "--keep-non-af", "0.5"] + SMALL_ARGS)
            log = (tmp / "run_keep" / "train.log").read_text()
            self.assertIn("importance weight 2.00", log)
            self.assertTrue((tmp / "run_keep" / "best.pth").exists())
            self.assertEqual(ckpt["scalars"], list(SCALARS))
            self.assertIn("scalar_mean", ckpt["model_state_dict"])
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(root), "--output", str(tmp / "p"),
                          "--kinds", "SNV", "--num-workers", "0", "--amp", "off"])
            self.assertTrue((tmp / "p" / "sample.tensors.SNV.metrics.json").exists())

    def test_ignore_reasons_run(self):
        """--ignore-reasons leaves those tensors out of training and validation, and predict (by default the
        checkpoint's reasons) out of the test, with label, reason and truth ids kept in the records."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8)
            train.main(["--tensors", str(root), "--output", str(tmp / "run_ignore"), "--epochs", "1", "--kinds", "SNV",
                        "--ignore-reasons", "residual_partial_somatic_truth"] + SMALL_ARGS)
            log = (tmp / "run_ignore" / "train.log").read_text()
            index = KindIndex(root / "SNV", tmp / "run_ignore" / "index_cache")
            tr, va = block_split(index, index.select(["chr2", "chr3"], labelled=False), 0.5, 21, 0)
            partial = index.meta["reasons"].index("residual_partial_somatic_truth")
            n_tr = int(((index.arrays["label"][tr] >= 0) & (index.arrays["reason"][tr] == partial)).sum())
            n_va = int(((index.arrays["eval_label"][va] >= 0) & (index.arrays["reason"][va] == partial)).sum())
            self.assertGreater(n_tr + n_va, 0)
            self.assertIn(f"left out (--ignore-reasons residual_partial_somatic_truth): {n_tr:,} training, {n_va:,} validation", log)
            predict.main(["--checkpoint", str(tmp / "run_ignore" / "best.pth"), "--tensors", str(root), "--output",
                          str(tmp / "p_ignore"), "--kinds", "SNV", "--num-workers", "0", "--amp", "off"])
            with gzip.open(tmp / "p_ignore" / "sample.tensors.SNV.predictions.ndjson.gz", "rt") as f:
                recs = [json.loads(line) for line in f]
            part = [r for r in recs if r["reason"] == "residual_partial_somatic_truth"]
            self.assertTrue(part and all(r["ignored"] and not r["in_test"] and r["test_label"] is None and r["label"] == 1
                                         and r["truth_ids"] for r in part))           # left out, information kept
            self.assertFalse(any(r["ignored"] for r in recs if r["reason"] != "residual_partial_somatic_truth"))
            rep = json.loads((tmp / "p_ignore" / "sample.tensors.SNV.metrics.json").read_text())
            self.assertEqual((rep["ignored"], rep["ignored_reasons"]), (len(part), ["residual_partial_somatic_truth"]))
            ckpt_i = torch.load(tmp / "run_ignore" / "best.pth", weights_only=False)
            self.assertEqual(sum(d["train"] for d in ckpt_i["data"]),
                             int((index.arrays["label"][tr] >= 0).sum()) - n_tr)


if __name__ == "__main__":
    unittest.main()

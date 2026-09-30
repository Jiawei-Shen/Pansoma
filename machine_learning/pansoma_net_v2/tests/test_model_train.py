import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

from .. import combine, predict, train
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


class RelabelTest(unittest.TestCase):
    """--reason-labels (inside the confident region) and --label-overrides (any tensor; exempt from --ignore-reasons)."""

    def test_relabel_in_place(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "sample" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8)
            index = KindIndex(root / "SNV", Path(tmp) / "cache")
            a = {k: v.copy() for k, v in index.arrays.items()}
            reasons, cands = index.meta["reasons"], index.candidates()
            off = np.flatnonzero(a["reason"] == reasons.index("off_reference_no_truth_match"))
            region = a["in_region"].astype(bool)
            self.assertTrue(off.size and region[off].any())
            partial = int(np.flatnonzero(a["reason"] == reasons.index("residual_partial_somatic_truth"))[0])
            overrides = {cands[partial]: -1, cands[int(off[0])]: 1}
            train.relabel([(index, None)], overrides, {"off_reference_no_truth_match": 0})
            label, evl = index.arrays["label"], index.arrays["eval_label"]
            rest = off[1:]
            self.assertTrue((label[rest[region[rest]]] == 0).all() and (evl[rest[region[rest]]] == 0).all())
            self.assertTrue((label[rest[~region[rest]]] == a["label"][rest[~region[rest]]]).all())  # outside: unchanged
            self.assertEqual((int(label[partial]), int(evl[partial])), (-1, -1))
            self.assertEqual(int(label[off[0]]), 1)
            self.assertEqual(int(evl[off[0]]), 1 if region[off[0]] else -1)
            other = np.setdiff1d(np.arange(len(label)), np.concatenate([off, [partial]]))
            self.assertTrue((label[other] == a["label"][other]).all() and (evl[other] == a["eval_label"][other]).all())

    def test_overrides_are_exempt_from_ignore_reasons(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "sample" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8)
            index = KindIndex(root / "SNV", tmp / "check_cache")
            reasons, cands = index.meta["reasons"], index.candidates()
            tr, _ = block_split(index, index.select(["chr2", "chr3"], labelled=False), 0.5, 21, 0)
            partial_code = reasons.index("residual_partial_somatic_truth")
            promoted = [int(k) for k in tr if index.arrays["reason"][k] == partial_code][:1]
            self.assertTrue(promoted)
            (tmp / "overrides.tsv").write_text("candidate_id\tlabel\n" + "".join(f"{cands[k]}\t1\n" for k in promoted))
            out = tmp / "run"
            train.main(["--tensors", str(root), "--output", str(out), "--epochs", "1", "--kinds", "SNV",
                        "--ignore-reasons", "residual_partial_somatic_truth", "--label-overrides", str(tmp / "overrides.tsv"),
                        "--reason-labels", "off_reference_no_truth_match=0"] + SMALL_ARGS)
            train.relabel([(index, None)], {cands[k]: 1 for k in promoted}, {"off_reference_no_truth_match": 0})
            a = index.arrays
            expected = int(((a["label"][tr] >= 0) & ((a["reason"][tr] != partial_code) | np.isin(tr, promoted))).sum())
            ckpt = torch.load(out / "best.pth", weights_only=False)
            self.assertEqual(sum(d["train"] for d in ckpt["data"]), expected)
            self.assertIn("--label-overrides 1 of 1 candidate ids", (out / "train.log").read_text())


class TrainPredictTest(unittest.TestCase):
    """Small end-to-end runs on CPU: train 2 epochs (chr1 left out, validation from node blocks of chr2/chr3),
    resume to 3, predict chr1; and one run with the scalars."""

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
            self.assertEqual(resumed["stats"], ckpt["stats"])  # statistics are not refitted

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
            edited = predict.OffrefSite90([(index, chr1_pos[k:k + 1])], 10)[0][0].numpy()
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
            # SNV and INDEL together against the whole truth VCF (combine): every truth allele once, whichever set
            both = combine.main([str(pred)])
            self.assertEqual(len(both), 1)
            c = both[0]
            with gzip.open(pred / "sample.tensors.INDEL.predictions.ndjson.gz", "rt") as f:
                indel = [json.loads(line) for line in f]
            self.assertFalse(any("offref_call" in r for r in indel))  # the rescue is for SNV only
            records += indel
            called = [r for r in records if r["in_test"] and r["pred"] == "somatic"]
            everything = wanted | other
            hit = set().union(*[set(r["truth_ids"]) for r in called if r["test_label"] == 1]) & everything
            false_calls = sum(r["test_label"] in (0, 2) for r in called)
            self.assertEqual((c["truth_alleles"], c["tp"], c["fp"], c["other_kind_tp"]),
                             (len(everything), len(hit), false_calls, 0))
            self.assertAlmostEqual(c["recall"], len(hit) / len(everything))
            self.assertEqual(set(c["per_truth_kind"]), {"SNP", "DEL", "INS"})
            on_truth = [r for r in called if r["test_label"] == 1 and set(r["truth_ids"]) & everything]
            self.assertEqual(c["repeated_calls"], len(on_truth) - len(hit))  # one truth id per tensor here
            self.assertEqual(len(KindIndex(root / "SNV", pred / "index_cache")), len(truth["SNV"]))

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
            train.main(["--tensors", str(root), "--output", str(tmp / "run_truth"), "--epochs", "2", "--kinds", "SNV",
                        "--select", "truth_f1"] + SMALL_ARGS)
            best = torch.load(tmp / "run_truth" / "best.pth", weights_only=False)
            self.assertEqual(best["somatic_threshold"], best["val"]["truth"]["best"]["threshold"])
            self.assertEqual(json.loads((tmp / "run_truth" / "args.json").read_text())["select"], "truth_f1")
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


if __name__ == "__main__":
    unittest.main()

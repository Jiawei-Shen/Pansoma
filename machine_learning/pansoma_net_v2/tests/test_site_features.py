"""Site features (site_features.py): loading and checks, the fixed transform, the dataset's extra scalars, the model's
site-feature input, the post-PoN validation report (metrics.postpon_report), train / predict --site-features, train
--best-by postpon and the call rule, --non-reasons / --non-reasons-in-bed, the feature files' sha256 in the checkpoint
and predict's warning, the builder, and graph_vcf --threshold postpon and its default call rule."""
import contextlib
import gzip
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset

from .. import graph_vcf, metrics, predict, site_features, train
from ..data import SCALARS, EpochSampler, KindIndex, TensorDataset, block_split
from ..model import PansomaNetV2
from .fixtures import GAP, make_site_features, make_tensor_set
from .test_encode import STATS, batch
from .test_model_train import SMALL, SMALL_ARGS, SPEC

F = len(site_features.SITE_FEATURES)


class LoadTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        tmp = Path(self.tmp.name)
        self.root = tmp / "sample" / "tensors"
        make_tensor_set(self.root, SPEC, shard_size=8, seed=6)
        self.index = KindIndex(self.root / "SNV", tmp / "cache")  # the SNV fixture has off-reference tensors
        self.dir = tmp / "features"
        self.data = make_site_features(self.dir, self.index, seed=1)
        self.path = self.dir / "sample.SNV.npz"

    def tearDown(self):
        self.tmp.cleanup()

    def rewrite(self, **changes):
        np.savez(self.path, **{**self.data, **changes})

    def test_set_name_and_path(self):
        self.assertEqual(site_features.set_name(self.index), "sample")  # <set>/tensors/<KIND>
        self.assertEqual(site_features.feature_path(self.dir, self.index), self.path)
        self.assertTrue(self.path.exists())

    def test_loads_and_attaches(self):
        got = site_features.load(self.dir, self.index)
        for k in ("popaf_exact", "popaf_pos", "pon_tag", "hp_len"):
            np.testing.assert_array_equal(got[k], self.data[k])
        self.assertEqual(got["sha256"], hashlib.sha256(self.path.read_bytes()).hexdigest())
        self.assertEqual(str(got["candidates_sha1"]), site_features.candidates_sha1(self.index))
        fresh = KindIndex(self.root / "SNV", Path(self.tmp.name) / "cache")  # from the cache: ids not read yet
        site_features.load(self.dir, fresh)
        self.assertIsNone(fresh._candidates)  # the check does not keep millions of ids in memory
        site_features.attach(self.index, self.dir, model_input=False)
        self.assertNotIn("site_features", self.index.arrays)  # the tags only
        np.testing.assert_array_equal(self.index.arrays["pon_tag"], self.data["pon_tag"])
        site_features.attach(self.index, self.dir)
        self.assertEqual(self.index.arrays["site_features"].shape, (len(self.index), F))

    def test_refuses_files_that_disagree_with_the_index(self):
        n = len(self.index)
        cases = [
            (dict(candidates_sha1=np.array("0" * 40)), "candidates_sha1"),
            (dict(chrom_counts=self.data["chrom_counts"][::-1].copy()), "chrom_counts"),
            (dict(chroms=np.array(self.index.meta["chroms"][::-1])), "chromosomes"),
            (dict(popaf_pos=self.data["popaf_pos"][:-1]), "popaf_pos has"),
            (dict(popaf_exact=np.where(np.arange(n) == 0, np.nan, self.data["popaf_exact"]).astype(np.float32)),
             "not finite"),
            (dict(popaf_exact=np.full(n, 2.0, np.float32)), "above 1"),
            (dict(meta=np.array(json.dumps(dict(pon_rule=dict(site_features.pon_rule(), indel_min_af=0.001))))),
             "PoN tags of the rule"),
        ]
        off = np.flatnonzero(self.index.arrays["off_reference"])
        self.assertGreater(len(off), 0)
        shifted = self.data["pon_tag"].copy()
        shifted[off] = 0  # a file in another order: an off-reference tensor with a tag
        cases.append((dict(pon_tag=shifted), "off-reference tensors have"))
        on = np.flatnonzero(~self.index.arrays["off_reference"])
        untagged = self.data["pon_tag"].copy()
        untagged[on[0]] = -1  # the other way: a tensor with a GRCh38 allele and no tag
        cases.append((dict(pon_tag=untagged), "tensors with a GRCh38 allele"))
        no_af = self.data["popaf_exact"].copy()
        no_af[on[1]] = -1
        cases.append((dict(popaf_exact=no_af), "tensors with a GRCh38 allele"))
        cases.append((dict(vcf_pos=self.data["vcf_pos"] + 100000), "vcf_pos agrees"))
        for change, message in cases:
            self.rewrite(**change)
            with self.assertRaisesRegex(ValueError, message):
                site_features.load(self.dir, self.index)
        for key in ("hp_len", "candidates_sha1"):
            np.savez(self.path, **{k: v for k, v in self.data.items() if k != key})
            with self.assertRaisesRegex(ValueError, f"missing {key}"):
                site_features.load(self.dir, self.index)
        self.path.unlink()
        with self.assertRaisesRegex(ValueError, "no site features"):
            site_features.load(self.dir, self.index)

    def test_a_relabel_keeps_the_file_valid(self):
        manifest = self.root / "SNV" / "labels.manifest.json"
        manifest.write_text(manifest.read_text() + " ")
        site_features.load(self.dir, self.index)  # the features do not depend on the labels
        self.rewrite(labels_sha=np.array("0" * 64))  # older files carry it: not checked
        site_features.load(self.dir, self.index)

    def test_pon_rule_and_script(self):
        current = hashlib.sha256(site_features.PON_SCRIPT.read_bytes()).hexdigest()
        self.rewrite(meta=np.array(json.dumps(dict(pon_rule=site_features.pon_rule(), pon_script_sha256=current))))
        data = site_features.load(self.dir, self.index)
        self.assertEqual(site_features.pon_check(data), [])
        self.rewrite(meta=np.array(json.dumps(dict(pon_script_sha256="1" * 64))))  # no rule recorded: sha256 only
        data = site_features.load(self.dir, self.index)
        self.assertRegex(site_features.pon_check(data)[0], r"built with .* sha256 111111111111.*\(rule not recorded\)")
        rule = site_features.pon_rule()
        self.assertEqual((rule["indel_min_af"], rule["indel_min_af_homopolymer"], rule["indel_homopolymer"]),
                         (site_features.pon_module().INDEL_MIN_AF, 0.01, 8))

    def test_transform(self):
        data = dict(popaf_exact=np.array([0, 1e-4, 1, -1, 0.01], np.float32),
                    popaf_pos=np.array([1e-6, 0, 0.5, -1, 0.01], np.float32),
                    pon_tag=np.array([0, 1, 1, -1, 0], np.int8), hp_len=np.array([0, 10, 30, 0, 8], np.int16))
        t = site_features.transform(data)
        self.assertEqual((t.shape, t.dtype), ((5, F), np.float32))
        self.assertTrue(((t >= 0) & (t <= 1)).all())
        col = {name: t[:, k] for k, name in enumerate(site_features.SITE_FEATURES)}
        np.testing.assert_allclose(col["log_popaf_exact"], [0, 1 + np.log10(1e-4 + 1e-6) / 6, 1, 0,
                                                            1 + np.log10(0.01 + 1e-6) / 6], atol=1e-6)
        self.assertAlmostEqual(float(col["log_popaf_pos"][0]), 1 + np.log10(2e-6) / 6, places=6)
        np.testing.assert_array_equal(col["pon_tag"], [0, 1, 1, 0, 0])
        np.testing.assert_allclose(col["hp_len"], [0, 0.5, 1, 0, 0.4])
        np.testing.assert_array_equal(col["has_grch38"], [1, 1, 1, 0, 1])

    def test_dataset_appends_the_features_to_the_scalars(self):
        positions = self.index.select(None, labelled=False)
        plain = TensorDataset([(self.index, positions)])
        self.assertEqual(plain.scalars.shape[1], len(SCALARS))
        site_features.attach(self.index, self.dir)
        ds = TensorDataset([(self.index, positions)])
        self.assertEqual(ds.scalars.shape, (len(positions), len(SCALARS) + F))
        np.testing.assert_array_equal(ds.scalars[:, :len(SCALARS)], plain.scalars)
        np.testing.assert_array_equal(ds.scalars[:, len(SCALARS):], site_features.transform(self.data)[positions])
        np.testing.assert_array_equal(ds[3][2].numpy(), ds.scalars[3])
        other = KindIndex(self.root / "INDEL", Path(self.tmp.name) / "cache")  # without features
        with self.assertRaisesRegex(ValueError, "site features for some tensor sets only"):
            TensorDataset([(self.index, positions), (other, other.select(None, labelled=False))])


class ModelTest(unittest.TestCase):
    def test_site_features_reach_the_logits_and_round_trip(self):
        x, blocks = batch(4)
        scalars = torch.cat([torch.randn(4, len(SCALARS)), torch.rand(4, F)], 1)
        for n in (0, len(SCALARS)):
            model = PansomaNetV2(3, stats=STATS, scalars=n, site_features=F, **SMALL).eval()
            self.assertEqual((model.config["site_features"], model.scalar_mlp[0].in_features), (F, n + F))
            if n:
                model.set_scalar_stats(torch.zeros(n), torch.full((n,), 2.0))
            payload = dict(format="pansoma_net_v2", config=model.config, model_state_dict=model.state_dict())
            again = PansomaNetV2.from_checkpoint(payload).eval()
            with torch.no_grad():
                logits = model(x, blocks, scalars)
                self.assertTrue(torch.equal(again(x, blocks, scalars), logits))
                changed = scalars.clone()
                changed[:, -F:] += 1  # the site features reach the logits
                self.assertFalse(torch.equal(model(x, blocks, changed), logits))
                if not n:  # a model without --scalars ignores the SCALARS columns
                    changed = scalars.clone()
                    changed[:, :len(SCALARS)] += 1
                    self.assertTrue(torch.equal(model(x, blocks, changed), logits))
                with self.assertRaisesRegex(ValueError, "site features"):
                    model(x, blocks, scalars[:, :len(SCALARS)])
        plain = PansomaNetV2(3, stats=STATS, scalars=len(SCALARS), **SMALL)
        self.assertNotIn("site_features", plain.config)  # configs without site features do not change
        self.assertEqual(plain.scalar_mlp[0].in_features, len(SCALARS))


def report_of(labels, probs, pon):
    return metrics.postpon_report(np.asarray(labels), np.asarray(probs, np.float64), np.asarray(pon))


def probs_of(p_somatic):
    p = np.asarray(p_somatic, np.float64)
    return np.stack([1 - p, p, np.zeros_like(p)], 1)


class PostPonTest(unittest.TestCase):
    def test_without_tags_it_is_the_plain_threshold(self):
        rng = np.random.default_rng(3)
        labels = rng.choice([-1, 0, 1, 2], 400)
        probs = rng.dirichlet([1, 1, 1], 400)
        for pon in (np.zeros(400, np.int8), np.where(rng.random(400) < 0.2, -1, 0)):  # -1: callable too
            r = report_of(labels, probs, pon)
            keep = labels >= 0
            t, p, rec, f1 = metrics.best_somatic_threshold(labels[keep], probs[keep])
            self.assertEqual((r["threshold"], r["precision"], r["recall"], r["f1"]), (t, p, rec, f1))
            self.assertEqual(r["at_recall"], metrics.somatic_at_recall(labels[keep], probs[keep]))
            self.assertEqual((r["tagged"], r["ceiling"]), (0, 1.0))
            self.assertEqual({k: v and {f: x for f, x in v.items() if f != "recall_callable"}
                              for k, v in r["at_recall_callable"].items()}, r["at_recall"])  # nothing tagged

    def test_tagged_tensors_are_never_called(self):
        # p:       0.95 0.9 0.8 0.7 0.6 0.5 0.4 0.3
        labels = [0, 1, 0, 1, 1, 0, 1, -1]
        p = [0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
        plain = report_of(labels, probs_of(p), [0] * 8)
        self.assertAlmostEqual(plain["threshold"], 0.4)  # every one: P 4/7, R 1, F1 0.727
        r = report_of(labels, probs_of(p), [1, 0, 1, 0, 0, 0, 1, 1])  # the two top negatives and one positive tagged
        self.assertEqual((r["threshold"], r["calls"], r["tagged"]), (0.6, 3, 3))  # tensor -1 is not counted
        self.assertAlmostEqual(r["precision"], 1.0)
        self.assertAlmostEqual(r["recall"], 0.75)  # the tagged positive is a miss
        self.assertAlmostEqual(r["ceiling"], 0.75)
        self.assertEqual(r["tagged_by_label"], {"non": 2, "somatic": 1, "germline": 0})
        self.assertEqual(r["at_recall"]["0.5"]["threshold"], 0.7)
        self.assertIsNone(r["at_recall"]["0.8"])  # out of reach: one of four positives is tagged
        # relative to the three callable positives: 0.8 needs 3 of them (p 0.6), 0.5 two (p 0.7)
        self.assertEqual(r["callable_positives"], 3)
        c = r["at_recall_callable"]
        self.assertEqual((c["0.5"]["threshold"], c["0.8"]["threshold"], c["0.95"]["threshold"]), (0.7, 0.6, 0.6))
        self.assertEqual((c["0.8"]["calls"], c["0.8"]["recall_callable"]), (3, 1.0))
        self.assertAlmostEqual(c["0.8"]["recall"], 0.75)  # recall and F1 stay over every positive
        self.assertAlmostEqual(c["0.8"]["f1"], 2 * 0.75 / 1.75)
        r = report_of(labels, probs_of(p), [0, 1, 0, 1, 1, 0, 1, 0])  # every positive tagged
        self.assertEqual((r["recall"], r["f1"], r["ceiling"]), (0.0, 0.0, 0.0))
        self.assertEqual(r["at_recall_callable"], {k: None for k in ("0.5", "0.8", "0.9", "0.95")})
        self.assertEqual(report_of([0, 0], probs_of([0.2, 0.3]), [0, 0])["threshold"], 0.5)  # no positives

    def test_truth_level(self):
        # tensors: two of truth 7 (one tagged), one of truth 8 (tagged), negatives (one tagged with the top score)
        labels = np.array([1, 1, 1, 0, 0, 2])
        p = probs_of([0.9, 0.6, 0.8, 0.95, 0.5, 0.7])
        pon = np.array([1, 0, 1, 1, 0, 0])
        matches = [{7}, {7}, {8}, set(), set(), set()]
        r = metrics.postpon_report(labels, p, pon, ({7, 8, 9}, matches, set()))
        tr = r["truth"]
        self.assertEqual((tr["truth_alleles"], tr["with_tensor"]), (3, 1))  # truth 7 by its untagged tensor only
        self.assertEqual(tr["negatives"], 2)  # the tagged negative is no false call
        expect = metrics.truth_report({7, 8, 9}, matches, np.array([-1, 1, -1, -1, 0, 2]), p, r["threshold"], set())
        self.assertEqual(tr, expect)
        self.assertEqual(tr["best"]["threshold"], 0.6)  # truth 7 at p 0.6, after the germline 0.7: P 1/2, R 1/3


class ConstantModel(torch.nn.Module):
    """Logits from the first scalar column (a test stand-in for the network in train.evaluate)."""

    def forward(self, x, blocks, scalars, part):
        s = scalars[:, 0]
        return torch.stack([torch.zeros_like(s), s, torch.zeros_like(s)], 1)


class Rows(Dataset):
    def __init__(self, score, label):
        self.score, self.label = score, label

    def __len__(self):
        return len(self.label)

    def __getitem__(self, i):
        return (torch.zeros(1), torch.zeros(4), torch.tensor([self.score[i]], dtype=torch.float32), int(self.label[i]), 0)


class EvaluateTest(unittest.TestCase):
    def test_evaluate_aligns_the_tags_with_the_samples(self):
        rng = np.random.default_rng(5)
        n = 300
        label, score, pon = rng.choice([0, 1, 2], n), rng.normal(size=n), rng.choice([-1, 0, 1], n)
        ds = Rows(score, label)
        sampler = EpochSampler(n, None, False)
        loader = torch.utils.data.DataLoader(ds, batch_size=7, sampler=sampler)

        class Run:
            device, ddp = torch.device("cpu"), False
        out = train.evaluate(ConstantModel(), loader, torch.nn.CrossEntropyLoss(), Run(), False, None, pon)
        probs = torch.softmax(torch.stack([torch.zeros(n), torch.tensor(score, dtype=torch.float32), torch.zeros(n)], 1),
                              1).numpy()
        expect = metrics.postpon_report(label, probs, pon)
        self.assertEqual(out["somatic_threshold_postpon"], expect["threshold"])
        self.assertEqual(out["somatic_at_recall_postpon"], expect["at_recall"])
        self.assertEqual(out["somatic_at_recall_postpon_callable"], expect["at_recall_callable"])
        self.assertEqual(out["postpon"]["tagged"], int((pon == 1).sum()))
        self.assertIn("post-PoN F1", train.describe(out))

    def test_shuffled(self):
        """train.evaluate with a shuffled sampler (order != identity, as after a DDP gather): the tags must follow the
        dataset index, not the iteration order."""
        rng = np.random.default_rng(11)
        n = 500
        label, score = rng.choice([0, 1, 2], n, p=[0.6, 0.2, 0.2]), rng.normal(size=n)
        pon = np.where(label == 1, rng.choice([0, 1], n, p=[0.7, 0.3]), rng.choice([-1, 0, 1], n))
        sampler = EpochSampler(n, None, True, seed=3)
        order = list(sampler)
        self.assertNotEqual(order, list(range(n)))
        self.assertEqual(order, list(sampler))  # deterministic, as evaluate assumes

        class Run:
            device, ddp = torch.device("cpu"), False
        out = train.evaluate(ConstantModel(), torch.utils.data.DataLoader(Rows(score, label), batch_size=9,
                                                                          sampler=sampler),
                             torch.nn.CrossEntropyLoss(), Run(), False, None, pon)
        probs = torch.softmax(torch.stack([torch.zeros(n), torch.tensor(score, dtype=torch.float32), torch.zeros(n)], 1),
                              1).numpy()
        expect = metrics.postpon_report(label, probs, pon)
        self.assertEqual(out["somatic_threshold_postpon"], expect["threshold"])
        self.assertEqual(out["postpon"]["calls"], expect["calls"])
        self.assertEqual(out["somatic_at_recall_postpon"], expect["at_recall"])
        self.assertEqual(out["somatic_at_recall_postpon_callable"], expect["at_recall_callable"])
        # and a wrong alignment would differ (the test can see a misalignment)
        wrong = metrics.postpon_report(label, probs, pon[np.asarray(order)])
        self.assertNotEqual((wrong["threshold"], wrong["calls"]), (expect["threshold"], expect["calls"]))


class TrainPredictTest(unittest.TestCase):
    """train / predict --site-features on the small sets, --best-by postpon and the call rule."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.t = Path(self.tmp.name)
        self.root = self.t / "sample" / "tensors"
        make_tensor_set(self.root, SPEC, shard_size=8, seed=6)

    def tearDown(self):
        self.tmp.cleanup()

    def features(self, name, pon="random"):
        d = self.t / name
        data = {}
        for kind in ("SNV", "INDEL"):
            data[kind] = make_site_features(d, KindIndex(self.root / kind, self.t / "fcache"), seed=2, pon=pon)
        return d, data

    def validation_tags(self, data):
        """PoN tags of the validation tensors of the training run (block split of chr2 / chr3, SMALL_ARGS)."""
        tags = []
        for kind in ("SNV", "INDEL"):
            index = KindIndex(self.root / kind, self.t / "fcache")
            _, va = block_split(index, index.select(["chr2", "chr3"], labelled=False), 0.5, 21, 0)
            va = va[index.arrays["eval_label"][va] >= 0]
            tags.append(data[kind]["pon_tag"][va])
        return np.concatenate(tags)

    def test_site_features_train_resume_predict(self):
        feat, data = self.features("features")
        out = self.t / "run"
        train.main(["--tensors", str(self.root), "--output", str(out), "--epochs", "1", "--site-features", str(feat)]
                   + SMALL_ARGS)
        ckpt = torch.load(out / "best.pth", weights_only=False)
        self.assertEqual((ckpt["config"]["site_features"], ckpt["config"]["scalars"]), (F, 0))
        rec = ckpt["site_features"]
        self.assertEqual((rec["directory"], rec["version"], rec["names"]),
                         (str(feat.resolve()), site_features.VERSION, list(site_features.SITE_FEATURES)))
        self.assertEqual({k: v["candidates_sha1"] for k, v in rec["sets"].items()},
                         {f"sample.{k}": str(data[k]["candidates_sha1"]) for k in ("SNV", "INDEL")})
        self.assertEqual(ckpt["call_rule"], "recall")  # --best-by f1: graph_vcf keeps the recall-0.9 threshold
        shas = {f"sample.{k}": hashlib.sha256((feat / f"sample.{k}.npz").read_bytes()).hexdigest() for k in ("SNV", "INDEL")}
        self.assertEqual(ckpt["feature_files"], dict(site_features=shas))
        val = ckpt["val"]
        tags = self.validation_tags(data)
        self.assertEqual((val["postpon"]["tensors"], val["postpon"]["tagged"]), (len(tags), int((tags == 1).sum())))
        self.assertGreater(val["postpon"]["tagged"], 0)
        self.assertEqual(val["somatic_threshold_postpon"], val["postpon"]["threshold"])
        self.assertEqual(set(val["somatic_at_recall_postpon"]), {"0.5", "0.8", "0.9", "0.95"})
        self.assertEqual(set(val["somatic_at_recall_postpon_callable"]), {"0.5", "0.8", "0.9", "0.95"})
        self.assertEqual(val["truth_threshold_postpon"], val["postpon"]["truth"]["best"]["threshold"])
        row = json.loads((out / "metrics.jsonl").read_text().splitlines()[-1])
        self.assertEqual(row["val"]["somatic_threshold_postpon"], val["somatic_threshold_postpon"])
        self.assertIn("post-PoN F1", (out / "train.log").read_text())
        # resume keeps the input; a run without the features cannot resume it
        train.main(["--tensors", str(self.root), "--output", str(out), "--epochs", "2", "--site-features", str(feat),
                    "--resume", str(out / "last.pth")] + SMALL_ARGS)
        self.assertEqual(torch.load(out / "last.pth", weights_only=False)["epoch"], 2)
        with self.assertRaisesRegex(SystemExit, "differs from the arguments"):
            train.main(["--tensors", str(self.root), "--output", str(self.t / "r2"), "--epochs", "3", "--resume",
                        str(out / "last.pth")] + SMALL_ARGS)

        pred = self.t / "pred"
        predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(self.root), "--output", str(pred),
                      "--chroms", "chr1", "--num-workers", "0", "--amp", "off", "--kinds", "INDEL"])
        with gzip.open(pred / "sample.tensors.INDEL.predictions.ndjson.gz", "rt") as f:
            records = [json.loads(line) for line in f]
        index = KindIndex(self.root / "INDEL", self.t / "fcache")
        chr1 = index.select(["chr1"], labelled=False)
        self.assertEqual([r["pon_tag"] for r in records], data["INDEL"]["pon_tag"][chr1].tolist())
        report = json.loads((pred / "sample.tensors.INDEL.metrics.json").read_text())
        self.assertEqual(report["site_features"], dict(directory=str(feat.resolve()), version=site_features.VERSION))
        files = report["feature_files"]
        self.assertEqual((files["sha256"], files["changed"]), (shas["sample.INDEL"], False))
        self.assertEqual(report["feature_file_warnings"], [])
        rules = report["rules"]
        best = torch.load(out / "best.pth", weights_only=False)
        self.assertEqual(rules["val_postpon"]["threshold"], best["val"]["somatic_threshold_postpon"])
        self.assertEqual(rules["val_postpon_truth"]["threshold"], best["val"]["truth_threshold_postpon"])
        self.assertEqual(rules["threshold_t_pon"]["threshold"], best["somatic_threshold"])
        self.assertIn("pipeline_pon", rules)
        self.assertNotIn("val_recall_0.9_postpon_callable", rules)
        for rule in ("val_postpon", "pipeline_pon", "threshold_t_pon"):
            t = rules[rule]["threshold"]
            self.assertEqual(rules[rule]["per_tensor"]["calls"],
                             sum(r["in_test"] and r["pon_tag"] != 1 and r["p_somatic"] >= t for r in records))
        # the predictions depend on the features: other values, other probabilities
        changed, _ = self.features("changed")
        for kind in ("SNV", "INDEL"):  # same tags, other AFs
            path = changed / f"sample.{kind}.npz"
            with np.load(path) as z:
                arrays = {k: z[k] for k in z.files}
            arrays["popaf_exact"] = np.where(arrays["popaf_exact"] >= 0, 0.9, -1).astype(np.float32)
            np.savez(path, **arrays)
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(self.root), "--output",
                          str(self.t / "p2"), "--chroms", "chr1", "--num-workers", "0", "--amp", "off", "--kinds",
                          "INDEL", "--site-features", str(changed)])
        with gzip.open(self.t / "p2" / "sample.tensors.INDEL.predictions.ndjson.gz", "rt") as f:
            other = [json.loads(line) for line in f]
        self.assertNotEqual([r["p_somatic"] for r in records], [r["p_somatic"] for r in other])
        # a training set's file changed since training: a warning on stderr and in .metrics.json
        self.assertIn("warning: site-features file", err.getvalue())
        self.assertIn("sample.INDEL has changed since the checkpoint was trained", err.getvalue())
        report = json.loads((self.t / "p2" / "sample.tensors.INDEL.metrics.json").read_text())
        changed_file = report["feature_files"]
        self.assertEqual((changed_file["changed"], changed_file["checkpoint_sha256"]), (True, shas["sample.INDEL"]))
        self.assertEqual(len(report["feature_file_warnings"]), 1)
        # a set without its file stops predict before anything is written
        missing = self.t / "missing"
        missing.mkdir()
        with self.assertRaisesRegex(SystemExit, "no site features"):
            predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(self.root), "--output",
                          str(self.t / "p3"), "--num-workers", "0", "--amp", "off", "--site-features", str(missing)])
        self.assertEqual(list((self.t / "p3").glob("*.predictions.ndjson.gz")), [])
        # graph_vcf reads the post-PoN threshold and the call rule of the checkpoint
        cp = graph_vcf.read_checkpoint(out / "best.pth")
        self.assertEqual(graph_vcf.postpon_threshold(cp, "x"), best["val"]["somatic_threshold_postpon"])
        self.assertEqual(cp["call_rule"], "recall")

    def test_with_scalars(self):
        """--scalars and --site-features together: the scalar statistics are the SCALARS' only, the features follow."""
        feat, _ = self.features("features")
        out = self.t / "run"
        train.main(["--tensors", str(self.root), "--output", str(out), "--epochs", "1", "--kinds", "INDEL", "--scalars",
                    "--site-features", str(feat)] + SMALL_ARGS)
        ckpt = torch.load(out / "best.pth", weights_only=False)
        self.assertEqual((ckpt["config"]["scalars"], ckpt["config"]["site_features"]), (len(SCALARS), F))
        self.assertEqual(tuple(ckpt["model_state_dict"]["scalar_mean"].shape), (len(SCALARS),))
        self.assertEqual(tuple(ckpt["model_state_dict"]["scalar_mlp.0.weight"].shape), (64, len(SCALARS) + F))
        predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(self.root), "--output", str(self.t / "p"),
                      "--chroms", "chr1", "--num-workers", "0", "--amp", "off", "--kinds", "INDEL"])
        self.assertTrue((self.t / "p" / "sample.tensors.INDEL.metrics.json").exists())

    def test_missing_or_stale_features_stop_training(self):
        feat, data = self.features("features")
        (feat / "sample.SNV.npz").unlink()
        with self.assertRaisesRegex(SystemExit, "no site features"):
            train.main(["--tensors", str(self.root), "--output", str(self.t / "r"), "--epochs", "1",
                        "--site-features", str(feat)] + SMALL_ARGS)
        train.main(["--tensors", str(self.root), "--output", str(self.t / "r_indel"), "--epochs", "1", "--kinds",
                    "INDEL", "--site-features", str(feat)] + SMALL_ARGS)  # the INDEL file alone is enough for INDEL
        np.savez(feat / "sample.INDEL.npz", **{**data["INDEL"], "candidates_sha1": np.array("0" * 40)})
        with self.assertRaisesRegex(SystemExit, "stale"):
            train.main(["--tensors", str(self.root), "--output", str(self.t / "r2"), "--epochs", "1", "--kinds", "INDEL",
                        "--site-features", str(feat)] + SMALL_ARGS)

    def test_without_features_and_without_tags(self):
        feat, _ = self.features("features")
        # a run without features records none; predict refuses --site-features and graph_vcf --threshold postpon
        train.main(["--tensors", str(self.root), "--output", str(self.t / "plain"), "--epochs", "1"] + SMALL_ARGS)
        plain = torch.load(self.t / "plain" / "best.pth", weights_only=False)
        self.assertNotIn("somatic_threshold_postpon", plain["val"])
        self.assertIsNone(plain["feature_files"])
        self.assertEqual(plain["call_rule"], "recall")
        with self.assertRaisesRegex(SystemExit, "no post-PoN validation threshold"):
            graph_vcf.postpon_threshold(graph_vcf.read_checkpoint(self.t / "plain" / "best.pth"), "x")
        with self.assertRaisesRegex(SystemExit, "trained without site features"):
            predict.main(["--checkpoint", str(self.t / "plain" / "best.pth"), "--tensors", str(self.root), "--output",
                          str(self.t / "p1"), "--num-workers", "0", "--amp", "off", "--site-features", str(feat)])
        predict.main(["--checkpoint", str(self.t / "plain" / "best.pth"), "--tensors", str(self.root), "--output",
                      str(self.t / "p"), "--chroms", "chr1", "--num-workers", "0", "--amp", "off", "--kinds", "INDEL"])
        with gzip.open(self.t / "p" / "sample.tensors.INDEL.predictions.ndjson.gz", "rt") as f:
            self.assertFalse(any("pon_tag" in json.loads(line) for line in f))
        # without tags the post-PoN threshold is the plain one
        none, _ = self.features("none", pon="none")
        train.main(["--tensors", str(self.root), "--output", str(self.t / "run0"), "--epochs", "1", "--site-features",
                    str(none)] + SMALL_ARGS)
        val = torch.load(self.t / "run0" / "best.pth", weights_only=False)["val"]
        self.assertEqual(val["somatic_threshold_postpon"], val["threshold"])
        self.assertAlmostEqual(val["postpon"]["f1"], val["thresholded"]["somatic"]["f1"], places=12)

    def test_best_by_postpon_and_the_call_rule(self):
        feat, _ = self.features("features")
        with self.assertRaisesRegex(SystemExit, "needs the PoN tags"):
            train.main(["--tensors", str(self.root), "--output", str(self.t / "r0"), "--epochs", "1", "--best-by",
                        "postpon"] + SMALL_ARGS)
        out = self.t / "run"
        train.main(["--tensors", str(self.root), "--output", str(out), "--epochs", "3", "--site-features", str(feat),
                    "--best-by", "postpon"] + SMALL_ARGS)
        rows = [json.loads(line) for line in (out / "metrics.jsonl").read_text().splitlines()]
        best, epoch = (-1.0, float("inf")), None
        for r in rows:  # the highest post-PoN F1, ties to the lower loss
            score, loss = r["val"]["postpon"]["f1"], r["val"]["loss"]
            if score > best[0] or (score == best[0] and loss < best[1]):
                best, epoch = (score, loss), r["epoch"]
        ckpt = torch.load(out / "best.pth", weights_only=False)
        self.assertEqual((ckpt["epoch"], ckpt["best"]["select"], ckpt["best"]["score"]), (epoch, "postpon", best[0]))
        self.assertEqual(ckpt["somatic_threshold"], ckpt["val"]["threshold"])  # predict's default is unchanged
        self.assertEqual((ckpt["call_rule"], ckpt["call_recall"]), ("postpon-recall", 0.8))
        cp = graph_vcf.read_checkpoint(out / "best.pth")
        self.assertEqual((cp["call_rule"], cp["call_recall"], cp["postpon"]),
                         ("postpon-recall", 0.8, ckpt["val"]["somatic_threshold_postpon"]))
        self.assertEqual(cp["at_recall_postpon_callable"], ckpt["val"]["somatic_at_recall_postpon_callable"])
        at = ckpt["val"]["somatic_at_recall_postpon_callable"].get("0.8")
        if at:  # reachable on this fixture: the call rule's threshold
            self.assertEqual(graph_vcf.postpon_recall_threshold(cp, "best.pth", 0.8), at["threshold"])
        else:
            with self.assertRaisesRegex(SystemExit, "no post-PoN validation threshold at recall 0.8"):
                graph_vcf.postpon_recall_threshold(cp, "best.pth", 0.8)
        # --call-recall sets the recall of the call rule
        train.main(["--tensors", str(self.root), "--output", str(self.t / "r5"), "--epochs", "1", "--site-features",
                    str(feat), "--best-by", "postpon", "--call-recall", "0.5"] + SMALL_ARGS)
        c5 = torch.load(self.t / "r5" / "best.pth", weights_only=False)
        self.assertEqual((c5["call_rule"], c5["call_recall"]), ("postpon-recall", 0.5))
        # resumed with the default --best-by f1: post-PoN scores are not compared with F1 ones
        train.main(["--tensors", str(self.root), "--output", str(out), "--epochs", "4", "--site-features", str(feat),
                    "--resume", str(out / "last.pth")] + SMALL_ARGS)
        self.assertIn("best epoch chosen afresh", (out / "train.log").read_text())
        ckpt = torch.load(out / "best.pth", weights_only=False)
        self.assertEqual((ckpt["epoch"], ckpt["best"]["select"], ckpt["call_rule"], ckpt["call_recall"]),
                         (4, "f1", "recall", None))
        # checkpoints written before the call rule was recorded: --best-by postpon counts as postpon-recall at 0.8
        self.assertEqual(graph_vcf.call_rule_of(dict(args=dict(best_by="postpon"))), ("postpon-recall", 0.8))
        self.assertEqual(graph_vcf.call_rule_of(dict(args=dict(best_by="f1"))), ("recall", None))
        self.assertEqual(graph_vcf.call_rule_of(dict()), ("recall", None))  # every SNV model
        self.assertEqual(graph_vcf.call_rule_of(dict(call_rule="recall", args=dict(best_by="postpon"))),
                         ("recall", None))
        self.assertEqual(graph_vcf.call_rule_of(dict(call_rule="postpon-recall")), ("postpon-recall", 0.8))
        self.assertEqual(graph_vcf.call_rule_of(dict(call_rule="postpon-recall", call_recall=0.5)),
                         ("postpon-recall", 0.5))
        self.assertEqual(graph_vcf.call_rule_of(dict(call_rule="postpon")), ("postpon", None))


class NonReasonsTest(unittest.TestCase):
    """--non-reasons / --non-reasons-in-bed: label reasons as non inside the benchmark only, in train and predict."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.t = Path(self.tmp.name)
        self.root = self.t / "sample" / "tensors"
        make_tensor_set(self.root, SPEC, shard_size=8, seed=6)

    def tearDown(self):
        self.tmp.cleanup()

    def index(self, kind="SNV"):
        return KindIndex(self.root / kind, self.t / "cache")

    def test_in_somatic_bed(self):
        index = self.index()
        a = index.arrays
        inside = train.in_somatic_bed(index)
        expect = (a["pos0"] >= 0) & ~((a["pos0"] >= GAP[0]) & (a["pos0"] < GAP[1]))
        np.testing.assert_array_equal(inside, expect)  # the fixture's somatic BED leaves GAP out
        self.assertTrue(((a["pos0"] >= GAP[0]) & (a["pos0"] < GAP[1])).any())

    def test_as_non_inside_the_benchmark_only(self):
        index = self.index()
        a = index.arrays
        code = index.meta["reasons"].index("outside_confident_region")
        mask = a["reason"] == code
        in_gap = (a["pos0"] >= GAP[0]) & (a["pos0"] < GAP[1])
        self.assertTrue((mask & in_gap).any() and (mask & ~in_gap).any())  # both sides of the benchmark
        n, inside = train.as_non(index, ["outside_confident_region"])
        self.assertEqual((n, inside), (int(mask.sum()), int((mask & ~in_gap).sum())))
        np.testing.assert_array_equal(a["label"][mask & ~in_gap], 0)  # inside the somatic BED: non
        np.testing.assert_array_equal(a["eval_label"][mask & ~in_gap], 0)
        np.testing.assert_array_equal(a["label"][mask & in_gap], -1)  # outside the benchmark: still unused
        np.testing.assert_array_equal(a["eval_label"][mask & in_gap], -1)
        before = {k: index.arrays[k].copy() for k in ("label", "eval_label")}
        other = self.index()  # the in-BED variant of the same reason gives the same tensors here
        self.assertEqual(train.as_non_in_bed(other, ["outside_confident_region"]), inside)
        for k in ("label", "eval_label"):
            np.testing.assert_array_equal(other.arrays[k], before[k])
        # a reason inside the confident region: every tensor of it (the region is inside the somatic BED here)
        index = self.index()
        code = index.meta["reasons"].index("below_snv_min_af")
        mask = index.arrays["reason"] == code
        n, inside = train.as_non(index, ["below_snv_min_af"])
        self.assertEqual(n, inside)
        np.testing.assert_array_equal(index.arrays["label"][mask], 0)
        # the label files are not touched
        self.assertEqual(int((self.index().arrays["label"][mask] == -1).sum()), int(mask.sum()))

    def test_train_and_predict_agree(self):
        out = self.t / "run"
        train.main(["--tensors", str(self.root), "--output", str(out), "--epochs", "1", "--kinds", "SNV",
                    "--non-reasons", "below_snv_min_af", "--non-reasons-in-bed", "outside_confident_region"] + SMALL_ARGS)
        log = (out / "train.log").read_text()
        self.assertIn("inside the benchmark trained and validated as non", log)
        self.assertIn("inside the somatic BED trained and validated as non", log)
        ckpt = torch.load(out / "best.pth", weights_only=False)
        self.assertEqual((ckpt["args"]["non_reasons"], ckpt["args"]["non_reasons_in_bed"]),
                         (["below_snv_min_af"], ["outside_confident_region"]))
        pred = self.t / "pred"
        predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(self.root), "--output", str(pred),
                      "--chroms", "chr1", "chr2", "chr3", "--num-workers", "0", "--amp", "off", "--kinds", "SNV"])
        with gzip.open(pred / "sample.tensors.SNV.predictions.ndjson.gz", "rt") as f:
            records = [json.loads(line) for line in f]
        index = self.index()
        positions = index.select(None, labelled=False)
        inside = train.in_somatic_bed(index)[positions]
        seen = set()
        for r, ins in zip(records, inside):
            if r["reason"] in ("below_snv_min_af", "outside_confident_region"):
                seen.add((r["reason"], bool(ins)))
                self.assertEqual(r["test_label"], 0 if ins else None, r)
        self.assertIn(("outside_confident_region", True), seen)
        self.assertIn(("outside_confident_region", False), seen)
        report = json.loads((pred / "sample.tensors.SNV.metrics.json").read_text())
        self.assertEqual(report["as_non"], sum(1 for r, ins in zip(records, inside)
                                               if r["reason"] in ("below_snv_min_af", "outside_confident_region") and ins))


class UnlabelledTest(unittest.TestCase):
    """A set without labels (a new sample): indexed from its summaries, predicted without test metrics."""

    def test_predict_an_unlabelled_set(self):
        import shutil
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            labelled_root = t / "labelled" / "tensors"
            make_tensor_set(labelled_root, SPEC, shard_size=8, seed=6)
            root = t / "new" / "tensors"
            shutil.copytree(labelled_root, root)
            for kind in ("SNV", "INDEL"):  # the summaries carry grch38 (merge --reference-path), the labels go
                for summary in sorted((root / kind).glob("*_variant_summary.ndjson")):
                    labels = (root / kind / summary.name.replace("_variant_summary", "_labels")).read_text().splitlines()
                    lines = [json.dumps(dict(json.loads(s), grch38=json.loads(lab)["grch38"]))
                             for s, lab in zip(summary.read_text().splitlines(), labels)]
                    summary.write_text("\n".join(lines) + "\n")
                for f in [root / kind / "labels.manifest.json"] + list((root / kind).glob("*_labels.*")):
                    f.unlink()
            (root / "somatic.recall.tsv").unlink()
            index = KindIndex(root / "SNV", t / "cache")
            known = KindIndex(labelled_root / "SNV", t / "cache")
            self.assertEqual(index.meta["labels"]["format"], "unlabelled")
            self.assertEqual(set(index.meta["reasons"]), {"unlabelled"})
            self.assertTrue((index.arrays["label"] == -1).all() and index.arrays["in_region"].all())
            np.testing.assert_array_equal(index.arrays["off_reference"], known.arrays["off_reference"])
            on = ~known.arrays["off_reference"]
            np.testing.assert_array_equal(index.arrays["pos0"][on], known.arrays["pos0"][on])
            self.assertEqual(index.candidates(), known.candidates())
            out = t / "run"
            train.main(["--tensors", str(labelled_root), "--output", str(out), "--epochs", "1", "--kinds", "SNV",
                        "--non-reasons", "below_snv_min_af"] + SMALL_ARGS)
            preds = {}
            for name, r in (("labelled", labelled_root), ("new", root)):
                predict.main(["--checkpoint", str(out / "best.pth"), "--tensors", str(r), "--output", str(t / name),
                              "--chroms", "chr1", "chr2", "--num-workers", "0", "--amp", "off", "--kinds", "SNV",
                              "--norm", "train"])
                preds[name] = read_ndjson(t / name / f"{name}.tensors.SNV.predictions.ndjson.gz")
            self.assertEqual([r["candidate_id"] for r in preds["new"]], [r["candidate_id"] for r in preds["labelled"]])
            self.assertEqual([r["p_somatic"] for r in preds["new"]], [r["p_somatic"] for r in preds["labelled"]])
            self.assertTrue(all(r["test_label"] is None and not r["in_test"] and r["reason"] == "unlabelled"
                                for r in preds["new"]))
            report = json.loads((t / "new" / "new.tensors.SNV.metrics.json").read_text())
            self.assertEqual(report["labels"]["format"], "unlabelled")
        with tempfile.TemporaryDirectory() as tmp:  # a summary without grch38 cannot place the tensors
            root = Path(tmp) / "x" / "tensors"
            make_tensor_set(root, SPEC, shard_size=8, seed=6)
            (root / "SNV" / "labels.manifest.json").unlink()
            with self.assertRaisesRegex(ValueError, "needs the summaries' grch38"):
                KindIndex(root / "SNV")


def read_ndjson(path):
    with gzip.open(path, "rt") as f:
        return [json.loads(line) for line in f]


class BuilderTest(unittest.TestCase):
    """site_features.build_chrom on a small FASTA, a summary and two databases (a stand-in for bcftools query)."""

    #         0         1         2         3         4         5         6
    #         0123456789012345678901234567890123456789012345678901234567890123
    SEQ = "GATCCGTAAAAAAAAAGCTGACCTGACTGCATGCCATTGCAGGTCAGTCAGTACCGTTGACGTA" + "C" * 36

    def setUp(self):
        import pysam
        self.tmp = tempfile.TemporaryDirectory()
        t = self.t = Path(self.tmp.name)
        (t / "ref.fa").write_text(">chr1\n" + self.SEQ + "\n")
        pysam.faidx(str(t / "ref.fa"))
        g = lambda pos0, ref, alt: dict(chrom="chr1", pos0=pos0, ref=ref, alt=alt)  # noqa: E731
        rows = [
            ("a", g(14, "A", "")),       # one A of the 9-A run (7..15) deleted at its end: left-aligned to POS 7
            ("b", g(31, "", "A")),       # an A inserted after the A at 30: left-aligned to pos0 30 (POS 30, C>CA)
            ("c", None),                 # off-reference
            ("d", g(44, "C", "")),       # a deletion without a database record with an AF of its own
        ]
        (t / "INDEL").mkdir()
        (t / "INDEL" / "chr1_variant_summary.ndjson").write_text("".join(
            json.dumps(dict(candidate_id=c, node_id=1, grch38=gr)) + "\n" for c, gr in rows))
        # database records as bcftools query prints them (POS, REF, ALT, AF): the left-aligned deletion with AF 0.02
        # (raw match: tagged, hp 9 >= 8 needs 0.01), the insertion right-shifted (only the normalised AF matches),
        # d's allele without an AF and another INDEL at its POS
        (t / "gnomad.vcf.gz.txt").write_text("7\tTA\tT\t0.02\n31\tA\tAA\t0.5\n")
        (t / "colors.vcf.gz.txt").write_text("7\tTA\tT,TAA\t0.001,0.3\n44\tTC\tT\t.\n44\tTCA\tT\t0.004\n")
        fake = t / "bcftools"
        fake.write_text('#!/bin/sh\nfor last; do :; done\ncat "$last.txt"\n')
        fake.chmod(0o755)

    def tearDown(self):
        self.tmp.cleanup()

    def test_build_chrom(self):
        import pysam
        t = self.t
        out = site_features.build_chrom(str(t / "INDEL"), "chr1", "chr1_variant_summary.ndjson", str(t / "ref.fa"),
                                        {"gnomad": str(t / "gnomad.vcf.gz"), "colors": str(t / "colors.vcf.gz")},
                                        str(t / "bcftools"))
        self.assertEqual(out["candidates"], ["a", "b", "c", "d"])
        self.assertEqual(out["vcf_pos"].tolist(), [7, 30, -1, 44])
        self.assertEqual([x.decode() for x in out["vcf_ref"]], ["TA", "C", "", "TC"])
        self.assertEqual([x.decode() for x in out["vcf_alt"]], ["T", "CA", "", "T"])
        pon = site_features.pon_module()
        with pysam.FastaFile(str(t / "ref.fa")) as fasta:
            hp = [pon.homopolymer_run(fasta, "chr1", 7, "TA", "T"), pon.homopolymer_run(fasta, "chr1", 30, "C", "CA"),
                  0, pon.homopolymer_run(fasta, "chr1", 44, "TC", "T")]
        self.assertEqual(out["hp_len"].tolist(), hp)
        self.assertEqual(hp[0], 9)
        np.testing.assert_allclose(out["rawaf_gnomad"], [0.02, 0, -1, 0], rtol=1e-6)
        np.testing.assert_allclose(out["popaf_gnomad"], [0.02, 0.5, -1, 0], rtol=1e-6)  # the right-shifted INS
        np.testing.assert_allclose(out["popaf_colors"], [0.001, 0, -1, 0], rtol=1e-6)  # d: no AF
        self.assertEqual(out["pon_tag_gnomad"].tolist(), [1, 0, -1, 0])  # 0.02 >= 0.01 in a 9-base homopolymer
        self.assertEqual(out["pon_tag_colors"].tolist(), [0, 0, -1, 0])  # 0.001 < 0.01
        # the largest INDEL AF at the left-aligned POS: a's position also holds the colors INS (0.3); d's two records
        np.testing.assert_allclose(out["popaf_pos"], [0.3, 0.5, -1, 0.004], rtol=1e-6)


class GraphVcfArgsTest(unittest.TestCase):
    def test_threshold_names(self):
        base = ["--predictions", "p", "--tensors", "t", "--kind", "INDEL", "--graph-index", "g", "--output", "o"]
        for name in ("postpon", "postpon-recall", "recall", "predict"):
            self.assertEqual(graph_vcf.parse_args(base + ["--threshold", name]).threshold, name)
        self.assertEqual(graph_vcf.parse_args(base + ["--threshold", "0.25"]).threshold, 0.25)
        self.assertIsNone(graph_vcf.parse_args(base).threshold)  # the default: the checkpoint's call rule
        self.assertIsNone(graph_vcf.parse_args(base).target_recall)  # resolved per rule (recall 0.9, call_recall)
        for name in ("postpon9", "postpon-truth", "postpon-recall-all"):
            with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
                graph_vcf.parse_args(base + ["--threshold", name])


def pipeline_available():
    from ..vcfio import import_pipeline
    import_pipeline()
    try:
        import indexed_gam_pipeline_v4.graph_index  # noqa: F401
        import indexed_gam_pipeline_v4.tools.graph_prep  # noqa: F401
        import pysam  # noqa: F401
        return True
    except ImportError:
        return False


@unittest.skipUnless(pipeline_available(), "indexed_gam_pipeline_v4 not importable (run from the repo, or PYTHONPATH)")
class GraphVcfPostPonTest(unittest.TestCase):
    def test_postpon_threshold_calls(self):
        from .test_vcf import Fixture, columns
        with tempfile.TemporaryDirectory() as tmp:
            f = Fixture(tmp)
            ckpt = torch.load(f.checkpoint, weights_only=False)
            ckpt["val"].update(somatic_threshold_postpon=0.85, truth_threshold_postpon=0.25)
            torch.save(ckpt, f.checkpoint)
            stats, out = f.graph("SNV", threshold="postpon")
            self.assertEqual(stats["threshold"], 0.85)
            self.assertIn("post-PoN", stats["call"])
            self.assertEqual({r[2]: r[6] for r in columns(out)}["7:0:SNP:C>A"], "LowQual")
            self.assertEqual(f.graph("SNV")[0]["threshold"], 0.5)  # no call rule: validation recall 0.9, unchanged
            self.assertEqual(f.graph("SNV", threshold="postpon", min_af=0.07)[0]["threshold"], 0.85)
            ckpt["call_rule"] = "postpon"  # an explicit postpon call rule: the post-PoN best-F1 threshold
            torch.save(ckpt, f.checkpoint)
            stats, _ = f.graph("SNV")
            self.assertEqual(stats["threshold"], 0.85)
            self.assertEqual(f.graph("SNV", threshold="recall")[0]["threshold"], 0.5)
            self.assertEqual(f.graph("SNV", threshold=0.3)[0]["threshold"], 0.3)
            # the INDEL recipe's checkpoints: the post-PoN threshold at recall 0.8 of the callable positives
            ckpt["val"]["somatic_at_recall_postpon_callable"] = {"0.5": dict(threshold=0.95), "0.8": dict(threshold=0.65)}
            ckpt.update(call_rule="postpon-recall", call_recall=0.8)
            torch.save(ckpt, f.checkpoint)
            stats, _ = f.graph("SNV")
            self.assertEqual(stats["threshold"], 0.65)
            self.assertIn("post-PoN recall 0.8 of the callable positives", stats["call"])
            self.assertEqual(f.graph("SNV", target_recall=0.5)[0]["threshold"], 0.95)  # --target-recall overrides
            self.assertEqual(f.graph("SNV", threshold="postpon")[0]["threshold"], 0.85)
            self.assertEqual(f.graph("SNV", threshold="recall")[0]["threshold"], 0.5)  # recall rule: 0.9 by default
            del ckpt["call_rule"], ckpt["call_recall"]  # a code_final checkpoint: --best-by postpon, no call rule
            ckpt["args"] = dict(best_by="postpon")
            torch.save(ckpt, f.checkpoint)
            self.assertEqual(f.graph("SNV")[0]["threshold"], 0.65)
            del ckpt["val"]["somatic_at_recall_postpon_callable"]
            torch.save(ckpt, f.checkpoint)
            with self.assertRaisesRegex(SystemExit, "no post-PoN validation threshold at recall 0.8"):
                f.graph("SNV")
            del ckpt["val"]["somatic_threshold_postpon"]
            ckpt["call_rule"] = "postpon"
            torch.save(ckpt, f.checkpoint)
            with self.assertRaisesRegex(SystemExit, "no post-PoN validation threshold"):
                f.graph("SNV")


if __name__ == "__main__":
    unittest.main()

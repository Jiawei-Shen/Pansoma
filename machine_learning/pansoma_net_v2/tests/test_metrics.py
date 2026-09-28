import unittest

import numpy as np

from .. import metrics


def brute_ap(positive, score):
    """sklearn's definition, written out: precision at each distinct threshold, weighted by the recall gained."""
    out, prev = 0.0, 0.0
    for t in sorted(set(score.tolist()), reverse=True):
        called = score >= t
        tp = (called & positive).sum()
        rec = tp / positive.sum()
        out += (rec - prev) * tp / called.sum()
        prev = rec
    return out


class MetricsTest(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.labels = rng.choice([0, 1, 2], 400, p=[0.7, 0.1, 0.2])
        logits = rng.normal(size=(400, 3)) + 2.0 * np.eye(3)[self.labels]
        self.probs = np.exp(logits) / np.exp(logits).sum(1, keepdims=True)

    def test_average_precision(self):
        pos, score = np.array([1, 0, 1, 0, 0], bool), np.array([0.9, 0.8, 0.7, 0.3, 0.1])
        self.assertAlmostEqual(metrics.average_precision(pos, score), (1 + 2 / 3) / 2)
        tied = np.array([0.9, 0.9, 0.9, 0.2, 0.2])  # ties count together whatever their order
        self.assertAlmostEqual(metrics.average_precision(pos, tied), brute_ap(pos, tied))
        self.assertAlmostEqual(metrics.average_precision(pos[::-1], tied[::-1]), brute_ap(pos, tied))
        som = self.labels == 1
        self.assertAlmostEqual(metrics.average_precision(som, self.probs[:, 1]), brute_ap(som, self.probs[:, 1]))

    def test_best_threshold_is_the_best_f1_of_the_thresholded_call(self):
        t, p, r, f1 = metrics.best_somatic_threshold(self.labels, self.probs)
        report = metrics.report(self.labels, self.probs, t)["thresholded"]["somatic"]
        self.assertAlmostEqual(report["f1"], f1)
        self.assertAlmostEqual(report["precision"], p)
        self.assertAlmostEqual(report["recall"], r)
        for other in np.unique(self.probs[:, 1])[::7]:
            self.assertLessEqual(metrics.report(self.labels, self.probs, other)["thresholded"]["somatic"]["f1"],
                                 f1 + 1e-12)

    def test_call_and_ignored_labels(self):
        probs = np.array([[0.5, 0.2, 0.3], [0.1, 0.3, 0.6], [0.2, 0.7, 0.1]])
        self.assertEqual(metrics.call(probs).tolist(), [0, 2, 1])
        self.assertEqual(metrics.call(probs, 0.25).tolist(), [0, 1, 1])   # somatic at p >= 0.25
        self.assertEqual(metrics.call(probs, 0.8).tolist(), [0, 2, 0])    # else the larger of non / germline
        labels = np.array([0, -1, 1])
        report = metrics.report(labels, probs, 0.5)
        self.assertEqual(report["tensors"], 2)
        self.assertEqual(sum(map(sum, report["thresholded"]["confusion"])), 2)


class TruthReportTest(unittest.TestCase):
    """Against brute force: every truth allele once (duplicates collapse), misses without a tensor, false
    positives from labels 0 and 2, found truth alleles of the other kind as true calls (once, not in the
    recall), labelled-1 tensors of truth alleles that do not count left out."""

    def setUp(self):
        rng = np.random.default_rng(5)
        n = 300
        self.labels = rng.choice([-1, 0, 1, 2], n, p=[0.1, 0.6, 0.15, 0.15])
        self.probs = rng.dirichlet([1, 1, 1], n)
        self.truth = set(range(40))                       # truth 30..39 have no tensor
        self.other = set(range(200, 206))                 # the other kind's truth (e.g. INDEL for an SNV model)
        self.matches = []
        for k in range(n):
            if self.labels[k] == 1:  # 100: a truth that does not count
                self.matches.append({int(rng.integers(0, 30))} if k % 5 else {100} if k % 10 else {200 + k % 6})
            else:
                self.matches.append(set())

    def brute(self, t):
        found = set()
        for k in range(len(self.labels)):
            if self.labels[k] == 1 and self.probs[k, 1] >= t:
                found |= self.matches[k] & (self.truth | self.other)
        fp = sum(1 for k in range(len(self.labels)) if self.labels[k] in (0, 2) and self.probs[k, 1] >= t)
        tp, hits = len(found & self.truth), len(found)
        prec = hits / (hits + fp) if hits + fp else 0.0
        rec = tp / len(self.truth)
        return tp, fp, 2 * prec * rec / (prec + rec) if prec + rec else 0.0, prec, rec, hits - tp

    def test_counts_best_threshold_and_ap(self):
        r = metrics.truth_report(self.truth, self.matches, self.labels, self.probs, 0.4, self.other)
        tp, fp, f1, prec, rec, other_tp = self.brute(0.4)
        at = r["at_threshold"]
        self.assertEqual((at["tp"], at["other_tp"], at["fp"], at["fn"]), (tp, other_tp, fp, 40 - tp))
        self.assertGreater(other_tp, 0)
        self.assertAlmostEqual(at["f1"], f1)
        self.assertAlmostEqual(at["precision"], prec)
        reachable = set().union(*[self.matches[k] for k in range(len(self.labels)) if self.labels[k] == 1])
        self.assertEqual(r["with_tensor"], len(reachable & self.truth))
        self.assertEqual(r["other_kind_truth_with_tensor"], len(reachable & self.other))
        self.assertAlmostEqual(r["ceiling"], len(reachable & self.truth) / 40)
        somatic = [self.matches[k] for k in range(len(self.labels)) if self.labels[k] == 1]
        self.assertEqual((r["somatic_tensors"], r["other_kind_tensors"], r["unmatched_tensors"]),
                         (sum(bool(m & self.truth) for m in somatic), sum(bool(m & self.other) for m in somatic),
                          sum(m == {100} for m in somatic)))
        thresholds = sorted(set(self.probs[:, 1].tolist()))
        self.assertAlmostEqual(r["best"]["f1"], max(self.brute(t)[2] for t in thresholds))
        self.assertAlmostEqual(r["best"]["f1"], self.brute(r["best"]["threshold"])[2])
        ap, prev = 0.0, 0.0                               # sklearn's step sum over the distinct scores
        for t in sorted(thresholds, reverse=True):
            _, _, _, p_t, r_t, _ = self.brute(t)
            ap += (r_t - prev) * p_t
            prev = r_t
        self.assertAlmostEqual(r["ap"], ap)
        no_other = metrics.truth_report(self.truth, self.matches, self.labels, self.probs, 0.4)
        self.assertEqual(no_other["at_threshold"]["tp"], tp)                  # the recall does not change
        self.assertEqual(no_other["at_threshold"]["other_tp"], 0)

    def test_duplicates_count_once(self):
        labels = np.array([1, 1, 1, 0])
        probs = np.array([[0, 0.9, 0.1], [0, 0.8, 0.2], [0, 0.7, 0.3], [0, 0.6, 0.4]])
        r = metrics.truth_report({"a", "b"}, [{"a"}, {"a"}, {"a"}, set()], labels, probs, 0.5)
        self.assertEqual((r["at_threshold"]["tp"], r["at_threshold"]["fp"], r["at_threshold"]["fn"]), (1, 1, 1))
        self.assertAlmostEqual(r["at_threshold"]["precision"], 0.5)   # 3 calls of one truth are one true positive
        self.assertEqual((r["somatic_tensors"], r["other_kind_tensors"]), (3, 0))
        # two tensors of one other-kind truth "c": one true call, not two, and not in the recall
        r = metrics.truth_report({"b"}, [{"c"}, {"c"}, {"b"}, set()], labels, probs, 0.5, other={"c"})
        at = r["at_threshold"]
        self.assertEqual((at["tp"], at["other_tp"], at["fp"], at["fn"]), (1, 1, 1, 0))
        self.assertAlmostEqual(at["precision"], 2 / 3)
        self.assertEqual((r["somatic_tensors"], r["other_kind_tensors"], r["unmatched_tensors"]), (1, 2, 0))
        r = metrics.truth_report({"b"}, [{"c"}, {"c"}, {"b"}, set()], labels, probs, 0.5)   # "c" counts for neither
        self.assertEqual((r["at_threshold"]["tp"], r["at_threshold"]["other_tp"], r["unmatched_tensors"]), (1, 0, 2))
        self.assertAlmostEqual(r["at_threshold"]["precision"], 0.5)


if __name__ == "__main__":
    unittest.main()

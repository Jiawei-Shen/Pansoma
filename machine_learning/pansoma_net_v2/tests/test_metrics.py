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
    positives from labels 0 and 2, labelled-1 tensors of truth alleles that do not count left out."""

    def setUp(self):
        rng = np.random.default_rng(5)
        n = 300
        self.labels = rng.choice([-1, 0, 1, 2], n, p=[0.1, 0.6, 0.15, 0.15])
        self.probs = rng.dirichlet([1, 1, 1], n)
        self.truth = set(range(40))                       # truth 30..39 have no tensor
        self.matches = []
        for k in range(n):
            if self.labels[k] == 1:
                self.matches.append({int(rng.integers(0, 30))} if k % 5 else {100})  # 100: a truth that does not count
            else:
                self.matches.append(set())

    def brute(self, t):
        found = set()
        for k in range(len(self.labels)):
            if self.labels[k] == 1 and self.probs[k, 1] >= t:
                found |= self.matches[k] & self.truth
        fp = sum(1 for k in range(len(self.labels)) if self.labels[k] in (0, 2) and self.probs[k, 1] >= t)
        tp = len(found)
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / len(self.truth)
        return tp, fp, 2 * prec * rec / (prec + rec) if prec + rec else 0.0, prec, rec

    def test_counts_best_threshold_and_ap(self):
        r = metrics.truth_report(self.truth, self.matches, self.labels, self.probs, 0.4)
        tp, fp, f1, prec, rec = self.brute(0.4)
        at = r["at_threshold"]
        self.assertEqual((at["tp"], at["fp"], at["fn"]), (tp, fp, 40 - tp))
        self.assertAlmostEqual(at["f1"], f1)
        reachable = set().union(*[self.matches[k] for k in range(len(self.labels)) if self.labels[k] == 1]) & self.truth
        self.assertEqual(r["with_tensor"], len(reachable))
        self.assertAlmostEqual(r["ceiling"], len(reachable) / 40)
        thresholds = sorted(set(self.probs[:, 1].tolist()))
        self.assertAlmostEqual(r["best"]["f1"], max(self.brute(t)[2] for t in thresholds))
        self.assertAlmostEqual(r["best"]["f1"], self.brute(r["best"]["threshold"])[2])
        ap, prev = 0.0, 0.0                               # sklearn's definition over the distinct scores
        for t in sorted(thresholds, reverse=True):
            _, _, _, p_t, r_t = self.brute(t)
            ap += (r_t - prev) * p_t
            prev = r_t
        self.assertAlmostEqual(r["ap"], ap)

    def test_duplicates_count_once(self):
        labels = np.array([1, 1, 1, 0])
        probs = np.array([[0, 0.9, 0.1], [0, 0.8, 0.2], [0, 0.7, 0.3], [0, 0.6, 0.4]])
        r = metrics.truth_report({"a", "b"}, [{"a"}, {"a"}, {"a"}, set()], labels, probs, 0.5)
        self.assertEqual((r["at_threshold"]["tp"], r["at_threshold"]["fp"], r["at_threshold"]["fn"]), (1, 1, 1))
        self.assertAlmostEqual(r["at_threshold"]["precision"], 0.5)   # 3 calls of one truth are one true positive
        self.assertEqual((r["somatic_tensors"], r["other_truth_tensors"]), (3, 0))
        r = metrics.truth_report({"b"}, [{"a"}, {"a"}, {"b"}, set()], labels, probs, 0.5)   # "a" does not count
        self.assertEqual((r["somatic_tensors"], r["other_truth_tensors"], r["negatives"]), (1, 2, 1))
        self.assertEqual((r["at_threshold"]["tp"], r["at_threshold"]["fp"]), (1, 1))      # neither a hit nor a false call


if __name__ == "__main__":
    unittest.main()

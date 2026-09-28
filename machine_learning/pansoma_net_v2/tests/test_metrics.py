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


if __name__ == "__main__":
    unittest.main()

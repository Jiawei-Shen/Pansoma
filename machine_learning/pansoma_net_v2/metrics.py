"""Metrics over (label, probabilities): per-class precision / recall / F1 of the argmax call, and for somatic and
germline the one-vs-rest average precision (PR-AUC). The somatic call can use a threshold: a tensor is somatic
when p_somatic >= t, otherwise the larger of p_non and p_germline. `best_somatic_threshold` picks the t of the
highest somatic F1 (on validation; the same t is then applied to the test set).

`truth_report` scores against the somatic truth VCF instead of per tensor: every truth allele counts once
(duplicates, e.g. several INDEL tensors matching one truth partially, collapse into its best tensor), truth
alleles without any tensor are misses, and a false positive is a somatic call on a tensor labelled 0 or 2."""
import numpy as np

CLASSES = ("non", "somatic", "germline")
SOMATIC = 1


def confusion(labels, pred):
    k = len(CLASSES)
    return np.bincount(labels.astype(np.int64) * k + pred, minlength=k * k).reshape(k, k)


def per_class(cm):
    cm = cm.astype(np.float64)
    out = {"accuracy": float(np.trace(cm) / max(cm.sum(), 1))}
    for k, name in enumerate(CLASSES):
        tp, fp, fn = cm[k, k], cm[:, k].sum() - cm[k, k], cm[k, :].sum() - cm[k, k]
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        out[name] = {"precision": float(prec), "recall": float(rec),
                     "f1": float(2 * prec * rec / (prec + rec)) if prec + rec else 0.0, "support": int(cm[k, :].sum())}
    return out


def average_precision(positive, score):
    """Area under the precision-recall curve as sklearn's average_precision_score: sum over the distinct scores of
    (recall step) x precision, tied scores taken together."""
    positive = np.asarray(positive, bool)
    n_pos = int(positive.sum())
    if n_pos == 0:
        return 0.0
    order = np.argsort(-score, kind="stable")
    s = score[order]
    tp = np.cumsum(positive[order])
    last = np.r_[s[1:] != s[:-1], True]
    tp, k = tp[last], np.flatnonzero(last) + 1
    precision, recall = tp / k, tp / n_pos
    return float(np.sum(np.diff(np.r_[0.0, recall]) * precision))


def best_somatic_threshold(labels, probs):
    """(threshold, precision, recall, f1) maximizing somatic F1 over the distinct p_somatic values."""
    positive = labels == SOMATIC
    n_pos = int(positive.sum())
    if n_pos == 0:
        return 0.5, 0.0, 0.0, 0.0
    score = probs[:, SOMATIC]
    order = np.argsort(-score, kind="stable")
    s, pos = score[order], positive[order]
    tp = np.cumsum(pos)
    k = np.arange(1, len(s) + 1)
    last = np.r_[s[1:] != s[:-1], True]  # a threshold can only fall between distinct scores
    prec, rec = tp / k, tp / n_pos
    f1 = np.where(prec + rec > 0, 2 * prec * rec / np.maximum(prec + rec, 1e-12), 0.0)
    f1 = np.where(last, f1, -1)
    i = int(np.argmax(f1))
    return float(s[i]), float(prec[i]), float(rec[i]), float(f1[i])


def call(probs, threshold=None):
    """Class calls: argmax, or somatic when p_somatic >= threshold and else the larger of non / germline."""
    if threshold is None:
        return probs.argmax(1)
    other = np.where(probs[:, 2] > probs[:, 0], 2, 0)
    return np.where(probs[:, SOMATIC] >= threshold, SOMATIC, other)


def report(labels, probs, threshold=None):
    """Metrics of labelled tensors (label >= 0): argmax, average precisions, and (if given) the thresholded call."""
    keep = labels >= 0
    labels, probs = labels[keep].astype(np.int64), probs[keep]
    out = {"tensors": int(len(labels)), "argmax": per_class(confusion(labels, call(probs)))}
    out["argmax"]["confusion"] = confusion(labels, call(probs)).tolist()
    out["somatic_ap"] = average_precision(labels == SOMATIC, probs[:, SOMATIC])
    out["germline_ap"] = average_precision(labels == 2, probs[:, 2])
    if threshold is not None:
        cm = confusion(labels, call(probs, threshold))
        out["threshold"] = float(threshold)
        out["thresholded"] = per_class(cm)
        out["thresholded"]["confusion"] = cm.tolist()
    return out


def truth_items(n_truth, truth_scores, negative_scores):
    """(positive, score) items: one per truth allele (its best tensor's p_somatic; -inf without a tensor) and
    one per negative tensor."""
    pos = np.full(n_truth, -np.inf)
    pos[:len(truth_scores)] = truth_scores
    score = np.concatenate([pos, np.asarray(negative_scores, np.float64)])
    return np.r_[np.ones(n_truth, bool), np.zeros(len(negative_scores), bool)], score


def truth_ap(n_truth, truth_scores, negative_scores):
    """Average precision with every truth allele in the recall denominator: truth alleles without a tensor are
    never retrieved, so they only lower the recall steps (average_precision would credit them at score -inf)."""
    if n_truth == 0:
        return 0.0
    score = np.concatenate([np.asarray(truth_scores, np.float64), np.asarray(negative_scores, np.float64)])
    positive = np.r_[np.ones(len(truth_scores), bool), np.zeros(len(negative_scores), bool)]
    order = np.argsort(-score, kind="stable")
    s, pos = score[order], positive[order]
    tp = np.cumsum(pos)
    last = np.r_[s[1:] != s[:-1], True]
    tp, k = tp[last], np.flatnonzero(last) + 1
    return float(np.sum(np.diff(np.r_[0.0, tp / n_truth]) * (tp / k)))


def truth_counts(n_truth, truth_scores, negative_scores, threshold):
    tp = int((np.asarray(truth_scores) >= threshold).sum())
    fp = int((np.asarray(negative_scores) >= threshold).sum())
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / n_truth if n_truth else 0.0
    return dict(tp=tp, fp=fp, fn=n_truth - tp, precision=prec, recall=rec,
                f1=2 * prec * rec / (prec + rec) if prec + rec else 0.0)


def best_truth_threshold(n_truth, truth_scores, negative_scores):
    """(threshold, counts) of the best truth-level F1 over the distinct finite scores."""
    positive, score = truth_items(n_truth, truth_scores, negative_scores)
    finite = np.isfinite(score)
    if not finite.any() or n_truth == 0:
        return 0.5, truth_counts(n_truth, truth_scores, negative_scores, 0.5)
    order = np.argsort(-score[finite], kind="stable")
    s, pos = score[finite][order], positive[finite][order]
    tp, fp = np.cumsum(pos), np.cumsum(~pos)
    last = np.r_[s[1:] != s[:-1], True]
    f1 = np.where(last, 2 * tp / np.maximum(2 * tp + fp + (n_truth - tp), 1), -1)
    t = float(s[int(np.argmax(f1))])
    return t, truth_counts(n_truth, truth_scores, negative_scores, t)


def truth_report(truth, matches, labels, probs, threshold):
    """Truth-level metrics. truth: the truth keys that count (any hashable); matches: per tensor the truth keys it
    stands for (tensors labelled 1); labels / probs: per tensor (evaluation labels; -1 left out). A truth allele
    scores the best p_somatic of its tensors; a labelled-1 tensor whose truth alleles do not count is left out
    (counted in other_truth_tensors: e.g. an SNV tensor that is a partial match of an INDEL truth)."""
    truth = set(truth)
    best = {}
    negatives = []
    somatic_tensors = other = 0
    for k in range(len(labels)):
        if labels[k] < 0:
            continue
        p = float(probs[k, SOMATIC])
        if labels[k] == SOMATIC:
            keys = set(matches[k]) & truth
            somatic_tensors += bool(keys)
            other += not keys
            for key in keys:
                best[key] = max(p, best.get(key, -np.inf))
        else:
            negatives.append(p)
    scores = np.array(list(best.values()), np.float64)
    out = dict(truth_alleles=len(truth), with_tensor=len(best), ceiling=len(best) / len(truth) if truth else 0.0,
               somatic_tensors=somatic_tensors, other_truth_tensors=other,
               ap=truth_ap(len(truth), scores, negatives), negatives=len(negatives),
               at_threshold=dict(threshold=float(threshold), **truth_counts(len(truth), scores, negatives, threshold)))
    t, counts = best_truth_threshold(len(truth), scores, negatives)
    out["best"] = dict(threshold=t, **counts)
    return out

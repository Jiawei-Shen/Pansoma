"""Metrics over (label, probabilities): per-class precision / recall / F1 of the argmax call, and for somatic and
germline the one-vs-rest average precision (PR-AUC). The somatic call can use a threshold: a tensor is somatic
when p_somatic >= t, otherwise the larger of p_non and p_germline. `best_somatic_threshold` picks the t of the
highest somatic F1 (on validation; the same t is then applied to the test set)."""
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

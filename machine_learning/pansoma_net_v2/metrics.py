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


def truth_curve(own, other, negatives):
    """Per distinct score from the highest: the score, and at a threshold of that score the truth alleles found
    (own kind, other kind) and the false calls. own / other: each truth allele's best score; negatives: one per
    tensor scored 0 or 2."""
    score = np.concatenate([np.asarray(own, np.float64), np.asarray(other, np.float64),
                            np.asarray(negatives, np.float64)])
    kind = np.r_[np.zeros(len(own), np.int8), np.ones(len(other), np.int8), np.full(len(negatives), 2, np.int8)]
    order = np.argsort(-score, kind="stable")
    s, k = score[order], kind[order]
    last = np.r_[s[1:] != s[:-1], True] if len(s) else np.zeros(0, bool)
    return s[last], np.cumsum(k == 0)[last], np.cumsum(k == 1)[last], np.cumsum(k == 2)[last]


def truth_prf(n_truth, tp, other_tp, fp):
    """Precision over the distinct truth alleles found (either kind) and the false calls; recall over the
    n_truth truth alleles of the model's kind."""
    tp, other_tp, fp = (np.asarray(v, np.float64) for v in (tp, other_tp, fp))
    hits = tp + other_tp
    prec = np.where(hits + fp > 0, hits / np.maximum(hits + fp, 1), 0.0)
    rec = tp / n_truth if n_truth else np.zeros_like(tp)
    f1 = np.where(prec + rec > 0, 2 * prec * rec / np.maximum(prec + rec, 1e-300), 0.0)
    return prec, rec, f1


def truth_counts(n_truth, own, other, negatives, threshold):
    tp = int((np.asarray(own) >= threshold).sum())
    other_tp = int((np.asarray(other) >= threshold).sum())
    fp = int((np.asarray(negatives) >= threshold).sum())
    prec, rec, f1 = (float(v) for v in truth_prf(n_truth, tp, other_tp, fp))
    return dict(tp=tp, other_tp=other_tp, fp=fp, fn=n_truth - tp, precision=prec, recall=rec, f1=f1)


def truth_ap(n_truth, own, other, negatives):
    """Average precision (sklearn's step sum) with every truth allele of the model's kind in the recall
    denominator: truth alleles without a tensor are never retrieved, so they only lower the recall steps."""
    if n_truth == 0:
        return 0.0
    _, tp, other_tp, fp = truth_curve(own, other, negatives)
    prec, rec, _ = truth_prf(n_truth, tp, other_tp, fp)
    return float(np.sum(np.diff(np.r_[0.0, rec]) * prec))


def best_truth_threshold(n_truth, own, other, negatives):
    """(threshold, counts) of the best truth-level F1 over the distinct scores."""
    s, tp, other_tp, fp = truth_curve(own, other, negatives)
    if len(s) == 0 or n_truth == 0:
        return 0.5, truth_counts(n_truth, own, other, negatives, 0.5)
    t = float(s[int(np.argmax(truth_prf(n_truth, tp, other_tp, fp)[2]))])
    return t, truth_counts(n_truth, own, other, negatives, t)


def truth_report(truth, matches, labels, probs, threshold, other=()):
    """Truth-level metrics: every truth allele counts once.
    truth: the truth keys of the model's kind (any hashable; the recall denominator); other: the truth keys of
    the other kind (an SNV tensor can partially match an INDEL truth); matches: per tensor the truth keys it
    stands for (tensors labelled 1); labels / probs: per tensor (evaluation labels; -1 left out).
    A truth allele scores the best p_somatic of its tensors, so duplicate tensors of one truth count once; truth
    alleles of the model's kind without a tensor are misses. Found other-kind truth alleles are true calls
    (precision) but not in the recall, which is over the model's kind. A false call is a call on a tensor scored
    0 or 2. A labelled-1 tensor whose truth alleles count for neither (not PASS, outside the BED) is left out."""
    truth = set(truth)
    other = set(other) - truth
    counted = truth | other
    best = {}
    negatives = []
    own_tensors = other_tensors = unmatched = 0
    for k in range(len(labels)):
        if labels[k] < 0:
            continue
        p = float(probs[k, SOMATIC])
        if labels[k] == SOMATIC:
            keys = set(matches[k]) & counted
            own_tensors += bool(keys & truth)
            other_tensors += bool(keys) and not keys & truth
            unmatched += not keys
            for key in keys:
                best[key] = max(p, best.get(key, -np.inf))
        else:
            negatives.append(p)
    own = np.array([v for key, v in best.items() if key in truth], np.float64)
    found_other = np.array([v for key, v in best.items() if key in other], np.float64)
    out = dict(truth_alleles=len(truth), with_tensor=len(own), ceiling=len(own) / len(truth) if truth else 0.0,
               other_kind_truth_with_tensor=len(found_other), somatic_tensors=own_tensors,
               other_kind_tensors=other_tensors, unmatched_tensors=unmatched,
               ap=truth_ap(len(truth), own, found_other, negatives), negatives=len(negatives),
               at_threshold=dict(threshold=float(threshold),
                                 **truth_counts(len(truth), own, found_other, negatives, threshold)))
    t, counts = best_truth_threshold(len(truth), own, found_other, negatives)
    out["best"] = dict(threshold=t, **counts)
    return out

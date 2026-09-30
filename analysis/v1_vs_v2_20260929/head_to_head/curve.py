"""Truth-level PR curve helpers: units = distinct truths (TP) and distinct false sites (FP), each with max score."""
import numpy as np


def curve(tp_scores, fp_scores, n_truth):
    """tp_scores: max score per truth that has any scored unit; fp_scores: max score per false site."""
    tp_scores = np.asarray(tp_scores, dtype=np.float64)
    fp_scores = np.asarray(fp_scores, dtype=np.float64)
    thr = np.unique(np.concatenate([tp_scores, fp_scores]))[::-1]
    tps = np.sort(tp_scores)
    fps = np.sort(fp_scores)
    tp = len(tps) - np.searchsorted(tps, thr, side="left")
    fp = len(fps) - np.searchsorted(fps, thr, side="left")
    rec = tp / n_truth
    prec = np.where(tp + fp > 0, tp / np.maximum(tp + fp, 1), 1.0)
    f1 = np.where(prec + rec > 0, 2 * prec * rec / np.maximum(prec + rec, 1e-12), 0)
    return dict(thr=thr, tp=tp, fp=fp, rec=rec, prec=prec, f1=f1)


def summary(c, n_truth, ceiling_n):
    out = dict(n_truth=n_truth, ceiling=ceiling_n / n_truth, ceiling_n=ceiling_n)
    for p in (0.05, 0.10, 0.15, 0.20):
        m = c["prec"] >= p
        if m.any():
            i = np.argmax(np.where(m, c["rec"], -1))
            out[f"R@P{p:.2f}"] = (round(float(c["rec"][i]), 4), int(c["tp"][i]), int(c["fp"][i]), float(c["thr"][i]))
        else:
            out[f"R@P{p:.2f}"] = None
    for r in (0.5, 0.7, 0.8, 0.9):
        m = c["rec"] >= r
        if m.any():
            i = np.argmax(np.where(m, c["prec"], -1))
            out[f"P@R{r:.1f}"] = (round(float(c["prec"][i]), 4), int(c["tp"][i]), int(c["fp"][i]), float(c["thr"][i]))
        else:
            out[f"P@R{r:.1f}"] = None
    i = int(np.argmax(c["f1"]))
    out["maxF1"] = (round(float(c["f1"][i]), 4), round(float(c["prec"][i]), 4), round(float(c["rec"][i]), 4),
                    int(c["tp"][i]), int(c["fp"][i]), float(c["thr"][i]))
    return out

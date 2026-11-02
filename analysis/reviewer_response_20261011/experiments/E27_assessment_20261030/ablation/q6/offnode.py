"""Off-node column mask (assessment 2026-10-30, R1.1 cross-node context ablation; not part of the released code).

With --mask-off-node DIR, every cell of a tensor that lies outside the candidate's own node visit is set to 0 in all 8
channels, i.e. made identical to a cell without evidence (indexed_gam_pipeline_v4 README: "a cell without evidence is
0 in every channel"). The site columns (insertion slots + span) are always kept. The window then holds only the
candidate node's bases, as a node-confined tensor would.

Own-node rule (per row; rows are in the candidate node's forward orientation; a cell consumes one graph base when its
graph base, channel 5, is 1-5):
  left of the anchor column a:  own <=> #consuming cells in [k, a) <= start
  slots [a, a + slots):          own
  right, k >= r0 = a + slots:    own <=> start + #(consuming or uncovered) cells in [r0, k) < node length
Checked against the per-column node coordinates of --debug-rows builds (q6/realcheck, q6/g8work).

DIR holds <SET>.<KIND>.side.npz (q6/side/build_side.py): per tensor of the set, in KindIndex order, node, start,
slots, span, node_len.
"""
from pathlib import Path

import numpy as np

DIR = None  # set by train / predict --mask-off-node before any TensorDataset is built


def side_arrays(index, positions):
    """(start, slots, span, node_len) of index's tensors at positions, or None without --mask-off-node."""
    if DIR is None:
        return None
    f = Path(DIR) / f"{index.dir.parent.parent.name}.{index.dir.name}.side.npz"
    with np.load(f) as z:
        if len(z["node"]) != len(index) or not (z["node"] == index.arrays["node"]).all():
            raise SystemExit(f"{f}: node order differs from the index of {index.dir}")
        return tuple(z[k][positions] for k in ("start", "slots", "span", "node_len"))


def own_mask(x, start, slots, span, node_len):
    """x int8 (B, 8, R, W) -> bool (B, R, W): True on the candidate's own node visit (and the site columns)."""
    B, _, R, W = x.shape
    a = W // 2
    start = np.asarray(start, np.int64)[:, None, None]
    node_len = np.asarray(node_len, np.int64)[:, None, None]
    slots = np.asarray(slots, np.int64)
    span = np.asarray(span, np.int64)
    cons = (x[:, 5] >= 1) & (x[:, 5] <= 5)
    cov = x[:, 4] > 0
    col = np.arange(W)[None, None, :]
    left = np.zeros((B, R, W), np.int64)
    left[:, :, :a] = cons[:, :, :a][:, :, ::-1].cumsum(2)[:, :, ::-1]  # #consuming in [k, a)
    cr = cons | ~cov
    cs = np.concatenate([np.zeros((B, R, 1), np.int64), cr.cumsum(2)], 2)  # cs[..., k] = # in [0, k)
    r0 = (a + slots)
    base = cs[np.arange(B)[:, None], np.arange(R)[None, :], r0[:, None]][:, :, None]  # # in [0, r0)
    before = cs[:, :, :W] - base  # # in [r0, k) for k >= r0
    r0 = r0[:, None, None]
    own = np.where(col < a, left <= start, np.where(col < r0, True, start + before < node_len))
    own |= (col >= a) & (col < r0 + span[:, None, None])  # the site columns
    return own


def apply(x, start, slots, span, node_len):
    """Zero, in place, every channel of the off-node cells of x (B, 8, R, W) or (8, R, W)."""
    single = x.ndim == 3
    xb = x[None] if single else x
    own = own_mask(xb, np.atleast_1d(start), np.atleast_1d(slots), np.atleast_1d(span), np.atleast_1d(node_len))
    xb *= own[:, None].astype(xb.dtype)
    return x

"""E28 (2): the nine final INDEL release models rerun with the off-path-excluded selection score logged. Per run: whether
the rerun reproduces the original per-epoch metrics (max |difference| of train loss and of the selection score, and the
first epoch where they differ), the epoch kept by the original selection (post-filter validation F1 with the off-path
negatives) and by the selection without them, and both scores. Writes select.md / select.json here."""
import json
from pathlib import Path
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
HERE = Path(__file__).resolve().parent
def load(f):
    return [json.loads(l) for l in open(f) if l.strip()] if f.exists() else []
def sel(v):  # selection score used by --best-by postpon
    return v["postpon"]["f1"] if isinstance(v.get("postpon"), dict) else v.get("postpon_f1")
rows = {}
for s in (0, 1, 2):
    for p in ("Illumina", "PacBio", "ONT"):
        run = f"release_s{s}_{p}_INDEL"
        o, n = load(R / run / "metrics.jsonl"), load(R / "selcheck_20261101" / run / "metrics.jsonl")
        r = dict(epochs_orig=len(o), epochs_rerun=len(n))
        if n:
            dl, ds, first = 0.0, 0.0, None
            for a, b in zip(o, n):
                x, y = abs(a["train_loss"] - b["train_loss"]), abs(sel(a["val"]) - sel(b["val"]))
                dl, ds = max(dl, x), max(ds, y)
                if first is None and (x > 0 or y > 0):
                    first = b["epoch"]
            r.update(max_dloss=dl, max_dsel=ds, first_diff_epoch=first)
            done = len(n) == 12
            sc = [(b["epoch"], sel(b["val"]), b["val"]["excl_offpath"]["postpon_f1"], b["val"]["excl_offpath"]["removed"]) for b in n]
            ea = max(sc, key=lambda t: (t[1], -t[0]))  # strict improvement keeps the earliest maximum
            eb = max(sc, key=lambda t: (t[2], -t[0]))
            r.update(complete=done, removed=sc[0][3], sel_epoch=ea[0], sel_score=ea[1], sel_score_excl=ea[2],
                     excl_epoch=eb[0], excl_score=eb[2], excl_score_with=eb[1],
                     per_epoch=[dict(epoch=e, with_=a, excl=b) for e, a, b, _ in sc])
        if o:
            oa = max(((a["epoch"], sel(a["val"])) for a in o), key=lambda t: (t[1], -t[0]))
            r.update(orig_sel_epoch=oa[0], orig_sel_score=oa[1])
        rows[run] = r
f = lambda x: "-" if x is None else (str(x) if isinstance(x, (int, bool)) else f"{x:.4f}")
L = ["| run | epochs (rerun / orig) | reproduces original (max Δloss; max Δselection; first differing epoch) | removed off-path negatives | original epoch (score) | rerun: selected epoch with off-path (score; score without) | rerun: selected epoch without off-path (score) | same epoch |",
     "|---|---|---|---|---|---|---|---|"]
for k, r in rows.items():
    if not r.get("epochs_rerun"):
        L.append(f"| {k} | 0 / {r['epochs_orig']} | running or pending | | | | | |"); continue
    rep = "yes, bit-identical" if r["max_dloss"] == 0 and r["max_dsel"] == 0 else f"no ({r['max_dloss']:.1e}; {r['max_dsel']:.4f}; epoch {r['first_diff_epoch']})"
    L.append(f"| {k} | {r['epochs_rerun']} / {r['epochs_orig']} | {rep} | {r['removed']} | {r.get('orig_sel_epoch')} ({f(r.get('orig_sel_score'))}) | "
             f"{r['sel_epoch']} ({f(r['sel_score'])}; {f(r['sel_score_excl'])}) | {r['excl_epoch']} ({f(r['excl_score'])}) | {r['sel_epoch'] == r['excl_epoch']} |")
(HERE / "select.md").write_text("\n".join(L) + "\n")
(HERE / "select.json").write_text(json.dumps(rows, indent=1))
print("\n".join(L))

"""E28 (2), test side: for each INDEL release rerun whose two selections keep different epochs, the chr1 (and HCC1395T
chr20) INDEL result of both checkpoints: benchmark regions, after the PoN, PASS ppr8 (post-PoN validation recall 0.8 of
the checkpoint). 'with' = best.pth (the original selection; for a bit-identical rerun, the original run's own
vcf_ppr8_* results), 'without' = best_excl_offpath.pth. Writes select_eval.md / select_eval.json here."""
import json
from pathlib import Path
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
HERE = Path(__file__).resolve().parent
sel = json.load(open(HERE / "select.json"))
def rep(d):
    f = d / "eval_bed/report.json"
    if not f.exists():
        return None
    t = json.load(open(f))["pon"]["at_threshold"]
    return dict(P=t["precision"], R=t["recall"], F1=t["f1"], TP=t["tp_call"], FP=t["fp"])
rows = []
for run, r in sel.items():
    if not r.get("complete") or r["sel_epoch"] == r["excl_epoch"]:
        continue
    p = run.split("_")[2]; co = "COLO829T_fiberseq" if p == "PacBio" else f"COLO829T_{p}"
    ident = r["max_dloss"] == 0 and r["max_dsel"] == 0
    for s, c in ((f"HCC1395T_{p}_d46", "chr1"), (f"HG008T_{p}_d46", "chr1"), (f"{co}_d46", "chr1"), (f"HCC1395T_{p}_d46", "chr20")):
        a = rep(R / run / f"vcf_ppr8_{c}_{s}") if ident else rep(R / "selcheck_20261101" / run / f"vcf_ppr8_best_{c}_{s}")
        b = rep(R / "selcheck_20261101" / run / f"vcf_ppr8_excl_{c}_{s}")
        rows.append(dict(run=run, identical=ident, epochs=(r["sel_epoch"], r["excl_epoch"]), set=s, chrom=c, with_=a, without=b))
f = lambda x: "-" if not x else f"{x['P']:.3f} / {x['R']:.3f} / {x['F1']:.3f} ({x['TP']} TP, {x['FP']} FP)"
L = ["| run (epoch with / without off-path) | test | with off-path negatives: P / R / F1 | without: P / R / F1 | Δ F1 |", "|---|---|---|---|---|"]
for x in rows:
    d = "-" if not (x["with_"] and x["without"]) else f"{x['without']['F1'] - x['with_']['F1']:+.3f}"
    L.append(f"| {x['run']}{' (bit-identical rerun)' if x['identical'] else ''} ({x['epochs'][0]} / {x['epochs'][1]}) | {x['chrom']} {x['set']} | {f(x['with_'])} | {f(x['without'])} | {d} |")
(HERE / "select_eval.md").write_text("\n".join(L) + "\n")
(HERE / "select_eval.json").write_text(json.dumps(rows, indent=1))
print("\n".join(L))

# Reproduction: reruns that are not bit-identical but keep the reported epoch under the original rule; test F1 of the
# rerun's best.pth (vcf_ppr8_best_*) against the reported model's (vcf_ppr8_* of the release run). -> select_repro.md
L = ["| run (epoch) | test | reported F1 | rerun F1 | Δ F1 |", "|---|---|---|---|---|"]
for run, r in sel.items():
    if not r.get("complete") or (r["max_dloss"] == 0 and r["max_dsel"] == 0) or r["sel_epoch"] != r["orig_sel_epoch"]:
        continue
    p = run.split("_")[2]; co = "COLO829T_fiberseq" if p == "PacBio" else f"COLO829T_{p}"
    for s, c in ((f"HCC1395T_{p}_d46", "chr1"), (f"HG008T_{p}_d46", "chr1"), (f"{co}_d46", "chr1"), (f"HCC1395T_{p}_d46", "chr20")):
        a = rep(R / run / f"vcf_ppr8_{c}_{s}"); b = rep(R / "selcheck_20261101" / run / f"vcf_ppr8_best_{c}_{s}")
        d = "-" if not (a and b) else f"{b['F1'] - a['F1']:+.3f}"
        L.append(f"| {run} ({r['sel_epoch']}) | {c} {s} | {a['F1'] if a else '-'} | {b['F1'] if b else '-'} | {d} |")
(HERE / "select_repro.md").write_text("\n".join(L) + "\n")
print("\n".join(L))

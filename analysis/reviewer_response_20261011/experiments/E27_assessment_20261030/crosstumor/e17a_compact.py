"""Compact E17a strata (final INDEL release models, ppr8, chr1, BED, after the PoN): per dataset x call set x stratum
truth N, TP (seed mean), FP (seed mean, call-side class), FN (seed mean), P (pooled), R (seed-pooled), F1, from
analysis/.../E17a_indel_strata/final_models.tsv."""
import csv
from collections import defaultdict
T = "/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E17a_indel_strata/final_models.tsv"
acc = defaultdict(lambda: [0, 0, 0, 0, None, set()])  # tpb, fn, tpc, fp, pon_tagged, runs
for r in csv.DictReader(open(T), delimiter="\t"):
    g = r["group"]
    if g not in ("release", "unseen", "ClairS-TO", "DeepSomatic-TO"):
        continue
    key = (r["dataset"], g, r["axis"], r["cls"])
    a = acc[key]
    a[0] += int(r["tp_baseline"]); a[1] += int(r["fn"]); a[2] += int(r["tp_call"]); a[3] += int(r["fp"])
    a[5].add(r["uid"])
    if r["truth_pon_tagged"] not in ("", "None"): a[4] = int(r["truth_pon_tagged"])
STRATA = [("all", "all"), ("type", "INS"), ("type", "DEL"), ("length", "1"), ("length", "2-4"), ("length", "5-10"),
          ("length", ">10"), ("context", "HP>=7"), ("context", "HP4-6"), ("context", "STR"), ("context", "none")]
DS = ["HG008T_Illumina", "HG008T_PacBio", "HG008T_ONT", "COLO829T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT"]
for g in ("release", "unseen", "DeepSomatic-TO", "ClairS-TO"):
    print(f"\n### {g}\n")
    print("| dataset | " + " | ".join(f"{a}:{c}" for a, c in STRATA) + " |")
    print("|---|" + "---|" * len(STRATA))
    for ds in DS:
        cells = []
        for ax, c in STRATA:
            a = acc.get((ds, g, ax, c))
            if not a:
                cells.append("-"); continue
            n = len(a[5]); tpb, fn, tpc, fp = a[0] / n, a[1] / n, a[2] / n, a[3] / n
            truth = tpb + fn
            P = tpc / (tpc + fp) if tpc + fp else float("nan"); Rr = tpb / truth if truth else float("nan")
            F = 2 * P * Rr / (P + Rr) if P == P and Rr == Rr and P + Rr else 0.0
            Ps = f"{P:.2f}" if P == P else "-"
            cells.append(f"N{truth:.0f}" + (f" (PoN {a[4]})" if a[4] is not None else "") + f" TP{tpb:.1f} FP{fp:.1f} FN{fn:.1f} P{Ps} R{Rr:.2f} F{F:.2f}")
        print(f"| {ds} | " + " | ".join(cells) + " |")

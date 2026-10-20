"""E16 per SMaHT stratum (R1.3): block-bootstrap intervals of chr1 SNV P / R / F1 after the PoN inside each SMaHT v2
stratum (easy / difficult / extreme), in the benchmark BED, and paired HPRC - GRCh38-only differences per stratum.

Same blocks, resamples and seed as e16_bootstrap.py (1-Mb blocks of the set's BED, 2,000 resamples, seed 20261015), so a
stratum's interval resamples the same blocks as the call set's overall interval; a record counts in the stratum that
contains the first base of its normalized reference interval (E10's per-record table,
`call_level/E10/<run>__<vcf dir>.v2.tsv.gz`, eval_bed rows). In-sample SNV call sets at the main threshold only.
Writes strata_ci.tsv, strata_paired.tsv and strata.md next to this script.
usage: /wanglab/jshen/anaconda3/bin/python e16_strata.py
"""
import csv
import gzip
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import callsets as C  # noqa: E402
import e16_bootstrap as E  # noqa: E402

E10 = Path("/scratch/jshen/data/reviewer_experiments_20261013/call_level/E10")
STRATA = ("easy", "difficult", "extreme")
PARTS = ("tp-baseline", "fn", "tp", "fp")
WIDTH = E.WIDTHS[0]


def stratum_positions(unit):
    """{stratum: {rtg part: [POS]}} from E10's per-record rows of the unit (eval_bed), or None."""
    run, vdir = unit.uid.split("/")
    f = E10 / f"{run}__{vdir}.v2.tsv.gz"
    if not f.exists():
        return None
    out = {s: {p: [] for p in PARTS} for s in STRATA}
    with gzip.open(f, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["eval"] == "eval_bed" and r["smaht"] in out and r["part"] in PARTS:
                out[r["smaht"]][r["part"]].append(int(r["pos"]))
    return out


def main():
    units = [u for u in C.discover(kinds=("SNV",), comparators=False)
             if u.role == "insample" and u.variant == "standard" and "bed" in u.rtg]
    bs = E.BlockSets(C.fasta_path(units), 2000, 20261015)
    rows, boot, point, missing = [], {}, {}, []
    for u in units:
        sp = stratum_positions(u)
        if sp is None:
            missing.append(u.uid)
            continue
        blocks = bs.get("bed", u.bed, WIDTH)
        for s in STRATA:
            m = E.block_counts(sp[s], blocks)
            tot = m.sum(axis=0)
            res = blocks.mult @ m
            boot[(u.uid, s)], point[(u.uid, s)] = (blocks.key, res), tot
            p, r, f = E.metrics(tot)
            bp, br, bf = E.metrics(res)
            row = dict(uid=u.uid, set=u.base_set, arm=u.arm, stratum=s, tp_baseline=int(tot[0]), fn=int(tot[1]),
                       tp_call=int(tot[2]), fp=int(tot[3]), pass_calls=int(tot[2] + tot[3]))
            for name, val, b in (("P", p, bp), ("R", r, br), ("F1", f, bf)):
                lo, hi, _, _ = E.ci(b)
                row.update({name: float(val), f"{name}_lo": lo, f"{name}_hi": hi})
            rows.append(row)
    prow = []
    by = {}
    for u in units:
        by.setdefault(u.base_set, {})[u.arm] = u
    for bset, arms in sorted(by.items()):
        lin = arms.get("linear")
        if lin is None:
            continue
        for arm in ("d9", "d46"):
            a = arms.get(arm)
            if a is None or (a.truth, a.bed) != (lin.truth, lin.bed):
                continue
            for s in STRATA:
                ka, kb = (a.uid, s), (lin.uid, s)
                if ka not in boot or kb not in boot:
                    continue
                assert boot[ka][0] == boot[kb][0]
                row = dict(set=bset, A=a.uid, B=lin.uid, arm=arm, stratum=s)
                for i, name in enumerate(("P", "R", "F1")):
                    d = E.metrics(boot[ka][1])[i] - E.metrics(boot[kb][1])[i]
                    lo, hi, _, _ = E.ci(d)
                    row.update({f"d{name}": float(E.metrics(point[ka])[i] - E.metrics(point[kb])[i]),
                                f"d{name}_lo": lo, f"d{name}_hi": hi})
                prow.append(row)
    E.write_tsv(HERE / "strata_ci.tsv", rows)
    E.write_tsv(HERE / "strata_paired.tsv", prow)
    f3 = lambda x: "-" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}"  # noqa: E731
    sg = lambda x: "-" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:+.3f}"  # noqa: E731
    L = ["# E16 per SMaHT stratum (chr1 SNV, after the PoN, in the benchmark BED; 1-Mb blocks, 2,000 resamples)", "",
         "## Per call set", "", "| set | arm | stratum | PASS calls | TP / FP / FN | P [95% CI] | R [95% CI] | F1 [95% CI] |",
         "|---|---|---|---:|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['set']} | {r['arm']} | {r['stratum']} | {r['pass_calls']:,} | {r['tp_baseline']} / {r['fp']} / {r['fn']} | "
                 f"{f3(r['P'])} [{f3(r['P_lo'])}, {f3(r['P_hi'])}] | {f3(r['R'])} [{f3(r['R_lo'])}, {f3(r['R_hi'])}] | "
                 f"{f3(r['F1'])} [{f3(r['F1_lo'])}, {f3(r['F1_hi'])}] |")
    L += ["", "## HPRC - GRCh38-only, paired (same resampled blocks)", "",
          "| set | arm | stratum | dP [95% CI] | dR [95% CI] | dF1 [95% CI] |", "|---|---|---|---|---|---|"]
    for r in prow:
        L.append(f"| {r['set']} | {r['arm']} | {r['stratum']} | {sg(r['dP'])} [{sg(r['dP_lo'])}, {sg(r['dP_hi'])}] | "
                 f"{sg(r['dR'])} [{sg(r['dR_lo'])}, {sg(r['dR_hi'])}] | {sg(r['dF1'])} [{sg(r['dF1_lo'])}, {sg(r['dF1_hi'])}] |")
    if missing:
        L += ["", "No E10 per-record table (not included): " + ", ".join(missing)]
    (HERE / "strata.md").write_text("\n".join(L) + "\n")
    print(f"{len(rows)} stratum rows, {len(prow)} paired rows, missing E10 tables: {len(missing)}")


if __name__ == "__main__":
    main()

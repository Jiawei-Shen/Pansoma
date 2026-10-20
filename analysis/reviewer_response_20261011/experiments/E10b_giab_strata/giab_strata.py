"""E10b (R1.3): chr1 P / R / F1 after the PoN inside GIAB v3.5 GRCh38 stratifications, per call set, and paired
HPRC - GRCh38-only (linear graph) differences with block-bootstrap intervals.

Strata: GIAB genome stratifications v3.5, GRCh38 (downloaded 2026-10-20 with the user's approval, md5-checked:
/scratch/jshen/data/GIAB_stratifications_v3.5/GRCh38@all/), the 12 of STRATA below; strata overlap, so each is scored
on its own (a record counts in every stratum that holds its position).
Records: E10's per-record tables (`call_level/E10/<run>__<vcf dir>.v2.tsv.gz`, eval_bed rows: rtg tp-baseline / fn /
tp / fp of the PASS calls after the PoN, in the benchmark BED); a record is in a stratum when its VCF position
(pos - 1, 0-based) lies in one of the stratum's chr1 intervals.
Bootstrap: E16's blocks, resamples and seed (1-Mb blocks of the set's BED, 2,000 resamples, seed 20261015), so the
paired difference resamples the same blocks for both arms. In-sample call sets at the main threshold (validation
recall 0.9), SNV and INDEL, as e16_strata.py.
Writes strata.md, strata_ci.tsv, strata_paired.tsv and summary.json next to this script.
usage: /wanglab/jshen/anaconda3/bin/python giab_strata.py
"""
import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "E16_bootstrap"))
import callsets as C  # noqa: E402
import e16_bootstrap as E  # noqa: E402

E10 = Path("/scratch/jshen/data/reviewer_experiments_20261013/call_level/E10")
GIAB = Path("/scratch/jshen/data/GIAB_stratifications_v3.5/GRCh38@all")
STRATA = [
    ("notinalldifficultregions", "Union/GRCh38_notinalldifficultregions.bed.gz", "not in any difficult region"),
    ("alldifficultregions", "Union/GRCh38_alldifficultregions.bed.gz", "all difficult regions"),
    ("alllowmapandsegdupregions", "Union/GRCh38_alllowmapandsegdupregions.bed.gz", "low mappability or segdup"),
    ("lowmappabilityall", "Mappability/GRCh38_lowmappabilityall.bed.gz", "low mappability"),
    ("segdups", "SegmentalDuplications/GRCh38_segdups.bed.gz", "segmental duplications"),
    ("AllTandemRepeatsandHomopolymers_slop5", "LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
     "tandem repeats and homopolymers"),
    ("AllHomopolymers_ge7bp_imperfectge11bp_slop5", "LowComplexity/GRCh38_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz",
     "homopolymers >= 7 bp"),
    ("AllTandemRepeats_le50bp_slop5", "LowComplexity/GRCh38_AllTandemRepeats_le50bp_slop5.bed.gz", "tandem repeats <= 50 bp"),
    ("satellites_slop5", "LowComplexity/GRCh38_satellites_slop5.bed.gz", "satellites"),
    ("L1H_gt500", "OtherDifficult/GRCh38_L1H_gt500.bed.gz", "L1H > 500 bp"),
    ("gclt25orgt65_slop50", "GCcontent/GRCh38_gclt25orgt65_slop50.bed.gz", "GC < 25% or > 65%"),
    ("CMRGv1.00_falselyduplicatedgenes", "FunctionalTechnicallyDifficult/GRCh38_CMRGv1.00_falselyduplicatedgenes.bed.gz",
     "CMRG falsely duplicated genes"),
]
PARTS = ("tp-baseline", "fn", "tp", "fp")
WIDTH = E.WIDTHS[0]


def chr1_intervals(path):
    """Merged chr1 intervals (starts, ends) of a bed.gz."""
    iv = []
    with gzip.open(path, "rt") as f:
        for line in f:
            c = line.split("\t", 3)
            if c[0] == "chr1":
                iv.append((int(c[1]), int(c[2])))
    iv.sort()
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    a = np.array(merged, dtype=np.int64).reshape(-1, 2)
    return a[:, 0], a[:, 1]


def inside(pos1, iv):
    starts, ends = iv
    p0 = np.asarray(pos1, dtype=np.int64) - 1
    i = np.searchsorted(starts, p0, side="right") - 1
    ok = i >= 0
    out = np.zeros(p0.shape, dtype=bool)
    out[ok] = p0[ok] < ends[i[ok]]
    return out


def records(unit):
    """{rtg part: np.array of positions} from E10's per-record table (eval_bed rows), or None."""
    run, vdir = unit.uid.split("/")
    f = E10 / f"{run}__{vdir}.v2.tsv.gz"
    if not f.exists():
        return None
    out = {p: [] for p in PARTS}
    with gzip.open(f, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["eval"] == "eval_bed" and r["part"] in out:
                out[r["part"]].append(int(r["pos"]))
    return {p: np.array(v, dtype=np.int64) for p, v in out.items()}


def main():
    ivs = {name: chr1_intervals(GIAB / rel) for name, rel, _ in STRATA}
    sizes = {name: int((iv[1] - iv[0]).sum()) for name, iv in ivs.items()}
    units = [u for u in C.discover(kinds=("SNV", "INDEL"), comparators=False)
             if u.role == "insample" and u.variant == "standard" and "bed" in u.rtg and not u.normds]
    bs = E.BlockSets(C.fasta_path(units), 2000, 20261015)
    rows, boot, point, missing = [], {}, {}, []
    for u in units:
        rec = records(u)
        if rec is None:
            missing.append(u.uid)
            continue
        blocks = bs.get("bed", u.bed, WIDTH)
        for name, _, _ in STRATA:
            pos = {p: rec[p][inside(rec[p], ivs[name])].tolist() for p in PARTS}
            m = E.block_counts(pos, blocks)
            tot = m.sum(axis=0)
            res = blocks.mult @ m
            boot[(u.uid, name)], point[(u.uid, name)] = (blocks.key, res), tot
            p, r, f = E.metrics(tot)
            bp, br, bf = E.metrics(res)
            row = dict(uid=u.uid, set=u.base_set, kind=u.kind, arm=u.arm, stratum=name, tp_baseline=int(tot[0]),
                       fn=int(tot[1]), tp_call=int(tot[2]), fp=int(tot[3]), pass_calls=int(tot[2] + tot[3]))
            for nm, val, b in (("P", p, bp), ("R", r, br), ("F1", f, bf)):
                lo, hi, _, _ = E.ci(b)
                row.update({nm: float(val), f"{nm}_lo": lo, f"{nm}_hi": hi})
            rows.append(row)
    by = {}
    for u in units:
        by.setdefault((u.base_set, u.kind), {})[u.arm] = u
    prow = []
    for (bset, kind), arms in sorted(by.items()):
        lin = arms.get("linear")
        if lin is None:
            continue
        for arm in ("d9", "d46"):
            a = arms.get(arm)
            if a is None or (a.truth, a.bed) != (lin.truth, lin.bed):
                continue
            for name, _, _ in STRATA:
                ka, kb = (a.uid, name), (lin.uid, name)
                if ka not in boot or kb not in boot:
                    continue
                assert boot[ka][0] == boot[kb][0]
                row = dict(set=bset, kind=kind, arm=arm, stratum=name, A=a.uid, B=lin.uid,
                           truth=int(point[ka][0] + point[ka][1]))
                for i, nm in enumerate(("P", "R", "F1")):
                    d = E.metrics(boot[ka][1])[i] - E.metrics(boot[kb][1])[i]
                    lo, hi, _, _ = E.ci(d)
                    row.update({f"d{nm}": float(E.metrics(point[ka])[i] - E.metrics(point[kb])[i]),
                                f"d{nm}_lo": lo, f"d{nm}_hi": hi})
                prow.append(row)
    E.write_tsv(HERE / "strata_ci.tsv", rows)
    E.write_tsv(HERE / "strata_paired.tsv", prow)

    f3 = lambda x: "-" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}"  # noqa: E731
    sg = lambda x: "-" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:+.3f}"  # noqa: E731
    label = {n: l for n, _, l in STRATA}
    summ = []
    L = ["# E10b: GIAB v3.5 stratifications (chr1, PASS after the PoN, in the benchmark BED; 1-Mb blocks, 2,000 resamples)",
         "", "## 1. Summary: HPRC graph - GRCh38-only (linear) graph, paired F1 difference", "",
         "Per stratum and variant kind: the 12 pairs (HG008T / COLO829T x Illumina / PacBio-fiberseq / ONT x v1.1 d9 / "
         "v2.1 d46, where both arms exist). Truth = truth records of the stratum in the BED (range over the six sets).", "",
         "| stratum | chr1 Mb | kind | pairs | truth per set | dF1 range | dF1 > 0 | CI above 0 | CI below 0 |",
         "|---|---:|---|---:|---|---|---:|---:|---:|"]
    for kind in ("SNV", "INDEL"):
        for name, _, _ in STRATA:
            pr = [r for r in prow if r["kind"] == kind and r["stratum"] == name and not np.isnan(r["dF1"])]
            if not pr:
                continue
            tr = [r["truth"] for r in pr]
            d = [r["dF1"] for r in pr]
            s = dict(stratum=name, label=label[name], kind=kind, chr1_mb=sizes[name] / 1e6, pairs=len(pr),
                     truth_min=min(tr), truth_max=max(tr), dF1_min=min(d), dF1_max=max(d),
                     positive=sum(x > 0 for x in d), ci_above=sum(r["dF1_lo"] > 0 for r in pr),
                     ci_below=sum(r["dF1_hi"] < 0 for r in pr))
            summ.append(s)
            L.append(f"| {label[name]} | {sizes[name] / 1e6:.1f} | {kind} | {s['pairs']} | {s['truth_min']}-{s['truth_max']} | "
                     f"{sg(s['dF1_min'])} to {sg(s['dF1_max'])} | {s['positive']} | {s['ci_above']} | {s['ci_below']} |")
    L += ["", "## 2. Paired differences (HPRC - linear graph)", "",
          "| set | kind | arm | stratum | truth | dP [95% CI] | dR [95% CI] | dF1 [95% CI] |", "|---|---|---|---|---:|---|---|---|"]
    for r in prow:
        L.append(f"| {r['set']} | {r['kind']} | {r['arm']} | {label[r['stratum']]} | {r['truth']} | "
                 f"{sg(r['dP'])} [{sg(r['dP_lo'])}, {sg(r['dP_hi'])}] | {sg(r['dR'])} [{sg(r['dR_lo'])}, {sg(r['dR_hi'])}] | "
                 f"{sg(r['dF1'])} [{sg(r['dF1_lo'])}, {sg(r['dF1_hi'])}] |")
    L += ["", "## 3. Per call set", "",
          "| set | kind | arm | stratum | PASS calls | TP / FP / FN | P [95% CI] | R [95% CI] | F1 [95% CI] |",
          "|---|---|---|---|---:|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['set']} | {r['kind']} | {r['arm']} | {label[r['stratum']]} | {r['pass_calls']:,} | "
                 f"{r['tp_baseline']} / {r['fp']} / {r['fn']} | {f3(r['P'])} [{f3(r['P_lo'])}, {f3(r['P_hi'])}] | "
                 f"{f3(r['R'])} [{f3(r['R_lo'])}, {f3(r['R_hi'])}] | {f3(r['F1'])} [{f3(r['F1_lo'])}, {f3(r['F1_hi'])}] |")
    if missing:
        L += ["", "No E10 per-record table (not included): " + ", ".join(missing)]
    (HERE / "strata.md").write_text("\n".join(L) + "\n")
    (HERE / "summary.json").write_text(json.dumps(dict(chr1_bp=sizes, summary=summ, units=[u.uid for u in units],
                                                       missing=missing), indent=1) + "\n")
    print("\n".join(L[:6 + len(summ) + 2]))
    print(f"{len(rows)} rows, {len(prow)} paired rows, missing {len(missing)}")


if __name__ == "__main__":
    main()

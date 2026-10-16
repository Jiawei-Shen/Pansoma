#!/usr/bin/env python3
"""E16 (R3.9): block-bootstrap confidence intervals of chr1 P / R / F1 after the PoN, per call set, and paired CIs of
F1 differences (reference arms; cross-sample vs in-sample; --norm dataset vs not; threshold / rescue variants vs the
standard call set; comparators vs Pansoma).

    /wanglab/jshen/anaconda3/bin/python e16_bootstrap.py [--resamples 2000] [--seed 20261015]

Input: every finished unit of callsets.discover() (rtg eval_bed / eval_nobed pon_calls: tp, tp-baseline, fp, fn).
Method:
  - chr1 is cut into fixed blocks of 1 Mb (main) and 5 Mb (sensitivity). Block set of a mode: eval_bed = the blocks
    that overlap the set's benchmark BED on chr1; eval_nobed = the blocks with at least one non-N GRCh38 base.
  - Each rtg record goes to the block of its POS: truth-side tp-baseline / fn and call-side tp / fp. Per block:
    TPb, FN, TPc, FP. P = TPc / (TPc + FP), R = TPb / (TPb + FN), F1 = 2PR / (P + R) (rtg's definitions).
  - 2,000 resamples of n blocks with replacement (n = size of the block set). The multiplicity matrix of a block set is
    drawn once (seed fixed per block set), so every call set on the same BED shares the same resamples: a paired
    difference is computed resample by resample.
  - 95 % percentile intervals; P(diff > 0) = share of resamples with A > B.
  - The point estimates are the totals; they are checked against report.json (at_threshold of 'pon').
Not covered: PASS calls without a GRCh38 position (linear_vcf's unplaced records) are outside rtg and outside these
CIs; seed-to-seed training variance needs repeated GPU training.
Writes ci.tsv, paired.tsv, summary.json, tables.md next to this script.
"""
import argparse
import csv
import json
import sys
import zlib
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import callsets as C  # noqa: E402

WIDTHS = (1_000_000, 5_000_000)
MODES = ("bed", "nobed")
CHROM = "chr1"


def chr1_geometry(fasta):
    """chr1 length and the per-1-Mb count of non-N bases."""
    import pysam
    f = pysam.FastaFile(fasta)
    seq = f.fetch(CHROM).upper()
    n = len(seq)
    step = WIDTHS[0]
    nonn = [step_len - seq.count("N", s, s + step_len)
            for s in range(0, n, step) for step_len in [min(step, n - s)]]
    return n, nonn


def read_bed_chr1(bed):
    out = []
    with open(bed) as f:
        for line in f:
            c = line.split("\t")
            if c[0] == CHROM:
                out.append((int(c[1]), int(c[2])))
    return sorted(out)


class Blocks:
    """Block set of one (mode, BED, width): block ids (window index on chr1) and the resample multiplicities."""

    def __init__(self, ids, width, key, resamples, seed):
        self.ids = list(ids)
        self.index = {b: i for i, b in enumerate(self.ids)}
        self.width = width
        self.key = key
        n = len(self.ids)
        rng = np.random.default_rng([seed, width // 1000, zlib.crc32(key.encode())])
        draws = rng.integers(0, n, size=(resamples, n))
        self.mult = np.zeros((resamples, n), dtype=np.int32)
        for r in range(resamples):
            self.mult[r] = np.bincount(draws[r], minlength=n)
        self.extra = 0  # records outside the block set (should stay 0)

    def block_of(self, pos):
        return (pos - 1) // self.width


class BlockSets:
    def __init__(self, fasta, resamples, seed):
        self.length, self.nonn = chr1_geometry(fasta)
        self.cache = {}
        self.resamples, self.seed = resamples, seed

    def get(self, mode, bed, width):
        key = f"{mode}|{bed if mode == 'bed' else '-'}|{width}"
        if key not in self.cache:
            nwin = (self.length + width - 1) // width
            if mode == "bed":
                ids = set()
                for s, e in read_bed_chr1(bed):
                    ids.update(range(s // width, (e - 1) // width + 1))
            else:
                per = width // WIDTHS[0]
                ids = {w for w in range(nwin) if sum(self.nonn[w * per:(w + 1) * per]) > 0}
            self.cache[key] = Blocks(sorted(ids), width, key, self.resamples, self.seed)
        return self.cache[key]


def positions(unit, mode):
    """{rtg file: [POS]} of chr1 records of the unit's pon_calls directory in this mode."""
    w = unit.rtg[mode]
    out = {}
    for k in C.RTG_FILES:
        out[k] = [r[1] for r in C.read_records(w / f"{k}.vcf.gz") if r[0] == CHROM]
    return out


def block_counts(pos, blocks):
    """n_blocks x 4 matrix (TPb, FN, TPc, FP)."""
    m = np.zeros((len(blocks.ids), 4), dtype=np.int64)
    for col, k in enumerate(("tp-baseline", "fn", "tp", "fp")):
        for p in pos[k]:
            b = blocks.block_of(p)
            i = blocks.index.get(b)
            if i is None:
                blocks.extra += 1
                continue
            m[i, col] += 1
    return m


def metrics(tot):
    """P, R, F1 arrays from totals (..., 4)."""
    tot = np.asarray(tot, dtype=float)
    tpb, fn, tpc, fp = tot[..., 0], tot[..., 1], tot[..., 2], tot[..., 3]
    with np.errstate(invalid="ignore", divide="ignore"):
        p = tpc / (tpc + fp)
        r = tpb / (tpb + fn)
        f = 2 * p * r / (p + r)
    return p, r, f


def ci(x):
    x = np.asarray(x, dtype=float)
    ok = x[~np.isnan(x)]
    if ok.size == 0:
        return (np.nan, np.nan, np.nan, 0)
    lo, hi = np.percentile(ok, [2.5, 97.5])
    return float(lo), float(hi), float(ok.std(ddof=1)) if ok.size > 1 else np.nan, int(x.size - ok.size)


# Known recipe differences of the d9 checkpoints (E3_reference_arms/README.md, caveat 1): pairs with these runs are
# not matched arms.
D9_NOTES = {
    "HG008_PacBio_SNV_7fc7396": "d9 trained on 1 GPU (global batch 1,024), linear / d46 on 4 GPUs (4,096): not a matched arm (E3 caveat 1)",
    "HG008_ONT_SNV_7fc7396": "d9 trained on 1 GPU (global batch 1,024), linear / d46 on 4 GPUs (4,096): not a matched arm (E3 caveat 1)",
}


def pairs(units):
    """(type, A, B, note): A - B differences of F1 on identical blocks."""
    by_uid = {u.uid: u for u in units}
    out = []
    std = [u for u in units if u.variant == "standard"]
    # reference arms: in-sample standard units of one base set and kind
    groups = defaultdict(dict)
    for u in std:
        if u.role == "insample":
            groups[(u.base_set, u.kind)][u.arm] = u
    for (bs, kind), arms in sorted(groups.items()):
        for a, b in (("d9", "linear"), ("d46", "linear"), ("d46", "d9")):
            if a in arms and b in arms:
                note = "; ".join(D9_NOTES[u.run] for u in (arms[a], arms[b]) if u.run in D9_NOTES)
                out.append(("arms", arms[a], arms[b], note))
        joint = [u for u in std if u.role == "joint" and u.set == bs and u.kind == kind]
        for j in joint:
            for b in ("linear", "d46"):
                if b in arms:
                    x, y = (j, arms[b]) if b == "linear" else (arms[b], j)
                    out.append(("arms (joint d9)", x, y, "d9 = joint model trained on 4 purities (not a matched arm)"))
    # cross-sample vs in-sample (same target set and kind)
    insample = {(u.set, u.kind): u for u in std if u.role == "insample"}
    for u in std:
        if u.role == "xs" and (u.set, u.kind) in insample:
            out.append(("cross-sample vs in-sample", u, insample[(u.set, u.kind)], ""))
    # --norm dataset (normds) vs the old checkpoints' prediction
    old = {(u.model, u.set, u.kind): u for u in std if u.role == "xs" and not u.normds}
    for u in std:
        if u.role == "xs" and u.normds and (u.model, u.set, u.kind) in old:
            out.append(("normds vs old xs", u, old[(u.model, u.set, u.kind)], ""))
    # variants vs the standard call set of the same predictions
    for u in units:
        if u.variant != "standard" and u.standard_uid in by_uid:
            out.append((f"{u.variant} vs standard", u, by_uid[u.standard_uid], ""))
    # comparators vs Pansoma in-sample units of the same base set and kind
    for c in units:
        if c.role == "comparator":
            for u in std:
                if u.role in ("insample", "joint") and u.base_set == C.base_set(c.set) and u.kind == c.kind:
                    out.append(("Pansoma vs comparator", u, c, ""))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--resamples", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261015)
    args = ap.parse_args()
    units = C.discover()
    if not units:
        raise SystemExit("no finished call sets")
    bs = BlockSets(C.fasta_path(units), args.resamples, args.seed)
    rows, check = [], []
    boot = {}  # (uid, mode, width) -> (blocks key, resampled totals B x 4)
    point = {}
    for u in units:
        for mode in MODES:
            if mode not in u.rtg:
                continue
            pos = positions(u, mode)
            for width in WIDTHS:
                blocks = bs.get(mode, u.bed, width)
                m = block_counts(pos, blocks)
                tot = m.sum(axis=0)
                res = blocks.mult @ m  # B x 4
                boot[(u.uid, mode, width)] = (blocks.key, res)
                point[(u.uid, mode, width)] = tot
                p, r, f = metrics(tot)
                bp, br, bf = metrics(res)
                row = dict(uid=u.uid, label=u.label, role=u.role, variant=u.variant, set=u.set, kind=u.kind,
                           arm=u.arm, model=u.model, normds=u.normds, mode=mode, block_mb=width // 1_000_000,
                           n_blocks=len(blocks.ids), tp_baseline=int(tot[0]), fn=int(tot[1]), tp_call=int(tot[2]),
                           fp=int(tot[3]), unplaced_pass=u.unplaced_pass)
                for name, val, b in (("P", p, bp), ("R", r, br), ("F1", f, bf)):
                    lo, hi, se, nan = ci(b)
                    row.update({name: float(val), f"{name}_lo": lo, f"{name}_hi": hi, f"{name}_width": hi - lo,
                                f"{name}_se": se, f"{name}_nan": nan})
                rows.append(row)
                if width == WIDTHS[0]:
                    at = (u.report.get(mode) or {}).get("pon", {}).get("at_threshold")
                    if at:
                        same = (at["tp_baseline"], at["fn"], at["tp_call"], at["fp"]) == tuple(int(x) for x in tot)
                        check.append(dict(uid=u.uid, mode=mode, ok=same, rtg=[int(x) for x in tot],
                                          report=[at["tp_baseline"], at["fn"], at["tp_call"], at["fp"]],
                                          f1_report=at["f1"]))
    # paired differences
    prow = []
    for typ, a, b, note in pairs(units):
        if (a.truth, a.bed) != (b.truth, b.bed):
            prow.append(dict(type=typ, A=a.uid, B=b.uid, A_label=a.label, B_label=b.label, mode="-", block_mb="-",
                             note=f"skipped: truth / BED differ ({a.truth} {a.bed} vs {b.truth} {b.bed})"))
            continue
        for mode in MODES:
            for width in WIDTHS:
                ka, kb = (a.uid, mode, width), (b.uid, mode, width)
                if ka not in boot or kb not in boot:
                    continue
                assert boot[ka][0] == boot[kb][0], (ka, kb)
                pa, ra, fa = metrics(boot[ka][1])
                pb, rb, fb = metrics(boot[kb][1])
                PA, RA, FA = metrics(point[ka])
                PB, RB, FB = metrics(point[kb])
                row = dict(type=typ, A=a.uid, B=b.uid, A_label=a.label, B_label=b.label, mode=mode,
                           block_mb=width // 1_000_000, note=note, F1_A=float(FA), F1_B=float(FB))
                for name, da, db, va, vb in (("F1", fa, fb, FA, FB), ("P", pa, pb, PA, PB), ("R", ra, rb, RA, RB)):
                    d = da - db
                    lo, hi, se, nan = ci(d)
                    ok = d[~np.isnan(d)]
                    row.update({f"d{name}": float(va - vb), f"d{name}_lo": lo, f"d{name}_hi": hi,
                                f"d{name}_se": se, f"d{name}_p_gt0": float((ok > 0).mean()) if ok.size else np.nan})
                prow.append(row)
    write_tsv(HERE / "ci.tsv", rows)
    write_tsv(HERE / "paired.tsv", prow)
    extra = {k: b.extra for k, b in bs.cache.items()}
    summary = dict(generated=datetime.now().isoformat(timespec="seconds"), resamples=args.resamples, seed=args.seed,
                   block_widths=WIDTHS, units=[u.uid for u in units],
                   block_sets={k: len(b.ids) for k, b in bs.cache.items()}, records_outside_blocks=extra,
                   report_check=check, ci=rows, paired=prow)
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    (HERE / "tables.md").write_text(tables(units, rows, prow, check, bs, args) + "\n")
    bad = [c for c in check if not c["ok"]]
    print(f"{len(units)} call sets, {len(rows)} CI rows, {len(prow)} paired rows; report check: "
          f"{len(check) - len(bad)} / {len(check)} match; records outside blocks: {sum(extra.values())}")
    for c in bad:
        print("MISMATCH", c)


def write_tsv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def cell(r, name):
    return f"{C.fmt(r[name])} [{C.fmt(r[name + '_lo'])}, {C.fmt(r[name + '_hi'])}]"


def dcell(r, name="F1"):
    return (f"{r['d' + name]:+.3f} [{r['d' + name + '_lo']:+.3f}, {r['d' + name + '_hi']:+.3f}]"
            if r.get("d" + name) == r.get("d" + name) and r.get("d" + name) is not None else "-")


def tables(units, rows, prow, check, bs, args):
    idx = {(r["uid"], r["mode"], r["block_mb"]): r for r in rows}
    lines = ["# E16: block-bootstrap CIs (chr1, PASS after the PoN)", "",
             f"Generated by `e16_bootstrap.py` ({datetime.now():%Y-%m-%d %H:%M}); {len(units)} call sets. "
             f"{args.resamples:,} resamples of 1-Mb blocks (main) and 5-Mb blocks (sensitivity); seed {args.seed}. "
             "95 % percentile intervals. Point values = rtg's totals (checked against report.json). Unplaced PASS "
             "calls (no GRCh38 position) are outside rtg and outside every interval; one training seed per model.", "",
             "Block sets: " + "; ".join(f"`{k}` {len(b.ids)} blocks" for k, b in sorted(bs.cache.items())), ""]
    for mode, title in (("bed", "inside the benchmark BED"), ("nobed", "whole chr1 (no BED)")):
        lines += [f"## 1{'a' if mode == 'bed' else 'b'}. Per call set, {title}", "",
                  "| call set | role | TP / FP / FN | P [95 % CI] | R [95 % CI] | F1 [95 % CI] (1 Mb) | F1 CI (5 Mb) | "
                  "unplaced PASS (not scored) |", "|---|---|---|---|---|---|---|---:|"]
        for u in units:
            r = idx.get((u.uid, mode, 1))
            if not r:
                continue
            r5 = idx.get((u.uid, mode, 5))
            c5 = f"[{C.fmt(r5['F1_lo'])}, {C.fmt(r5['F1_hi'])}]" if r5 else "-"
            lines.append(f"| {u.label} | {u.role} | {r['tp_baseline']} / {r['fp']} / {r['fn']} | {cell(r, 'P')} | "
                         f"{cell(r, 'R')} | {cell(r, 'F1')} | {c5} | "
                         f"{u.unplaced_pass if u.unplaced_pass is not None else '-'} |")
        lines.append("")
    groups = defaultdict(list)
    for r in prow:
        groups[r["type"]].append(r)
    order = ["arms", "arms (joint d9)", "cross-sample vs in-sample", "normds vs old xs", "valbest vs standard",
             "norescue vs standard", "Pansoma vs comparator"]
    n = 2
    for typ in order + sorted(set(groups) - set(order)):
        if typ not in groups:
            continue
        rs = groups[typ]
        lines += [f"## {n}. Paired F1 differences: {typ} (A - B, same resampled blocks)", "",
                  "| A | B | mode | F1 A | F1 B | dF1 [95 % CI] 1 Mb | P(dF1 > 0) | dF1 CI 5 Mb | dP [CI] | dR [CI] |",
                  "|---|---|---|---:|---:|---|---:|---|---|---|"]
        by = defaultdict(dict)
        for r in rs:
            by[(r["A"], r["B"], r["mode"])][r["block_mb"]] = r
        for (a, b, mode), d in by.items():
            r, r5 = d.get(1), d.get(5)
            if r is None:
                rr = d.get("-")
                if rr:
                    lines.append(f"| {rr['A_label']} | {rr['B_label']} | - | - | - | {rr['note']} | | | | |")
                continue
            c5 = f"[{r5['dF1_lo']:+.3f}, {r5['dF1_hi']:+.3f}]" if r5 else "-"
            if all(r[k] == 0 for k in ("dF1", "dF1_lo", "dF1_hi", "dP_lo", "dP_hi", "dR_lo", "dR_hi")):
                lines.append(f"| {r['A_label']} | {r['B_label']}{' (*)' if r.get('note') else ''} | {mode} | {C.fmt(r['F1_A'])} | {C.fmt(r['F1_B'])} | "
                             "identical rtg call sets (0 in every resample) | - | - | - | - |")
                continue
            lines.append(f"| {r['A_label']} | {r['B_label']}{' (*)' if r.get('note') else ''} | {mode} | {C.fmt(r['F1_A'])} | {C.fmt(r['F1_B'])} | "
                         f"{dcell(r)} | {C.fmt(r['dF1_p_gt0'])} | {c5} | {dcell(r, 'P')} | {dcell(r, 'R')} |")
        notes = sorted({r["note"] for r in rs if r.get("note") and r["mode"] != "-"})
        lines += [""] + [f"(*) {x}" for x in notes] + ([""] if notes else [])
        n += 1
    if not any(t in groups for t in ("arms",)):
        lines += [f"## {n}. Paired F1 differences: arms", "",
                  "No sample / platform / kind has two reference arms finished yet (GPU array 405185 is producing the "
                  "linear / d9 / d46 in-sample runs). Rerun `collect.sh`.", ""]
        n += 1
    lines += [f"## {n}. Block-size sensitivity (F1 CI width, 5-Mb blocks / 1-Mb blocks)", "",
              "| mode | kind | role | call sets | F1 CI width 1 Mb (min-max) | width 5 Mb (min-max) | ratio 5 / 1 Mb "
              "(median; min-max) |", "|---|---|---|---:|---|---|---|"]
    groups = defaultdict(list)
    for u in units:
        for mode in MODES:
            r1, r5 = idx.get((u.uid, mode, 1)), idx.get((u.uid, mode, 5))
            if r1 and r5 and r1["F1_width"] > 0:
                groups[(mode, u.kind, u.role)].append((r1["F1_width"], r5["F1_width"]))
    for (mode, kind, role), v in sorted(groups.items()):
        w1, w5 = [a for a, _ in v], [b for _, b in v]
        ratio = sorted(b / a for a, b in v)
        lines.append(f"| {mode} | {kind} | {role} | {len(v)} | {min(w1):.3f}-{max(w1):.3f} | {min(w5):.3f}-{max(w5):.3f} | "
                     f"{float(np.median(ratio)):.2f}; {ratio[0]:.2f}-{ratio[-1]:.2f} |")
    lines.append("")
    n += 1
    bad = [c for c in check if not c["ok"]]
    lines += [f"## {n}. Checks", "",
              f"- rtg record counts vs report.json `pon.at_threshold` (TPb, FN, TPc, FP): {len(check) - len(bad)} of "
              f"{len(check)} call-set / mode pairs identical." + ("" if not bad else " Mismatches: " +
                                                                    ", ".join(f"{c['uid']} {c['mode']}" for c in bad)),
              f"- Records outside the block set: {sum(b.extra for b in bs.cache.values())}.", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    main()

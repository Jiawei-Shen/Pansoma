"""E10: precision / recall per genomic stratum (R1.3: difficult regions).

GIAB stratifications are not local; local strata only:
  smaht    SMaHT v2 easy / difficult / extreme (GRCh38, /scratch/qfu/reference/hg38_DAC_stratification/
           SMaHT_{easy,difficult,extreme}_v2.bed; disjoint, they tile chr1). COLO829T's somatic BED is their union
           (COLO829T_truth/README.txt); for HG008T the strata are applied inside the HG008T benchmark BED.
  bench    HG008T only: inside vs outside the HG008T benchmark BED, from the evaluation without the BED (eval_nobed;
           outside the BED the truth is not complete, so precision there is a lower bound).
  repeat   GRCh38 repeat context of the allele (calls_lib.Repeats = the pipeline's RepeatContext; SNVs by its tract()):
           HP>=7 / STR / other, in the BED.
A record's stratum is the one containing the first base of its normalized reference interval.
Unplaced PASS calls (no GRCh38 position) per stratum: from the E9a projection TSV of the run when it exists
(<E9A_DIR>/<run>.projected.tsv: classes i / ii carry a GRCh38 position, class iii none); otherwise only their number.
Per-record rows: DATA/E10/<call set>.tsv.gz. Tables: tables.md, summary.json. usage: python region_strata.py
"""
import csv
import gzip
import json
import os
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "call_level_common"))
import calls_lib as L  # noqa: E402

OUT = L.DATA / "E10"
SMAHT = ("easy", "difficult", "extreme")
REPEAT = ("HP>=7", "STR", "other")
PARTS = ("tp-baseline", "fn", "tp", "fp")


class Strata:
    def __init__(self):
        self.smaht = {s: L.bed(str(L.SMAHT_DIR / f"SMaHT_{s}_v2.bed")) for s in SMAHT}
        self.rep = {}

    def smaht_of(self, x):
        for s in SMAHT:
            if self.smaht[s].contains_point(x):
                return s
        return "none"

    def repeats(self, fasta):
        if fasta not in self.rep:
            self.rep[fasta] = L.Repeats(fasta)
        return self.rep[fasta]


def process(cs, st):
    out = OUT / f"{cs.id}.v2.tsv.gz"
    srcs = [cs.path(ev, p) for ev in ("eval_bed", "eval_nobed") for p in PARTS]
    if out.exists() and all(out.stat().st_mtime > s.stat().st_mtime for s in srcs):
        return read_rows(out)
    m = L.manifest(cs.set, cs.kind)
    seq = L.chrom_seq(m["fasta"])
    rep = st.repeats(m["fasta"])
    bench = L.bed(m["truth_bed"])
    tk = L.truth_tensor_keys(cs.set, cs.kind)
    rows = []
    for ev in ("eval_bed", "eval_nobed"):
        for part in PARTS:
            for rec in L.read_records(cs.path(ev, part)):
                key = None
                for alt in rec.alts:
                    key = L.normalize(seq, rec.pos, rec.ref, alt)
                    if key:
                        break
                if key is None:
                    continue
                a, b = L.interval(*key)
                rows.append([ev, part, rec.pos, rec.ref, ",".join(rec.alts), st.smaht_of(a),
                             int(bench.overlaps(a, b)), L.repeat_group(rep.classify(rec.pos, rec.ref, rec.alts[0])),
                             rec.info_get("RGN") or "",
                             int(L.has_tensor(tk, rec)) if part in ("tp-baseline", "fn") else ""])
    L.write_tsv_gz(out, ["eval", "part", "pos", "ref", "alt", "smaht", "in_bench_bed", "repeat", "truth_rgn",
                         "has_tensor"], rows)
    return read_rows(out)


def read_rows(path):
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f]


def tally(rows, ev, axis, values):
    new = lambda: dict(tpb=0, fn=0, tpc=0, fp=0, truth_tensor=0, tpb_tensor=0)  # noqa: E731
    d = {v: new() for v in values}
    k = {"tp-baseline": "tpb", "fn": "fn", "tp": "tpc", "fp": "fp"}
    for r in rows:
        if r["eval"] != ev:
            continue
        v = r[axis]
        if axis == "in_bench_bed":
            v = "inside" if v == "1" else "outside"
        x = d.setdefault(v, new())
        x[k[r["part"]]] += 1
        if r.get("has_tensor") == "1":
            x["truth_tensor"] += 1
            x["tpb_tensor"] += r["part"] == "tp-baseline"
    for v, x in d.items():
        p, rr, f = L.prf(x["tpb"], x["fn"], x["tpc"], x["fp"])
        n = x["tpb"] + x["fn"]
        x.update(precision=p, recall=rr, f1=f, truth=n, calls=x["tpc"] + x["fp"],
                 tensor_share=x["truth_tensor"] / n if n else None,
                 recall_with_tensor=x["tpb_tensor"] / x["truth_tensor"] if x["truth_tensor"] else None)
    return d


def unplaced(cs, st, bench):
    """Unplaced PASS calls per SMaHT stratum / benchmark BED, from the E9a projection TSV of the run."""
    n = cs.stats.get("unplaced_pass")
    res = dict(unplaced_pass=n, projected=False)
    tsv = L.E9A_DIR / f"{cs.run}.projected.tsv"
    if cs.variant != "main" or cs.group != "insample" or not tsv.exists():
        return res
    if tsv.stat().st_mtime < cs.unplaced_vcf.stat().st_mtime:
        res["note"] = "projection older than the call set; not used"
        return res
    sm, bb, cl = Counter(), Counter(), Counter()
    with open(tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["filter"] != "PASS":
                continue
            cl[r["proj_class"]] += 1
            if r["proj_class"] in ("i", "ii") and r["grch38_start0"] not in ("", None):
                x = int(r["grch38_start0"])
                sm[st.smaht_of(x)] += 1
                bb["inside" if bench.contains_point(x) else "outside"] += 1
            else:
                sm["no position"] += 1
                bb["no position"] += 1
    res.update(projected=True, tsv=str(tsv), classes=dict(cl), smaht=dict(sm), bench=dict(bb), total=sum(cl.values()))
    return res


def intersect_bp(x, y):
    i = j = 0
    tot = 0
    while i < len(x.starts) and j < len(y.starts):
        a, b = max(x.starts[i], y.starts[j]), min(x.ends[i], y.ends[j])
        if a < b:
            tot += b - a
        if x.ends[i] < y.ends[j]:
            i += 1
        else:
            j += 1
    return tot


def row_cells(d, values):
    cells = []
    for v in values:
        x = d.get(v, {})
        cells += [f"{x.get('truth', 0)} / {x.get('calls', 0)}", x.get("recall"), x.get("precision"), x.get("f1"),
                  f"{L.fmt(x.get('tensor_share'))} / {L.fmt(x.get('recall_with_tensor'))}"]
    return cells


def main():
    t0 = time.time()
    sets = sorted(L.discover(groups=("insample", "joint"), variants=("main", "norescue", "valbest")), key=L.sort_key)
    only = os.environ.get("ONLY")
    if only:
        sets = [cs for cs in sets if only in cs.id]
    L.log(f"{len(sets)} call sets")
    st = Strata()
    summ = []
    for cs in sets:
        rows = process(cs, st)
        m = L.manifest(cs.set, cs.kind)
        bench = L.bed(m["truth_bed"])
        res = dict(id=cs.id, label=cs.label, group=cs.group, variant=cs.variant, sample=cs.sample,
                   platform=cs.platform, kind=cs.kind, arm=cs.arm, purity=cs.purity,
                   smaht=tally(rows, "eval_bed", "smaht", SMAHT), repeat=tally(rows, "eval_bed", "repeat", REPEAT),
                   unplaced=unplaced(cs, st, bench))
        if cs.sample == "HG008T":
            res["bench"] = tally(rows, "eval_nobed", "in_bench_bed", ("inside", "outside"))
        if cs.sample == "COLO829T":
            agree = Counter((r["truth_rgn"].lower(), r["smaht"]) for r in rows
                            if r["eval"] == "eval_bed" and r["part"] in ("tp-baseline", "fn"))
            res["rgn_vs_bed"] = {f"{a} -> {b}": n for (a, b), n in agree.items()}
        summ.append(res)
        L.log(f"{cs.id}: {len(rows)} records")

    out = ["# E10: precision / recall per genomic stratum (chr1)", "",
           f"Generated by `region_strata.py` ({time.strftime('%Y-%m-%d %H:%M')}). PASS calls at the run's threshold "
           "after the PoN (rtg pon_calls). Cells: `truth / calls` (TP + FP), R, P, F1, and `tensor share / R with tensor` "
           "(share of the stratum's truth with an exact tensor in the call set's own tensors, recall among those). "
           "GIAB stratifications are not local; SMaHT v2 easy / difficult / extreme are.", ""]
    # stratum sizes
    hb = L.bed("/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed")
    sizes = []
    for s in SMAHT:
        b = st.smaht[s]
        sizes.append([s, f"{b.bp / 1e6:.1f}", f"{intersect_bp(b, hb) / 1e6:.1f}", L.fmt(intersect_bp(b, hb) / b.bp)])
    out += ["## 1. SMaHT strata on chr1 and the HG008T benchmark BED", ""]
    out += L.md_table(["stratum", "chr1 Mb", "in the HG008T BED, Mb", "share in the HG008T BED"], sizes)
    out += ["", f"HG008T benchmark BED on chr1: {hb.bp / 1e6:.1f} Mb.", ""]
    sec = 2
    sm_head = ["call set"] + [f"{s} {c}" for s in SMAHT for c in ("truth / calls", "R", "P", "F1", "tensor share / R with tensor")]
    for title, pred in (("SMaHT strata, in the BED: in-sample", lambda s: s["group"] == "insample"),
                        ("SMaHT strata: joint model per purity (COLO829T fiberseq)", lambda s: s["group"] == "joint")):
        rs = [s for s in summ if pred(s)]
        out += [f"## {sec}. {title}", ""]
        sec += 1
        out += (L.md_table(sm_head, [[s["label"]] + row_cells(s["smaht"], SMAHT) for s in rs]) if rs else ["(none finished yet)"]) + [""]
    rs = [s for s in summ if "bench" in s]
    out += [f"## {sec}. HG008T: inside vs outside the benchmark BED (evaluation without the BED)", "",
            "Outside the BED the truth is incomplete: FPs there include unverified calls, so P is a lower bound.", ""]
    sec += 1
    out += (L.md_table(["call set"] + [f"{s} {c}" for s in ("inside", "outside") for c in ("truth / calls", "R", "P", "F1", "tensor share / R with tensor")],
                       [[s["label"]] + row_cells(s["bench"], ("inside", "outside")) for s in rs]) if rs else ["(none finished yet)"]) + [""]
    rs = summ
    out += [f"## {sec}. Repeat context (GRCh38; pipeline RepeatContext), in the BED", ""]
    sec += 1
    out += (L.md_table(["call set"] + [f"{s} {c}" for s in REPEAT for c in ("truth / calls", "R", "P", "F1", "tensor share / R with tensor")],
                       [[s["label"]] + row_cells(s["repeat"], REPEAT) for s in rs]) if rs else ["(none finished yet)"]) + [""]
    out += [f"## {sec}. Unplaced PASS calls per stratum (E9a projection; not scored by rtg)", "",
            "Classes i / ii of E9a carry a GRCh38 position (stratum of that position); class iii has none. Runs without "
            "a current E9a projection give only the number of unplaced PASS calls.", ""]
    sec += 1
    urows = []
    for s in summ:
        u = s["unplaced"]
        if u.get("projected"):
            urows.append([s["label"], u["unplaced_pass"], u["total"]] + [u["smaht"].get(x, 0) for x in SMAHT]
                         + [u["smaht"].get("no position", 0), u["bench"].get("inside", 0) if s["sample"] == "HG008T" else "-"])
        else:
            urows.append([s["label"], u["unplaced_pass"], "not projected", "-", "-", "-", "-", "-"])
    out += L.md_table(["call set", "unplaced PASS (linear_vcf)", "projected PASS rows", "easy", "difficult", "extreme",
                       "no position (iii)", "in HG008T BED"], urows) + [""]
    rs = [s for s in summ if "rgn_vs_bed" in s]
    if rs:
        out += [f"## {sec}. COLO829T: truth RGN annotation against the SMaHT BED assignment (truth records)", ""]
        out += L.md_table(["call set", "truth RGN -> SMaHT BED: n"], [[s["label"], ", ".join(f"{k}: {v}" for k, v in sorted(s["rgn_vs_bed"].items()))]
                                                       for s in [next(x for x in rs if x["kind"] == k) for k in ("SNV", "INDEL") if any(x["kind"] == k for x in rs)]]) + [""]
    (HERE / "tables.md").write_text("\n".join(out) + "\n")
    (HERE / "summary.json").write_text(json.dumps(dict(generated=time.strftime("%Y-%m-%d %H:%M"), call_sets=summ,
                                                       smaht_sizes=sizes), indent=1) + "\n")
    L.log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Scout prototype (read-only): Table R1.3a with the six-tumor models on three references. Prints to stdout only."""
import json, sys, time
from pathlib import Path
import numpy as np
T0 = time.time()
EXP = Path("/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments")
sys.path.insert(0, str(EXP / "E16_bootstrap")); sys.path.insert(0, str(EXP / "E10b_giab_strata"))
import callsets as C
import e16_bootstrap as E
import giab_strata as G
SMAHT_EXTREME = "/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_extreme_v2.bed"
GNAME = {"d9": "HPRC v1.1", "d46": "HPRC v2.1"}
GIAB_ST = [("alldifficultregions", "all difficult"), ("alllowmapandsegdupregions", "low mappability or segdup"),
           ("AllTandemRepeatsandHomopolymers_slop5", "tandem repeats and homopolymers"), ("notinalldifficultregions", "not difficult")]
HP7 = "AllHomopolymers_ge7bp_imperfectge11bp_slop5"
PARTS = ("tp-baseline", "fn", "tp", "fp")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs"); G6 = RUNS / "m4m6_graphs_20261110"; NEW = RUNS / "snv_newpipe_20261109"
SETS = [(t, p) for t in ("HG008T", "COLO829T", "HCC1395T") for p in ("Illumina", "PacBio", "ONT")]
def ds(t, p): return "COLO829T_fiberseq" if (t, p) == ("COLO829T", "PacBio") else f"{t}_{p}"
def lab(t, p): return f"{t} " + {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if t == "COLO829T" else "PacBio HiFi"}[p]
def vdir(kind, t, p, arm):
    r = "" if kind == "SNV" else ("pp_" if p == "ONT" else "ppr8_")
    if arm in ("linear", "d9"):
        return G6 / f"multicancer6_{p}_{kind}_{'linear_nositef' if arm == 'linear' else 'd9'}" / f"vcf_{r}chr1_{t}"
    if kind == "INDEL":
        return RUNS / f"release_s0_{p}_INDEL" / f"vcf_{r}chr1_{ds(t, p)}_d46"
    return NEW / f"multicancer6_{p}_SNV_newpipe" / f"vcf_chr1_{ds(t, p)}_d46"
def bed_chr1(path):
    iv = sorted((int(c[1]), int(c[2])) for c in (l.split("\t")[:3] for l in open(path)) if c[0] == "chr1")
    m = []
    for s, e in iv:
        if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
        else: m.append([s, e])
    a = np.array(m, dtype=np.int64).reshape(-1, 2); return a[:, 0], a[:, 1]
ivs = {"extreme": bed_chr1(SMAHT_EXTREME)}
ivs.update({name: G.chr1_intervals(G.GIAB / rel) for name, rel, _ in G.STRATA if name in dict(GIAB_ST) or name == HP7})
print("strata loaded %.0fs" % (time.time() - T0), file=sys.stderr)
bs = E.BlockSets(C.fasta_path(), 2000, 20261015)
print("blocksets %.0fs" % (time.time() - T0), file=sys.stderr)
def score(d, which):
    rep = json.loads((d / "eval_bed/report.json").read_text())
    w = d / "eval_bed/vcfeval" / which
    assert C._rtg_complete(w), w
    rec = {p: np.array([r[1] for r in C.read_records(w / f"{p}.vcf.gz") if r[0] == "chr1"], dtype=np.int64) for p in PARTS}
    a = rep["raw" if which == "raw_calls" else "pon"]["at_threshold"]
    assert (len(rec["tp-baseline"]), len(rec["fn"]), len(rec["fp"])) == (a["tp_baseline"], a["fn"], a["fp"]), (w, a)
    blocks = bs.get("bed", rep["bed"], G.WIDTH)
    out = {}
    for name, iv in ivs.items():
        m = E.block_counts({p: rec[p][G.inside(rec[p], iv)].tolist() for p in PARTS}, blocks)
        out[name] = (m.sum(axis=0), blocks.mult @ m, blocks.key)
    return out
def paired(a, b):
    assert a[2] == b[2]
    lo, hi, _, _ = E.ci(E.metrics(a[1])[2] - E.metrics(b[1])[2])
    return float(E.metrics(a[0])[2] - E.metrics(b[0])[2]), lo, hi
f3 = lambda x: "-" if np.isnan(x) else f"{x:.3f}"
ci = lambda t: f"{t[0]:+.3f} [{t[1]:+.3f}, {t[2]:+.3f}]"
res = {}
for kind in ("SNV", "INDEL"):
    for t, p in SETS:
        for arm in ("linear", "d9", "d46"):
            for w in ("raw_calls", "pon_calls"):
                res[(kind, t, p, arm, w)] = score(vdir(kind, t, p, arm), w)
    print(kind, "scored %.0fs" % (time.time() - T0), file=sys.stderr)
for w, title in (("raw_calls", "BEFORE population filtering"), ("pon_calls", "AFTER population filtering")):
    print(f"\n## {title}\n")
    print("| chr1 | graph | SMaHT extreme: PASS calls (precision), GRCh38-only -> HPRC | SMaHT extreme dF1 | " + " | ".join(f"GIAB {n} dF1" for _, n in GIAB_ST) + " | INDEL F1 HP>=7: GRCh38-only -> HPRC | truth SNV in extreme | truth INDEL in HP7 |")
    print("|---" * 10 + "|")
    for t, p in SETS:
        lin = res[("SNV", t, p, "linear", w)]
        for arm in ("d9", "d46"):
            a = res[("SNV", t, p, arm, w)]
            t0, t1 = lin["extreme"][0], a["extreme"][0]
            row = [f"{int(t0[2]+t0[3]):,} ({f3(E.metrics(t0)[0])}) -> {int(t1[2]+t1[3]):,} ({f3(E.metrics(t1)[0])})", ci(paired(a["extreme"], lin["extreme"]))]
            row += [ci(paired(a[st], lin[st])) for st, _ in GIAB_ST]
            i0 = res[("INDEL", t, p, "linear", w)][HP7][0]; i1 = res[("INDEL", t, p, arm, w)][HP7][0]
            row.append(f"{f3(float(E.metrics(i0)[2]))} -> {f3(float(E.metrics(i1)[2]))}")
            row.append(f"{int(t0[0]+t0[1])}/{int(t1[0]+t1[1])}"); row.append(f"{int(i0[0]+i0[1])}/{int(i1[0]+i1[1])}")
            print(f"| {lab(t, p)} | {GNAME[arm]} | " + " | ".join(row) + " |")
print("extra (records outside block set):", {k: b.extra for k, b in bs.cache.items()}, "nblocks", {k: len(b.ids) for k, b in bs.cache.items()})
print("total %.0fs" % (time.time() - T0), file=sys.stderr)

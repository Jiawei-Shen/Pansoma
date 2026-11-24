#!/usr/bin/env python3
"""E19 on the current models (user 2026-11-15: "is this analysis on the latest results?"; table layout of the SNV table
in the letter: Platform | Reference | FP | Precision | Recall | F1 | dF1, each 0 / 1 / >=2).

tables.md of this directory (2026-10-16) reads the 7fc7396 single-tumor runs of 2026-10-12 (old SNV and INDEL
recipes). This script repeats the W = 50 bp stratification with the same functions (germline_density.process: rtg
eval_bed pon_calls = PASS calls after the population filter in the benchmark BED; germline records of the normal
within 50 bp, the variant's own interval excluded; strata 0 / 1 / 2+) on the models of the current graph comparison
(user 2026-11-09 / 11-12); one seed (seed 0) unless stated. R = pansoma_net_v2_runs, G = R/m4m6_graphs_20261110, N = 6 (six-tumor model, every
tumor's chr1 held out) or 4 (four-cell-line model, HG008T and COLO829T unseen):

  INDEL  GRCh38-only  G/multicancer<N>_<P>_INDEL_linear_nositef   (no population site features as a model input)
         HPRC v1.1    G/multicancer<N>_<P>_INDEL_d9
         HPRC v2.1    R/release_s<k>_<P>_INDEL (N = 6), R/unseen_s<k>_<P>_INDEL (N = 4), seeds k = 0, 1, 2 (the
                      final batch). Two rows (user 2026-11-13: compare with the best v2.1 seed, not seed 0, and show
                      the mean): the best seed = highest whole-chr1 PASS F1 on that tumor and platform (as E33
                      indel_hcc1395_latest.py), and the mean of the three seeds per cell (a seed with no call in a
                      stratum has no precision / F1 there and is left out of that mean).
         PASS rule ppr8 for Illumina / PacBio HiFi / Fiber-seq, pp for ONT (user 2026-11-13).
  SNV    GRCh38-only  G/multicancer6_<P>_SNV_linear_nositef
         HPRC v1.1    G/multicancer6_<P>_SNV_d9
         HPRC v2.1    R/snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe
         six-tumor model only. A call set that does not exist yet is printed as "pending"; rerun to fill it.

The stratum of a truth variant and whether it has an exact candidate come from the tensor labels, not from the model.
The final-recipe single-tumor INDEL runs (ris_ppr8_20261026) cannot be used: for GRCh38-only and HPRC v1.1 only their
report.json is on this cluster (no call VCFs).
    /wanglab/jshen/anaconda3/bin/python strata_latest.py > strata_latest.md
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import germline_density as G  # noqa: E402

L = G.L
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
GR = R / "m4m6_graphs_20261110"
G.OUT = L.DATA / "E19_latest"
ARMS = (("linear", "GRCh38-only"), ("d9", "HPRC v1.1"), ("d46", "HPRC v2.1"))
TUMORS = ("HG008T", "COLO829T", "HCC1395T")
PLATS = ("Illumina", "PacBio", "ONT")
MODELS = {"INDEL": ((6, "six-tumor model"), (4, "four-cell-line model")), "SNV": ((6, "six-tumor model"),)}
SL = ("0", "1", "2+")


def pname(t, p):
    return "Fiber-seq" if (t, p) == ("COLO829T", "PacBio") else "PacBio HiFi" if p == "PacBio" else p


def vcf_dir(kind, n, t, p, arm, seed=0):
    """The call set's directory, or None while it does not exist."""
    if kind == "INDEL":
        rule = "pp" if p == "ONT" else "ppr8"
        if arm == "d46":
            vs = sorted((R / f"{'release' if n == 6 else 'unseen'}_s{seed}_{p}_INDEL").glob(f"vcf_{rule}_chr1_{t}_*_d46"))
        else:
            vs = [GR / f"multicancer{n}_{p}_INDEL_{'d9' if arm == 'd9' else 'linear_nositef'}" / f"vcf_{rule}_chr1_{t}"]
    elif arm == "d46":
        vs = sorted((R / "snv_newpipe_20261109" / f"multicancer{n}_{p}_SNV_newpipe").glob(f"vcf_chr1_{t}_*_d46"))
    else:
        vs = [GR / f"multicancer{n}_{p}_SNV_{'d9' if arm == 'd9' else 'linear_nositef'}" / f"vcf_chr1_{t}"]
    assert len(vs) <= 1, vs
    return vs[0] if vs and L.complete(vs[0]) else None


def call_set(kind, n, t, p, arm, seed=0):
    v = vcf_dir(kind, n, t, p, arm, seed)
    if v is None:
        return None
    lin = sorted(v.glob(f"*.{kind}.linear.vcf.gz"))
    assert len(lin) == 1, v
    st = lin[0].name[: -len(f".{kind}.linear.vcf.gz")]
    sample, platform, a, purity = L.set_parts(st[: -len("_snvaf05")] if st.endswith("_snvaf05") else st)
    assert sample == t and a == arm and purity == "P100", (st, t, arm)
    run = v.parent
    return L.CallSet(run=f"{run.parent.name}__{run.name}", vcf_dir=v, group="insample", variant="main", model=run.name,
                     model_sample=f"m{n}", model_platform=p, kind=kind, set=st, sample=t, platform=platform, arm=arm,
                     purity=purity, prefix=f"{st}.{kind}")


def stats(rows, stratum):
    """Counts of one stratum ('all' = every record): truth, truth with a tensor, TP-baseline, TP calls, FP."""
    c = dict(truth=0, tensor=0, tpb=0, tpc=0, fp=0)
    for r in rows:
        if stratum != "all" and G.stratum(int(r["n50"])) != stratum:
            continue
        if r["part"] in ("tp-baseline", "fn"):
            c["truth"] += 1
            c["tensor"] += r["has_tensor"] == "1"
            c["tpb"] += r["part"] == "tp-baseline"
        else:
            c["tpc" if r["part"] == "tp" else "fp"] += 1
    c["p"], c["r"], c["f"] = L.prf(c["tpb"], c["truth"] - c["tpb"], c["tpc"], c["fp"])
    return c


def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def f3(x, sign=False):
    return "–" if x is None else f"{x:+.3f}" if sign else f"{x:.3f}"


def tri(vals, sign=False):
    return " / ".join(f3(v, sign) if not isinstance(v, int) else str(v) for v in vals)


def main():
    res, pending, seedinfo = {}, [], {}
    for kind, models in MODELS.items():
        for n, _ in models:
            for t in TUMORS:
                for p in PLATS:
                    for arm, _ in ARMS:
                        seeds = (0, 1, 2) if (kind, arm) == ("INDEL", "d46") else (0,)
                        per = []
                        for seed in seeds:
                            cs = call_set(kind, n, t, p, arm, seed)
                            if cs is None:
                                break
                            rows = G.process(cs)
                            rep = cs.report("eval_bed")["pon"]["at_threshold"]
                            a = stats(rows, "all")
                            assert (a["tpb"], a["truth"] - a["tpb"], a["tpc"], a["fp"]) == (rep["tp_baseline"], rep["fn"], rep["tp_call"], rep["fp"]), cs.id
                            per.append(({s: stats(rows, s) for s in SL + ("all",)}, rows, seed))
                        if len(per) < len(seeds):
                            pending.append((kind, n, t, p, arm))
                            continue
                        best = max(per, key=lambda x: x[0]["all"]["f"] or 0.0)
                        res[kind, n, t, p, arm] = best[0], best[1]
                        if len(per) > 1:
                            seedinfo[kind, n, t, p, arm] = best[2], per
    out = ["# E19 on the current models (strata_latest.py)", "",
           "chr1, benchmark regions, PASS calls after the population filter. Stratum = germline variants of the normal "
           "within 50 bp of the variant (0 / 1 / ≥2); recall by the truth variant's stratum, precision by the call's "
           "stratum. Models, PASS rules and directories: the docstring of strata_latest.py. One seed (seed 0), except "
           "INDEL on HPRC v2.1: three seeds, shown as the best seed (highest whole-chr1 PASS F1 on that tumor and "
           "platform) and as the mean of the three. ΔF1 = HPRC − GRCh38-only of the same platform. pending = the call "
           "set does not exist yet.", ""]
    sec = 0
    for kind, models in MODELS.items():
        for n, mname in models:
            for t in TUMORS:
                sec += 1
                truth = next((res[k][0] for k in res if k[:3] == (kind, n, t)), None)
                cnt = "truth " + kind + "s 0 / 1 / ≥2: " + tri([truth[s]["truth"] for s in SL]) if truth else "no call set yet"
                out += [f"## {sec}. chr1 {kind}, {t}, {mname} ({cnt})", "",
                        "| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |",
                        "|---|---|---|---|---|---|---|"]
                for p in PLATS:
                    base = res.get((kind, n, t, p, "linear"))
                    for i, (arm, gname) in enumerate(ARMS):
                        first = pname(t, p) if i == 0 else ""
                        if (kind, n, t, p, arm) not in res:
                            out.append(f"| {first} | {gname} | pending | pending | pending | pending | {'–' if arm == 'linear' else 'pending'} |")
                            continue
                        S = res[kind, n, t, p, arm][0]
                        assert truth is None or all(S[s]["truth"] == truth[s]["truth"] for s in SL), (kind, n, t, p, arm)
                        if arm == "linear":
                            d = "–"
                        elif base is None:
                            d = "pending"
                        else:
                            d = tri([None if S[s]["f"] is None or base[0][s]["f"] is None else S[s]["f"] - base[0][s]["f"] for s in SL], True)
                        si = seedinfo.get((kind, n, t, p, arm))
                        out.append(f"| {first} | {gname}{f', best seed (seed {si[0]})' if si else ''} | {tri([S[s]['fp'] for s in SL])} | {tri([S[s]['p'] for s in SL])} | "
                                   f"{tri([S[s]['r'] for s in SL])} | {tri([S[s]['f'] for s in SL])} | {d} |")
                        if si:
                            m = {s: {k: mean([x[0][s][k] for x in si[1]]) for k in ("fp", "p", "r", "f")} for s in SL}
                            dm = "pending" if base is None else tri([None if m[s]["f"] is None or base[0][s]["f"] is None else m[s]["f"] - base[0][s]["f"] for s in SL], True)
                            out.append(f"|  | {gname}, mean of 3 seeds | " + " / ".join(f"{m[s]['fp']:.1f}" for s in SL) + f" | {tri([m[s]['p'] for s in SL])} | "
                                       f"{tri([m[s]['r'] for s in SL])} | {tri([m[s]['f'] for s in SL])} | {dm} |")
                out.append("")
    sec += 1
    out += [f"## {sec}. Truth variants with an exact candidate per stratum, and the whole-chr1 result", "",
            "| chr1 call set | Exact candidate (0 / 1 / ≥2) | chr1 all: P / R / F1 |", "|---|---|---|"]
    for (kind, n, t, p, arm), (S, _) in res.items():
        si = seedinfo.get((kind, n, t, p, arm))
        for seed, X in ([(x[2], x[0]) for x in si[1]] if si else [(None, S)]):
            tag = "" if seed is None else f", seed {seed}" + (" (best)" if seed == si[0] else "")
            out.append(f"| {kind}, {dict(MODELS[kind])[n]}, {t} {pname(t, p)}, {dict(ARMS)[arm]}{tag} | "
                       + " / ".join(f"{100 * X[s]['tensor'] / X[s]['truth']:.1f}%" for s in SL)
                       + f" | {tri([X['all'][k] for k in ('p', 'r', 'f')])} |")
    out.append("")
    sec += 1
    out += [f"## {sec}. Paired recall in stratum ≥2, HPRC vs GRCh38-only (same truth variants)", "",
            "| chr1 comparison | truth ≥2 | found by HPRC only / by GRCh38-only only | exact McNemar p |", "|---|---|---|---|"]
    for kind, models in MODELS.items():
        for n, mname in models:
            for t in TUMORS:
                for p in PLATS:
                    for arm, gname in ARMS[1:]:
                        if (kind, n, t, p, arm) not in res or (kind, n, t, p, "linear") not in res:
                            continue
                        ta, tb = G.truth_status(res[kind, n, t, p, arm][1], 50), G.truth_status(res[kind, n, t, p, "linear"][1], 50)
                        assert set(ta) == set(tb)
                        keys = [k for k in ta if ta[k][0] == "2+"]
                        x, y = sum(ta[k][1] and not tb[k][1] for k in keys), sum(tb[k][1] and not ta[k][1] for k in keys)
                        si = seedinfo.get((kind, n, t, p, arm))
                        out.append(f"| {kind}, {mname}, {t} {pname(t, p)}, {gname}{f' (best seed, seed {si[0]})' if si else ''} | {len(keys)} | {x} / {y} | {G.pfmt(L.binom_two_sided(x, y))} |")
    out += ["", "## Pending call sets", ""] + ([f"- {k} m{n} {t} {pname(t, p)} {dict(ARMS)[a]}" for k, n, t, p, a in pending] or ["none"])
    print("\n".join(out))


if __name__ == "__main__":
    main()

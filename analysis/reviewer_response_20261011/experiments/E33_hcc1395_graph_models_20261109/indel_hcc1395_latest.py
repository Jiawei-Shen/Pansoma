"""HCC1395T chr1 INDEL results on the latest models (2026-11-13): GRCh38-only vs HPRC v1.1 vs v2.1, four lines (m4)
vs six tumors (m6), against DeepSomatic-TO and ClairS-TO. chr1 only (user 2026-11-12: no HCC1395T chr20 rows).

Pansoma (all the final pre-registered INDEL recipe; PASS ppr8 for Illumina / PacBio HiFi and pp, the post-PoN validation
best-F1 threshold, for ONT (user 2026-11-13: ONT INDEL calls at pp, four lines and six tumors); one model per platform /
graph / tumor set):
- v2.1 (d46): unseen_s{0,1,2}_<P>_INDEL (m4) / release_s{0,1,2}_<P>_INDEL (m6), vcf_ppr8_chr1_HCC1395T_<P>_d46;
  one row per seed and the seed mean (the final batch's pre-registered summary).
- v1.1 (d9) and GRCh38-only (linear, no site features as input, user 2026-11-12): m4m6_graphs_20261110/
  multicancer{4,6}_<P>_INDEL_{d9,linear_nositef}, vcf_ppr8_chr1_HCC1395T. Also multicancer4_PacBio_INDEL_linear
  (GRCh38-only WITH site features, the arm before the nositef decision; its m6 twin was stopped at epoch 8).
Comparators: E8 score, eval_bed: pass_pon (PASS after our filter), pass (PASS, no filter), all_pon (curve).

Columns (rtg vcfeval, SEQC2 BED, truth 133 INDELs on chr1):
  PASS after filter P / R / F1 and TP / FP / FN (pon.at_threshold); PASS before filter (raw.at_threshold) P / R / F1;
  best F1 on the post-filter curve (pon.best: threshold-free, every record ranked); curve end = highest recall of the
  post-filter curve (pon.ceiling).
Tensor level (test_chr1 metrics.json; truth 124 = position in the BED): truth alleles with an exact tensor
  (candidate recall of the graph's tensors), somatic AP over tensors.
Writes indel_hcc1395_latest.md next to this script.
usage: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python indel_hcc1395_latest.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
G = R / "m4m6_graphs_20261110"
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
PLATS = ("Illumina", "PacBio", "ONT")
GNAME = {"linear_nositef": "GRCh38-only", "linear": "GRCh38-only + site features", "d9": "HPRC v1.1", "d46": "HPRC v2.1"}


def rule(p):
    """The PASS rule's VCF directory prefix: pp for ONT (user 2026-11-13), ppr8 for the other platforms."""
    return "pp" if p == "ONT" else "ppr8"


def runs(p, m, g, seed=0):
    """(run dir, vcf dir, tensor metrics) or None if not finished."""
    if g == "d46":
        d = R / f"{'unseen' if m == 4 else 'release'}_s{seed}_{p}_INDEL"
        v = d / f"vcf_{rule(p)}_chr1_HCC1395T_{p}_d46"
        t = d / f"test_chr1_HCC1395T_{p}_d46" / f"HCC1395T_{p}_d46.tensors.INDEL.metrics.json"
    else:
        d = G / f"multicancer{m}_{p}_INDEL_{g}"
        v = d / f"vcf_{rule(p)}_chr1_HCC1395T"
        t = d / "test_chr1" / f"HCC1395T_{p}_{g.split('_')[0]}.tensors.INDEL.metrics.json"
    rep = v / "eval_bed" / "report.json"
    return (d, rep, t) if rep.exists() else None


def pansoma(p, m, g, seed=0):
    x = runs(p, m, g, seed)
    if x is None:
        return None
    d, rep, t = x
    r = json.loads(rep.read_text())
    assert r["chromosomes"] == ["chr1"] and "HCC1395T" in r["calls"], rep
    a, b = r["pon"]["at_threshold"], r["raw"]["at_threshold"]
    tm = json.loads(t.read_text()) if t.exists() else None
    best_ep = int(next(l for l in (d / "train.log").read_text().splitlines() if l.startswith("done; best epoch")).split()[3])
    return dict(P=a["precision"], R=a["recall"], F1=a["f1"], TP=a["tp_baseline"], FP=a["fp"], FN=a["fn"],
                rP=b["precision"], rR=b["recall"], rF1=b["f1"], rFP=b["fp"],
                best=r["pon"]["best"]["f1"], bestP=r["pon"]["best"]["precision"], bestR=r["pon"]["best"]["recall"],
                end=r["pon"]["ceiling"], n=r["pon"]["baseline"], epoch=best_ep,
                cand=(tm["truth"]["with_tensor"], tm["truth"]["truth_alleles"]) if tm else None,
                ap=tm["somatic_ap"] if tm else None)


def comparator(c, p):
    e = json.loads((E8 / c / f"HCC1395T_{p}" / "INDEL" / "report.json").read_text())["eval_bed"]
    a, b, allp = e["pass_pon"]["at_threshold"], e["pass"]["at_threshold"], e["all_pon"]
    return dict(P=a["precision"], R=a["recall"], F1=a["f1"], TP=a["tp_baseline"], FP=a["fp"], FN=a["fn"],
                rP=b["precision"], rR=b["recall"], rF1=b["f1"], rFP=b["fp"],
                best=allp["best"]["f1"], bestP=allp["best"]["precision"], bestR=allp["best"]["recall"],
                end=allp["ceiling"], n=e["pass_pon"]["baseline"], epoch=None, cand=None, ap=None)


def f3(x):
    return "–" if x is None else f"{x:.3f}"


def row(name, x, seeds=None):
    if x is None:
        return f"| {name} | still training | | | | | | | |"
    sr = f" (s0–s2 {min(seeds):.3f}–{max(seeds):.3f})" if seeds else ""
    cand = f"{x['cand'][0]} / {x['cand'][1]}" if x["cand"] else "–"
    return (f"| {name} | {x['P']:.3f} / {x['R']:.3f} / **{x['F1']:.3f}**{sr} | {x['TP']} / {x['FP']} / {x['FN']} | "
            f"{x['rP']:.3f} / {x['rR']:.3f} / {x['rF1']:.3f} ({x['rFP']} FP) | {x['best']:.3f} ({x['bestP']:.3f} / {x['bestR']:.3f}) | "
            f"{x['end']:.3f} | {cand} | {f3(x['ap'])} | {x['epoch'] or '–'} |")


def main():
    L = ["# HCC1395T chr1 INDEL, latest models (generated by indel_hcc1395_latest.py)", "",
         "rtg vcfeval in the SEQC2 BED, chr1 truth 133 INDELs. PASS = ppr8 for Illumina / PacBio HiFi and pp for ONT "
         "(Pansoma) / the tool's PASS (comparators). "
         "Tensor columns: test_chr1 metrics (truth 124 by position in the BED).", ""]
    for p in PLATS:
        L += [f"## {p}", "",
              "| model | PASS after filter P / R / **F1** | TP / FP / FN | PASS before filter P / R / F1 | best F1 on curve (P / R) | "
              "curve end R | truth with exact tensor | somatic AP (tensors) | epoch |",
              "|---|---|---|---|---|---|---|---|---|"]
        for m in (4, 6):
            for g in ("linear_nositef", "linear", "d9", "d46"):
                if g == "linear" and not (m == 4 and p == "PacBio"):
                    continue
                name = f"Pansoma {'four lines' if m == 4 else 'six tumors'}, {GNAME[g]}"
                if g != "d46":
                    L.append(row(name, pansoma(p, m, g)))
                    continue
                # v2.1: every seed (user 2026-11-13: not seed 0 alone), then the pre-registered seed mean
                xs = [pansoma(p, m, g, s) for s in (0, 1, 2)]
                for s, x in enumerate(xs):
                    L.append(row(f"{name}, seed {s}", x))
                f1 = [x["F1"] for x in xs]
                k = max(range(3), key=lambda i: f1[i])  # user 2026-11-13: compare with the best seed
                L.append(f"| {name}, best seed (seed {k}) / mean of 3 | F1 **{f1[k]:.3f}** / {sum(f1) / 3:.3f} | | | | | | | |")
        for c, cn in (("deepsomatic", "DeepSomatic-TO"), ("clairs_to", "ClairS-TO")):
            L.append(row(cn, comparator(c, p)))
        L.append("")
    (HERE / "indel_hcc1395_latest.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

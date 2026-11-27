#!/usr/bin/env python3
"""Latest results of every current model in one file (user 2026-11-23). chr1 of the three test tumors (HG008T, COLO829T,
HCC1395T), benchmark regions (eval_bed), rtg vcfeval, PASS calls after the population filter unless a column says
"before": P / R / F1 from <run>/<vcf dir>/eval_bed/report.json ("pon" / "raw" -> "at_threshold").

Models (R = /scratch/jshen/data/pansoma_net_v2_runs), seed 0 unless a seed is named:
  six-tumor (m6) and four-line (m4) x GRCh38-only / HPRC v1.1 / HPRC v2.1 x Illumina / PacBio HiFi / ONT x SNV / INDEL
    GRCh38-only  m4m6_graphs_20261110/multicancer{4,6}_<P>_<KIND>_linear_nositef  (no population site features as input)
    HPRC v1.1    m4m6_graphs_20261110/multicancer{4,6}_<P>_<KIND>_d9
    HPRC v2.1    SNV m6 snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe; SNV m4 Illumina / PacBio
                 multicancer4_<P>_SNV_sitef_af05_lnf25w1 (rule-2 VCFs *_r095; COLO829T on the AF 0.05 set), m4 ONT
                 snv_newpipe_20261109/multicancer4_ONT_SNV_newpipe; INDEL m6 release_s{0,1,2}_<P>_INDEL, m4
                 unseen_s{0,1,2}_<P>_INDEL
    PASS: SNV rule 2 (post-filter validation recall 0.95 of the callable positives); INDEL ppr8 (Illumina / PacBio HiFi),
    pp (ONT).
  six-tumor SNV without site features as a model input on v1.1 / v2.1: trained on RIS, runs not on this cluster; the
    numbers are read from analysis/m6_nositef_graphs_20261118/results.md.
  COLO829T single-tumor SNV: snv_newpipe_20261109/COLO829_<P>_SNV_<graph>_newpipe[_s1|_s2], the path-count ablation
    (_nopc_) and the GRCh38-only models without site features (_linear_newpipe_nositef).
  Comparators: E8 score (DeepSomatic-TO, ClairS-TO), the tool's PASS after the same filter.
The truth count (TP + FN) of every cell of a row is asserted equal.
usage: python3 collect.py > results.md"""
import json
import re
import statistics
from pathlib import Path

R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
G = R / "m4m6_graphs_20261110"
N = R / "snv_newpipe_20261109"
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
REPO = Path(__file__).resolve().parents[2]
TUMORS = ["HCC1395T", "COLO829T", "HG008T"]
PLATS = ["Illumina", "PacBio", "ONT"]


def ds(tumor, p):
    return "COLO829T_fiberseq" if (tumor, p) == ("COLO829T", "PacBio") else f"{tumor}_{p}"


def pname(tumor, p):
    return {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if tumor == "COLO829T" else "PacBio HiFi"}[p]


def rep(d, stage="pon"):
    f = Path(d) / "eval_bed" / "report.json"
    return json.loads(f.read_text())[stage]["at_threshold"] if f.exists() else None


def full(d):
    f = Path(d) / "eval_bed" / "report.json"
    return json.loads(f.read_text()) if f.exists() else None


def rule(kind, p):
    return "" if kind == "SNV" else ("pp_" if p == "ONT" else "ppr8_")


def run_dir(n, kind, p, graph, seed=0):
    """(run directory, vcf directory name pattern) of one multi-tumor model."""
    if graph in ("linear_nositef", "d9"):
        return G / f"multicancer{n}_{p}_{kind}_{graph}", "vcf_{r}chr1_{tumor}"
    if kind == "INDEL":
        return R / f"{'release' if n == 6 else 'unseen'}_s{seed}_{p}_INDEL", "vcf_{r}chr1_{ds}_d46"
    if n == 6 or p == "ONT":
        return N / f"multicancer{n}_{p}_SNV_newpipe", "vcf_chr1_{ds}_d46"
    return R / f"multicancer4_{p}_SNV_sitef_af05_lnf25w1", "vcf_chr1_{ds}_d46{af}_r095"


def cell_dir(n, kind, p, graph, tumor, seed=0):
    run, pat = run_dir(n, kind, p, graph, seed)
    af = "_snvaf05" if tumor == "COLO829T" and p != "ONT" else ""
    return run / pat.format(r=rule(kind, p), tumor=tumor, ds=ds(tumor, p), af=af)


def prf(a, bold=False):
    if a is None:
        return "-"
    f1 = f"**{a['f1']:.3f}**" if bold else f"{a['f1']:.3f}"
    return f"{a['precision']:.3f} / {a['recall']:.3f} / {f1}"


def truth(a):
    return a["tp_baseline"] + a["fn"]


def comparator(caller, kind, tumor, p):
    e = json.loads((E8 / caller / ds(tumor, p) / kind / "report.json").read_text())
    return e["eval_bed"]["pass_pon"]["at_threshold"]


def best_epoch(run):
    log = Path(run) / "train.log"
    if not log.exists():
        return "-"
    m = re.findall(r"^done; best epoch (\d+) \(postpon ([0-9.]+)", log.read_text(), re.M)
    return f"{m[-1][0]} ({m[-1][1]})" if m else "-"


def nositef_ris():
    """F1 / P / R of the RIS models from analysis/m6_nositef_graphs_20261118/results.md, keyed (tumor, platform, graph)."""
    out = {}
    for line in (REPO / "analysis/m6_nositef_graphs_20261118/results.md").read_text().splitlines():
        c = [x.strip() for x in line.strip("|").split("|")]
        if len(c) == 14 and c[2] in ("HPRC v1.1", "HPRC v2.1") and c[0] in TUMORS:
            out[(c[0], c[1], c[2])] = dict(precision=float(c[9]), recall=float(c[10]), f1=float(c[11]),
                                           tp_baseline=int(c[6]), fn=int(c[8]), fp=int(c[7]))
    return out


RIS = nositef_ris()
print("# Latest results, 2026-11-23 (generated by `collect.py`; its docstring names every run)\n")
print("Every number: chr1 of the test tumor, benchmark regions, rtg vcfeval, PASS calls, P / R / F1. \"after\" = after the "
      "population filter. Seed 0 unless a seed is named. HCC1395T chr1 is held out of every multi-tumor model; HG008T and "
      "COLO829T are unseen tumors for the four-line models and have chr1 held out of the six-tumor models.")

# ---------------------------------------------------------------- 1. six-tumor v2.1 vs comparators
print("\n## 1. Six-tumor models on HPRC v2.1 against DeepSomatic-TO and ClairS-TO, chr1 after the population filter\n")
for kind in ("SNV", "INDEL"):
    head = "Pansoma six-tumor, chr1" if kind == "SNV" else "Pansoma six-tumor seed 0, chr1"
    extra = "" if kind == "SNV" else " chr1 F1 seeds 0 / 1 / 2 | chr1 F1 mean ± SD |"
    print(f"\n### 1{'a' if kind == 'SNV' else 'b'}. {kind} (PASS: {'rule 2' if kind == 'SNV' else 'ppr8; ONT pp'})\n")
    print(f"| tumor | platform | chr1 truth | {head} | DeepSomatic-TO, chr1 | ClairS-TO, chr1 |{extra}")
    print("|---|---|---|---|---|---|" + ("---|---|" if kind == "INDEL" else ""))
    wins = 0
    for tumor in TUMORS:
        for p in PLATS:
            a = rep(cell_dir(6, kind, p, "d46", tumor))
            d, c = comparator("deepsomatic", kind, tumor, p), comparator("clairs_to", kind, tumor, p)
            assert truth(a) == truth(d) == truth(c), (kind, tumor, p)
            best = max(a["f1"], d["f1"], c["f1"])
            wins += a["f1"] == best
            row = f"| {tumor} | {pname(tumor, p)} | {truth(a):,} | {prf(a, a['f1'] == best)} | {prf(d, d['f1'] == best)} | {prf(c, c['f1'] == best)} |"
            if kind == "INDEL":
                f1 = [rep(cell_dir(6, kind, p, "d46", tumor, s))["f1"] for s in range(3)]
                row += f" {' / '.join(f'{x:.3f}' for x in f1)} | {statistics.mean(f1):.3f} ± {statistics.stdev(f1):.3f} |"
            print(row)
    print(f"\nPansoma has the highest chr1 F1 in {wins} of 9 cells (seed 0).")

# ---------------------------------------------------------------- 2 / 3. graph comparison
GRAPHS = [("linear_nositef", "GRCh38-only"), ("d9", "HPRC v1.1"), ("d46", "HPRC v2.1")]
for n, sec, name in ((6, 2, "Six-tumor"), (4, 3, "Four-line")):
    print(f"\n## {sec}. {name} models on the three references, chr1 after the population filter\n")
    print("GRCh38-only models get no population site features as a model input (user 2026-11-12); the v1.1 / v2.1 models do.")
    for kind in ("SNV", "INDEL"):
        ris = n == 6 and kind == "SNV"
        cols = "| tumor | platform | chr1 truth | GRCh38-only, chr1 | HPRC v1.1, chr1 | HPRC v2.1, chr1 |"
        if ris:
            cols += " v1.1 no site features, chr1 | v2.1 no site features, chr1 |"
        if kind == "INDEL":
            cols += " v2.1 chr1 F1 seeds 0 / 1 / 2 (mean) |"
        print(f"\n### {sec}{'a' if kind == 'SNV' else 'b'}. {kind}\n\n{cols}")
        print("|" + "---|" * (cols.count("|") - 1))
        tally = {g: 0 for g, _ in GRAPHS}
        for tumor in TUMORS:
            for p in PLATS:
                cells = [rep(cell_dir(n, kind, p, g, tumor)) for g, _ in GRAPHS]
                assert len({truth(a) for a in cells}) == 1, (n, kind, tumor, p, [truth(a) for a in cells])
                best = max(a["f1"] for a in cells)
                for (g, _), a in zip(GRAPHS, cells):
                    tally[g] += a["f1"] == best
                row = f"| {tumor} | {pname(tumor, p)} | {truth(cells[0]):,} | " + " | ".join(prf(a, a["f1"] == best) for a in cells) + " |"
                if ris:
                    for g in ("HPRC v1.1", "HPRC v2.1"):
                        a = RIS[(tumor, pname(tumor, p), g)]
                        assert truth(a) == truth(cells[0]), (tumor, p, g)
                        row += f" {prf(a)} |"
                if kind == "INDEL":
                    f1 = [rep(cell_dir(n, kind, p, "d46", tumor, s))["f1"] for s in range(3)]
                    row += f" {' / '.join(f'{x:.3f}' for x in f1)} ({statistics.mean(f1):.3f}) |"
                print(row)
        print("\nHighest chr1 F1 of the three references (seed 0, bold): " + ", ".join(f"{nm} {tally[g]} of 9" for g, nm in GRAPHS) + ".")
    if n == 4:
        o = {t: rep(R / "multicancer4_ONT_SNV" / f"vcf_chr1_{t}_ONT_d46_normfix") for t in TUMORS}
        if all(o.values()):
            print("\nFour-line ONT SNV on v2.1: the cell above is the current-recipe model (multicancer4_ONT_SNV_newpipe), which "
                  "fails on chr1 (false positives in paralog-like regions). The original-recipe model (multicancer4_ONT_SNV, "
                  "normfix), still used in Supplementary Fig. R1.6 / Table R1.7b, has chr1 P / R / F1 after the filter: "
                  + "; ".join(f"{t} {prf(o[t])}" for t in TUMORS) + ".")

# ---------------------------------------------------------------- 4. before the filter
print("\n## 4. Six-tumor models, chr1 PASS before the population filter (F1; FP count in brackets)\n")
print("| kind | tumor | platform | GRCh38-only, chr1 before | HPRC v1.1, chr1 before | HPRC v2.1, chr1 before |")
print("|---|---|---|---|---|---|")
for kind in ("SNV", "INDEL"):
    for tumor in TUMORS:
        for p in PLATS:
            cells = [rep(cell_dir(6, kind, p, g, tumor), "raw") for g, _ in GRAPHS]
            print(f"| {kind} | {tumor} | {pname(tumor, p)} | " + " | ".join(f"{a['f1']:.3f} ({a['fp']:,} FP)" for a in cells) + " |")

# ---------------------------------------------------------------- 5. training summary
print("\n## 5. Selected epoch (post-filter validation F1 in brackets) and PASS threshold of each multi-tumor model\n")
print("| model | kind | platform | GRCh38-only: epoch, PASS thr | HPRC v1.1: epoch, PASS thr | HPRC v2.1: epoch, PASS thr |")
print("|---|---|---|---|---|---|")
for n, name in ((6, "six-tumor"), (4, "four-line")):
    for kind in ("SNV", "INDEL"):
        for p in PLATS:
            cs = []
            for g, _ in GRAPHS:
                f = full(cell_dir(n, kind, p, g, "HCC1395T"))
                cs.append(f"{best_epoch(run_dir(n, kind, p, g)[0])}, {f['threshold']:.3f}")
            print(f"| {name} | {kind} | {pname('HCC1395T', p)} | " + " | ".join(cs) + " |")

# ---------------------------------------------------------------- 6. COLO829T single-tumor
print("\n## 6. COLO829T single-tumor SNV models of the current pipeline, COLO829T chr1 after the population filter\n")
print("Illumina / Fiber-seq on the AF 0.05 sets, ONT at AF 0.08. \"no path count\" = the path-count channel ablation; "
      "\"GRCh38-only, no site features\" = the arm without the population site features as a model input.\n")
print("| platform | model | COLO829T chr1 seed 0 | seed 1 | seed 2 | chr1 F1 mean |")
print("|---|---|---|---|---|---|")
VARIANTS = [("HPRC v2.1", "d46_newpipe"), ("HPRC v1.1", "d9_newpipe"), ("GRCh38-only (with site features)", "linear_newpipe"),
            ("GRCh38-only, no site features", "linear_newpipe_nositef"), ("HPRC v2.1, no path count", "nopc_d46_newpipe"),
            ("HPRC v1.1, no path count", "nopc_d9_newpipe")]
for p in ("Illumina", "fiberseq", "ONT"):
    for label, v in VARIANTS:
        cells = []
        for s in range(3):
            run = N / f"COLO829_{p}_SNV_{v}{'' if s == 0 else f'_s{s}'}"
            ds_ = sorted(run.glob(f"vcf_chr1_COLO829T_{p}*")) if run.exists() else []
            cells.append(rep(ds_[0]) if ds_ else None)
        if not any(cells):
            continue
        f1 = [c["f1"] for c in cells if c]
        mean = f"{statistics.mean(f1):.3f}" if len(f1) > 1 else "-"
        print(f"| {'Fiber-seq' if p == 'fiberseq' else p} | {label} | " + " | ".join(prf(c) for c in cells) + f" | {mean} |")
    for c, cn in (("deepsomatic", "DeepSomatic-TO"), ("clairs_to", "ClairS-TO")):
        a = json.loads((E8 / c / f"COLO829T_{p}" / "SNV" / "report.json").read_text())["eval_bed"]["pass_pon"]["at_threshold"]
        print(f"| {'Fiber-seq' if p == 'fiberseq' else p} | {cn} | {prf(a)} | | | |")

# ---------------------------------------------------------------- 7. six-tumor SNV seeds, 8. path-count (RIS, pulled 2026-11-26)
def ris_set(tumor, p):
    """The test set name of the RIS six-tumor v2.1 SNV runs (vcf_chr1_<SET>)."""
    return f"{ds(tumor, p)}_d46" + ("_snvaf05" if p != "ONT" and tumor != "HG008T" else "")


SEED0 = R / "snv_newpipe_20261109"
print("\n## 7. Six-tumor SNV models on HPRC v2.1, three seeds, chr1 after the population filter\n")
print("Seed 0 = the model of section 1a; seeds 1 and 2 = `m6_seeds_20261123/multicancer6_<P>_SNV_newpipe_s<SEED>` (trained and "
      "tested on RIS, only `--seed` changed).\n")
print("| tumor | platform | chr1 seed 0 | chr1 seed 1 | chr1 seed 2 | chr1 F1 mean | chr1 F1 SD |")
print("|---|---|---|---|---|---|---|")
for tumor in ("HCC1395T", "COLO829T", "HG008T"):
    for p in ("Illumina", "PacBio", "ONT"):
        cells = [rep(SEED0 / f"multicancer6_{p}_SNV_newpipe" / f"vcf_chr1_{ds(tumor, p)}_d46")]
        cells += [rep(R / "m6_seeds_20261123" / f"multicancer6_{p}_SNV_newpipe_s{s}" / f"vcf_chr1_{ris_set(tumor, p)}") for s in (1, 2)]
        assert len({truth(c) for c in cells if c}) <= 1, (tumor, p)
        f1 = [c["f1"] for c in cells if c]
        ms = (f"{statistics.mean(f1):.3f}", f"{statistics.stdev(f1):.3f}") if len(f1) > 1 else ("-", "-")
        print(f"| {tumor} | {pname(tumor, p)} | " + " | ".join(prf(c) for c in cells) + f" | {ms[0]} | {ms[1]} |")

print("\n## 8. Six-tumor SNV models on HPRC v2.1 with and without the path-count channel, chr1 after the population filter\n")
print("With = the model of section 1a (seed 0); without = `m6_nopathcount_20261124/multicancer6_<P>_SNV_newpipe_nopc` "
      "(RIS, seed 0, `--drop-planes path_count`, otherwise the same recipe).\n")
print("| tumor | platform | chr1 truth | with path-count, chr1 | chr1 TP / FP | without path-count, chr1 | chr1 TP / FP | chr1 F1 change |")
print("|---|---|---|---|---|---|---|---|")
for tumor in ("HCC1395T", "COLO829T", "HG008T"):
    for p in ("Illumina", "PacBio", "ONT"):
        w = rep(SEED0 / f"multicancer6_{p}_SNV_newpipe" / f"vcf_chr1_{ds(tumor, p)}_d46")
        n = rep(R / "m6_nopathcount_20261124" / f"multicancer6_{p}_SNV_newpipe_nopc" / f"vcf_chr1_{ris_set(tumor, p)}")
        assert w is None or n is None or truth(w) == truth(n), (tumor, p)
        tf = lambda a: "-" if a is None else f"{a['tp_baseline']:,} / {a['fp']:,}"
        d = "-" if None in (w, n) else f"{n['f1'] - w['f1']:+.3f}"
        print(f"| {tumor} | {pname(tumor, p)} | {truth(w) if w else '-'} | {prf(w)} | {tf(w)} | {prf(n)} | {tf(n)} | {d} |")

#!/usr/bin/env python3
"""Letter ver4: the tables that can be filled from existing evaluation reports of the current six-tumor models
(no new computation). Writes tables.json (cell strings, in the letter's layouts) and tables.md next to this script.

chr1 of the test tumor, benchmark regions, rtg vcfeval, PASS calls: <run>/<vcf dir>/eval_bed/report.json,
"raw" = before the population filter, "pon" = after it, at_threshold.
  SNV   HPRC v2.1    R/snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe/vcf_chr1_<SET>_d46          (PASS = rule 2)
        HPRC v1.1    R/m4m6_graphs_20261110/multicancer6_<P>_SNV_d9/vcf_chr1_<TUMOR>
        GRCh38-only  R/m4m6_graphs_20261110/multicancer6_<P>_SNV_linear_nositef/vcf_chr1_<TUMOR>      (no site features)
        v1.1 / v2.1 without site features: trained on RIS; after-filter values parsed from
                     analysis/m6_nositef_graphs_20261118/results.md (their run directories are not on this cluster)
  INDEL HPRC v2.1    R/release_s{0,1,2}_<P>_INDEL/vcf_{ppr8|pp}_chr1_<SET>_d46   (six-tumor; pp on ONT)
                     R/unseen_s{0,1,2}_<P>_INDEL/...                              (four-tumor)
        v1.1 / GRCh38-only  R/m4m6_graphs_20261110/multicancer6_<P>_INDEL_{d9,linear_nositef}/vcf_{ppr8|pp}_chr1_<TUMOR>
  SNV single-tumor seeds  R/snv_newpipe_20261109/COLO829_<P>_SNV_{d46,d9}_newpipe[_s1|_s2]
Other inputs: E19_germline_density/strata_latest.md (sections 1 and 7), graph_absorbed_somatic_hcc1395t_20261112/
tables_v21_patient_frame.md. Truth counts (TP + FN) are asserted equal across the arms of every row.
usage: python3 make_tables.py"""
import json
import statistics
from pathlib import Path

R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
G = R / "m4m6_graphs_20261110"
N = R / "snv_newpipe_20261109"
HERE = Path(__file__).resolve().parent
REPO = Path("/scratch/jshen/Github/Pansoma")
RR = REPO / "analysis/reviewer_response_20261011"
TUMORS = ["HG008T", "COLO829T", "HCC1395T"]
PLATS = ["Illumina", "PacBio", "ONT"]
GRAPHS = ["linear_nositef", "d9", "d46"]


def ds(tumor, p):
    return "COLO829T_fiberseq" if (tumor, p) == ("COLO829T", "PacBio") else f"{tumor}_{p}"


def pname(tumor, p):
    return {"Illumina": "Illumina", "ONT": "ONT", "PacBio": "Fiber-seq" if tumor == "COLO829T" else "PacBio HiFi"}[p]


def exact(a):
    """P / R / F1 from the counts (report.json stores them rounded to 4 decimals), as rtg: precision on the calls
    (tp_call), recall on the truth (tp_baseline)."""
    tpc = a.get("tp_call", a["tp_baseline"])
    p = tpc / (tpc + a["fp"]) if tpc + a["fp"] else 0.0
    r = a["tp_baseline"] / (a["tp_baseline"] + a["fn"])
    b = dict(a)
    b.update(precision=p, recall=r, f1=2 * p * r / (p + r) if p + r else 0.0)
    assert abs(b["f1"] - a["f1"]) < 6e-4 and abs(p - a["precision"]) < 6e-4 and abs(r - a["recall"]) < 6e-4, (a, b)
    return b


def rep(d, stage="pon"):
    return exact(json.loads((Path(d) / "eval_bed" / "report.json").read_text())[stage]["at_threshold"])


def rule(kind, p):
    return "" if kind == "SNV" else ("pp_" if p == "ONT" else "ppr8_")


def cell_dir(n, kind, p, graph, tumor, seed=0):
    r = rule(kind, p)
    if graph in ("linear_nositef", "d9"):
        return G / f"multicancer{n}_{p}_{kind}_{graph}" / f"vcf_{r}chr1_{tumor}"
    if kind == "INDEL":
        return R / f"{'release' if n == 6 else 'unseen'}_s{seed}_{p}_INDEL" / f"vcf_{r}chr1_{ds(tumor, p)}_d46"
    assert n == 6
    return N / f"multicancer6_{p}_SNV_newpipe" / f"vcf_chr1_{ds(tumor, p)}_d46"


def prf(a):
    return f"{a['precision']:.3f} / {a['recall']:.3f} / {a['f1']:.3f}"


def truth(a):
    return a["tp_baseline"] + a["fn"]


def mean3(cells):
    return {k: statistics.mean(c[k] for c in cells) for k in ("precision", "recall", "f1")}


def seeds_cell(cells, k):
    v = [c[k] for c in cells]
    return f"{' / '.join(f'{x:.3f}' for x in v)} ({statistics.mean(v):.3f} ± {statistics.stdev(v):.3f})"


def table_after(path, heading):
    lines = Path(path).read_text().splitlines()
    i = next(k for k, x in enumerate(lines) if x.startswith(heading))
    out = []
    for x in lines[i + 1:]:
        if x.startswith("|"):
            out.append([y.strip().replace("**", "") for y in x.strip().strip("|").split("|")])
        elif out:
            break
    return [r for r in out if not set("".join(r)) <= set("-: ")]


def nositef_ris():
    out = {}
    for line in (REPO / "analysis/m6_nositef_graphs_20261118/results.md").read_text().splitlines():
        c = [x.strip() for x in line.strip("|").split("|")]
        if len(c) == 14 and c[2] in ("HPRC v1.1", "HPRC v2.1") and c[0] in TUMORS:
            tp, fp, fn = int(c[6]), int(c[7]), int(c[8])
            out[(c[0], c[1], c[2])] = exact(dict(precision=float(c[9]), recall=float(c[10]), f1=float(c[11]),
                                                 tp_baseline=tp, fn=fn, fp=fp))
    assert len(out) == 18, len(out)
    return out


EN = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators")


def comp_nofilter(caller, kind, tumor, p):
    """Comparator PASS calls without its own population filter (ClairS-TO rerun without the panel-of-normals step,
    DeepSomatic-TO with its PON filter lifted), from the counts."""
    d = "score_nopon_rerun" if caller == "clairs_to" else "score_nopon"
    e = json.loads((EN / d / caller / ds(tumor, p) / kind / "report.json").read_text())
    return exact(e["eval_bed"]["pass"]["at_threshold"])


def dfmt(x):
    """Signed difference with three decimals, Unicode minus, no signed zero."""
    if round(x, 3) == 0:
        return "0.000"
    return f"{x:+.3f}".replace("-", "−")


T = {}
md = ["# Letter ver4: tables from existing reports of the current models (make_tables.py)", ""]


def emit(key, title, header, rows, note=""):
    T[key] = dict(title=title, header=header, rows=rows, note=note)
    md.extend([f"## {key}: {title}", ""] + ([note, ""] if note else []))
    md.append("| " + " | ".join(header) + " |")
    md.append("|" + "---|" * len(header))
    md.extend("| " + " | ".join(str(c).replace("\n", "<br>") for c in r) + " |" for r in rows)
    md.append("")


# ------------------------------------------------------------------ no population filter, three references
rows, facts = [], {}
for tumor in TUMORS:
    for p in PLATS:
        c = [rep(cell_dir(6, "SNV", p, g, tumor), "raw") for g in GRAPHS]
        assert len({truth(x) for x in c}) == 1, (tumor, p)
        cl, dv = comp_nofilter("clairs_to", "SNV", tumor, p), comp_nofilter("deepsomatic", "SNV", tumor, p)
        assert truth(cl) == truth(dv) == truth(c[0]), (tumor, p)
        rows.append([f"{tumor} {pname(tumor, p)}"] + [prf(x) for x in c] + [prf(cl), prf(dv)])
        facts[f"{tumor} {p}"] = [x["f1"] for x in c]
emit("nofilter_snv", "chr1 SNVs, no population filter, six-tumor models on three references (P / R / F1)",
     ["chr1 SNVs, no population filter (P / R / F1)", "Pansoma GRCh38-only", "Pansoma HPRC v1.1", "Pansoma HPRC v2.1", "ClairS-TO", "DeepSomatic-TO"], rows,
     "Seed 0. Comparator columns: the same call sets as in ver3, P / R / F1 recomputed from the counts.")
T["nofilter_snv"]["f1"] = facts

rows, facts = [], {}
for tumor in TUMORS:
    for p in PLATS:
        lin, d9 = (rep(cell_dir(6, "INDEL", p, g, tumor), "raw") for g in ("linear_nositef", "d9"))
        s = [rep(cell_dir(6, "INDEL", p, "d46", tumor, k), "raw") for k in range(3)]
        assert len({truth(x) for x in [lin, d9] + s}) == 1, (tumor, p)
        cl, dv = comp_nofilter("clairs_to", "INDEL", tumor, p), comp_nofilter("deepsomatic", "INDEL", tumor, p)
        assert truth(cl) == truth(dv) == truth(lin), (tumor, p)
        rows.append([f"{tumor} {pname(tumor, p)}", prf(lin), prf(d9), prf(mean3(s)), prf(cl), prf(dv)])
        facts[f"{tumor} {p}"] = [lin["f1"], d9["f1"], statistics.mean(x["f1"] for x in s)]
emit("nofilter_indel", "chr1 INDELs, no population filter, six-tumor models on three references (P / R / F1)",
     ["chr1 INDELs, no population filter (P / R / F1)", "Pansoma GRCh38-only", "Pansoma HPRC v1.1", "Pansoma HPRC v2.1", "ClairS-TO", "DeepSomatic-TO"], rows,
     "GRCh38-only and HPRC v1.1: seed 0; HPRC v2.1: mean of three seeds. PASS = ppr8 (ONT: pp).")
T["nofilter_indel"]["f1"] = facts

# ------------------------------------------------------------------ Response Table 3.1: three references after the filter
RIS = nositef_ris()
rows, t31 = [], {}
for tumor in TUMORS:
    for p in PLATS:
        lin, d9, d46 = (rep(cell_dir(6, "SNV", p, g, tumor)) for g in GRAPHS)
        n1, n2 = RIS[(tumor, pname(tumor, p), "HPRC v1.1")], RIS[(tumor, pname(tumor, p), "HPRC v2.1")]
        assert len({truth(x) for x in (lin, d9, d46, n1, n2)}) == 1, (tumor, p, [truth(x) for x in (lin, d9, d46, n1, n2)])
        for name, a, key in (("GRCh38-only", lin, "linear_nositef"), ("HPRC v1.1", d9, "d9"),
                             ("HPRC v1.1, no site features", n1, "d9_nositef"), ("HPRC v2.1", d46, "d46"),
                             ("HPRC v2.1, no site features", n2, "d46_nositef")):
            d = "–" if key == "linear_nositef" else dfmt(a["f1"] - lin["f1"])
            rows.append([tumor, pname(tumor, p), name, f"{a['precision']:.3f}", f"{a['recall']:.3f}", f"{a['f1']:.3f}", d])
            t31[f"{tumor}|{p}|{key}"] = dict(p=a["precision"], r=a["recall"], f1=a["f1"], tp=a["tp_baseline"], fp=a["fp"], fn=a["fn"])
emit("t31_three_refs", "chr1 SNV after the population filter, six-tumor models on three references",
     ["tumor", "platform", "reference", "precision", "recall", "F1", "ΔF1, HPRC − GRCh38-only"], rows,
     "Seed 0. GRCh38-only models take no site features. 'no site features' HPRC rows: RIS models, results.md. The paired "
     "95 % intervals of ΔF1 (block bootstrap) are added by the build step from the regions item for the rows whose call sets "
     "are on this cluster.")
T["t31_three_refs"]["cells"] = t31

# ------------------------------------------------------------------ Response Table 3.2: per-seed results
rows, sds = [], {"six": {}, "four": {}}
for tumor in TUMORS:
    for p in PLATS:
        for n, name in ((6, "6 tumors"), (4, "4 tumors")):
            c = [rep(cell_dir(n, "INDEL", p, "d46", tumor, k)) for k in range(3)]
            assert len({truth(x) for x in c}) == 1
            rows.append([tumor, pname(tumor, p), name] + [seeds_cell(c, k) for k in ("precision", "recall", "f1")])
            sds["six" if n == 6 else "four"][f"{tumor} {p}"] = statistics.stdev(x["f1"] for x in c)
emit("t32_indel_seeds", "Per-seed chr1 INDEL results of the final models (ONT at pp)",
     ["tumor", "platform", "training set", "precision", "recall", "F1"], rows,
     "PASS = ppr8 on Illumina / PacBio HiFi / Fiber-seq, pp on ONT. Seeds 0 / 1 / 2 (mean ± sample SD).")
T["t32_indel_seeds"]["sd"] = sds

rows, snv_sd = [], {}
for p, pn in (("Illumina", "Illumina"), ("fiberseq", "Fiber-seq")):
    for g, gn, suf in (("d46", "HPRC v2.1", "_d46"), ("d9", "HPRC v1.1", "")):
        c = []
        for k in range(3):
            run = N / f"COLO829_{p}_SNV_{g}_newpipe{'' if k == 0 else f'_s{k}'}"
            c.append(rep(run / f"vcf_chr1_COLO829T_{p}{suf}"))
        assert len({truth(x) for x in c}) == 1
        rows.append(["COLO829T", pn, f"COLO829T only, SNV, {gn}"] + [seeds_cell(c, k) for k in ("precision", "recall", "f1")])
        snv_sd[f"{pn} {gn}"] = statistics.stdev(x["f1"] for x in c)
emit("t32_snv_seeds", "Per-seed chr1 SNV results of the COLO829T single-tumor models of the current pipeline",
     ["tumor", "platform", "training set", "precision", "recall", "F1"], rows, "PASS = rule 2. Seeds 0 / 1 / 2 (mean ± sample SD).")
T["t32_snv_seeds"]["sd"] = snv_sd

# ------------------------------------------------------------------ Response Tables 2.7 / 2.8: nearby germline variants
ST = RR / "experiments/E19_germline_density/strata_latest.md"
for key, sec, title in (("t27_snv_nearby", "## 7.", "chr1 SNV, HG008T, six-tumor model, by nearby germline variants"),
                        ("t28_indel_nearby", "## 1.", "chr1 INDEL, HG008T, six-tumor model, by nearby germline variants")):
    t = table_after(ST, sec)
    rows = [[c.replace("+0.", "+0.").replace("-0.", "−0.") if i == 6 else c for i, c in enumerate(r)] for r in t[1:]]
    emit(key, title, t[0], rows, f"Copied from E19_germline_density/strata_latest.md section {sec[3:]}")
head7 = next(x for x in ST.read_text().splitlines() if x.startswith("## 7."))
head1 = next(x for x in ST.read_text().splitlines() if x.startswith("## 1."))
T["t27_snv_nearby"]["section_title"], T["t28_indel_nearby"]["section_title"] = head7, head1

# ------------------------------------------------------------------ absorbed INDELs on HPRC v2.1
t = table_after(REPO / "analysis/graph_absorbed_somatic_hcc1395t_20261112/tables_v21_patient_frame.md", "## 1. Repeat class")
NAMES = ["Absorbed by the graph (union across platforms)", "Homopolymer ≥ 7 bp", "Short tandem repeat (unit 2–6 bp, ≥ 3 copies)",
         "Homopolymer + short tandem repeat, total", "Other repeats", "Substitution within repetitive sequence", "Others",
         "Non-repetitive sequence"]
assert len(t) == 9 and t[1][0].startswith("absorbed") and t[8][0] == "non-repeat", t
emit("absorbed_v21", "Sequence context of the somatic truth INDELs absorbed by the HPRC v2.1 graph (chr1–22)",
     ["Truth INDELs (chr1–22)", "HG008T", "COLO829T", "HCC1395T"], [[NAMES[i]] + r[1:] for i, r in enumerate(t[1:])],
     "Copied from graph_absorbed_somatic_hcc1395t_20261112/tables_v21_patient_frame.md section 1.")

# ------------------------------------------------------------------ numbers quoted in sentences
F = {}
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")


def comp(caller, kind, tumor, p):
    e = json.loads((E8 / caller / ds(tumor, p) / kind / "report.json").read_text())
    return exact(e["eval_bed"]["pass_pon"]["at_threshold"])


def rng(v, nd=3):
    return f"{min(v):.{nd}f}–{max(v):.{nd}f}"


raw = {g: [rep(cell_dir(6, "SNV", p, g, t), "raw")["f1"] for t in TUMORS for p in PLATS] for g in GRAPHS}
F["snv_raw_f1_range"] = {g: rng(raw[g]) for g in GRAPHS}
F["snv_raw_f1_v21_by_tumor"] = {t: rng([rep(cell_dir(6, "SNV", p, "d46", t), "raw")["f1"] for p in PLATS]) for t in TUMORS}
c31 = T["t31_three_refs"]["cells"]
d_with = [c31[f"{t}|{p}|{g}"]["f1"] - c31[f"{t}|{p}|linear_nositef"]["f1"] for t in ("HG008T", "COLO829T") for p in PLATS for g in ("d9", "d46")]
d_nosf = [c31[f"{t}|{p}|{g}_nositef"]["f1"] - c31[f"{t}|{p}|linear_nositef"]["f1"] for t in ("HG008T", "COLO829T") for p in PLATS for g in ("d9", "d46")]
F["t31_gain_hg_colo_with_sitef"] = rng(d_with)
F["t31_gain_hg_colo_no_sitef"] = rng(d_nosf)
F["t31_hcc_ont_v21_no_sitef"] = round(c31["HCC1395T|ONT|d46_nositef"]["f1"] - c31["HCC1395T|ONT|linear_nositef"]["f1"], 3)
F["t31_hcc_delta"] = {f"{p}|{g}": round(c31[f"HCC1395T|{p}|{g}"]["f1"] - c31[f"HCC1395T|{p}|linear_nositef"]["f1"], 3) for p in PLATS for g in ("d9", "d46")}
# four-line and six-tumor SNV on v2.1 against the comparators (R1.7)
m4 = {}
for t in TUMORS:
    for p in PLATS:
        if p == "ONT":
            d4 = N / "multicancer4_ONT_SNV_newpipe" / f"vcf_chr1_{ds(t, p)}_d46"
        else:
            d4 = R / f"multicancer4_{p}_SNV_sitef_af05_lnf25w1" / f"vcf_chr1_{ds(t, p)}_d46{'_snvaf05' if t == 'COLO829T' else ''}_r095"
        a4, a6 = rep(d4), rep(cell_dir(6, "SNV", p, "d46", t))
        dsv, cl = comp("deepsomatic", "SNV", t, p), comp("clairs_to", "SNV", t, p)
        assert truth(a4) == truth(a6) == truth(dsv) == truth(cl), (t, p)
        m4[f"{t} {p}"] = dict(four=round(a4["f1"], 3), six=round(a6["f1"], 3), deepsomatic=round(dsv["f1"], 3), clairs=round(cl["f1"], 3))
F["snv_v21_vs_comparators"] = m4
# INDEL six-tumor mean of three seeds against the comparators (R3.10)
wins = []
for t in ("HG008T", "COLO829T"):
    for p in PLATS:
        f6 = statistics.mean(rep(cell_dir(6, "INDEL", p, "d46", t, k))["f1"] for k in range(3))
        if f6 > max(comp("deepsomatic", "INDEL", t, p)["f1"], comp("clairs_to", "INDEL", t, p)["f1"]):
            wins.append(f"{t} {pname(t, p)}")
F["indel_six_above_both_hg_colo"] = wins
sd = T["t32_indel_seeds"]["sd"]
allsd = sorted(list(sd["six"].values()) + list(sd["four"].values()))
F["indel_sd_18"] = dict(range=rng(allsd), median=round(statistics.median(allsd), 3))
ssd = sorted(T["t32_snv_seeds"]["sd"].values())
F["snv_sd_4"] = dict(range=rng(ssd), median=round(statistics.median(ssd), 3))
s27 = T["t27_snv_nearby"]["rows"]
dd = [[float(x.replace("−", "-")) for x in r[6].split(" / ")] for r in s27 if r[6] != "–"]
F["t27_delta"] = dict(none=rng([x[0] for x in dd]), one_or_more=rng([y for x in dd for y in x[1:]]))
T["facts"] = F
md.extend(["## facts (numbers quoted in sentences)", "", "```", json.dumps(F, indent=1, ensure_ascii=False), "```", ""])

(HERE / "tables.json").write_text(json.dumps(T, indent=1, ensure_ascii=False) + "\n")
(HERE / "tables.md").write_text("\n".join(md) + "\n")
print("wrote", HERE / "tables.json", HERE / "tables.md", "tables:", ", ".join(T))

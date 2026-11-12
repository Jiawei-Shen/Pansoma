#!/usr/bin/env python3
"""Letter Supplementary Table R1.6a on the CURRENT pipelines (user 2026-11-10: "我们现在用的是AF 0.05的 ... 更新"):
chr1, benchmark regions, complete truth; HG008T / COLO829T / HCC1395T x SNV / INDEL.

Current pipelines (pansoma-net-v2 plan, decided 2026-11-09 / 10-26):
  SNV, PacBio HiFi / Fiber-seq / Illumina: tensors at SNV AF 0.05 (<SET>_snvaf05), site features, rule 2 (call recall 0.95)
      COLO829T  snv_newpipe_20261109/COLO829_<P>_SNV_d46_newpipe{,_s1,_s2}      (single-tumor, three seeds)
      HCC1395T  four-line multicancer4_<P>_SNV_sitef_af05_lnf25w1 (vcf_*_r095) and six-tumor
                snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe              (no single-tumor HCC1395T model)
      HG008T    NOT rebuilt / retrained (user 11-09: no truth below AF 0.1): ris_7fc7396_20261012/HG008_<P>_SNV_d46_normbq{,_s1,_s2}
  SNV, ONT: original recipe (AF 0.08, no site features, PASS = validation recall 0.9)
      HG008T / COLO829T  ris_7fc7396_20261012/<S>_ONT_SNV_d46_normbq{,_s1,_s2}
      HCC1395T           multicancer{4,6}_ONT_SNV/vcf_chr1_HCC1395T_ONT_d46_normfix
  INDEL (all platforms, INDEL floor 0.08 unchanged): final release models, release_s{0,1,2}_<P>_INDEL/vcf_ppr8_chr1_<SET>
Runs still training are left out until their report exists; the table lists what each cell contains.

Candidate rows (exact candidate; exact and not population-tagged), PASS truth with VCF POS in the BED, on the set each
model is tested on: tmp/assessment_20261030/denominators/{label_status.json (HG008T / COLO829T d46),
hcc1395t/label_status.json (HCC1395T d46), snvaf05/label_status.test.json (the four SNV AF 0.05 sets)}.
Truth row: rtg baseline of the Pansoma reports. Comparator rows: ../nopon_comparators/r16a_rows.json (unchanged inputs).
Writes r16a_current.md / .json here.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
RIS, NEW = RUNS / "ris_7fc7396_20261012", RUNS / "snv_newpipe_20261109"
LS = Path("/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/denominators")
TUMORS = {"HG008T": ("Illumina", "PacBio", "ONT"), "COLO829T": ("Illumina", "fiberseq", "ONT"), "HCC1395T": ("Illumina", "PacBio", "ONT")}
S = {"HG008T": "HG008", "COLO829T": "COLO829"}
SEEDS = ("", "_s1", "_s2")


def snv_runs(t, p):
    """(label, report path, candidate set) for every model of tumor t, platform p (the current SNV pipeline)."""
    if t == "HG008T" or p == "ONT" and t == "COLO829T":
        return [(f"{S[t]}_{p}_SNV_d46_normbq{x}", RIS / f"{S[t]}_{p}_SNV_d46_normbq{x}/vcf_chr1/eval_bed/report.json") for x in SEEDS], f"{t}_{p}_d46"
    if t == "COLO829T":
        return [(f"COLO829_{p}_SNV_d46_newpipe{x}", NEW / f"COLO829_{p}_SNV_d46_newpipe{x}/vcf_chr1_COLO829T_{p}_d46/eval_bed/report.json") for x in SEEDS], f"{t}_{p}_d46_snvaf05"
    if p == "ONT":
        return [(f"multicancer{g}_ONT_SNV (normfix)", RUNS / f"multicancer{g}_ONT_SNV/vcf_chr1_HCC1395T_ONT_d46_normfix/eval_bed/report.json") for g in (4, 6)], f"{t}_{p}_d46"
    return [(f"multicancer4_{p}_SNV_sitef_af05_lnf25w1 (r095)", RUNS / f"multicancer4_{p}_SNV_sitef_af05_lnf25w1/vcf_chr1_HCC1395T_{p}_d46_r095/eval_bed/report.json"),
            (f"multicancer6_{p}_SNV_newpipe", NEW / f"multicancer6_{p}_SNV_newpipe/vcf_chr1_HCC1395T_{p}_d46/eval_bed/report.json")], f"{t}_{p}_d46_snvaf05"


def indel_runs(t, p):
    m = "PacBio" if p == "fiberseq" else p  # the PacBio HiFi release model is the one applied to Fiber-seq
    return [(f"release_s{s}_{m}_INDEL (ppr8)", RUNS / f"release_s{s}_{m}_INDEL/vcf_ppr8_chr1_{t}_{p}_d46/eval_bed/report.json") for s in (0, 1, 2)], f"{t}_{p}_d46"


def rng(v):
    lo, hi = f"{min(v):.3f}", f"{max(v):.3f}"
    return lo if lo == hi else f"{lo}–{hi}"


ls = {}
for f in ("label_status.json", "hcc1395t/label_status.json", "snvaf05/label_status.test.json"):
    ls.update(json.loads((LS / f).read_text()))
comp = json.loads((HERE.parent / "nopon_comparators/r16a_rows.json").read_text())
out, used, pending = {}, [], []
for t, plats in TUMORS.items():
    for kind, fn in (("SNV", snv_runs), ("INDEL", indel_runs)):
        raw, pon, n, ex, un, csets = [], [], set(), [], [], []
        for p in plats:
            runs, cset = fn(t, p)
            for lab, f in runs:
                if not f.exists():
                    pending.append(f"{t} {kind} {p}: {lab}")
                    continue
                r = json.loads(f.read_text())
                a, b = r["raw"]["at_threshold"], r["pon"]["at_threshold"]
                # recall from the counts (the report's recall is already rounded to 4 decimals: no double rounding)
                raw.append(a["tp_baseline"] / (a["tp_baseline"] + a["fn"])); pon.append(b["tp_baseline"] / (b["tp_baseline"] + b["fn"]))
                n.add(a["tp_baseline"] + a["fn"])
                used.append(f"{t} {kind} {p}: {lab}")
            c = ls[f"{cset}|{kind}"]["chr1"]
            ex.append(c["all"].get("exact", 0) / c["all"]["n"]); un.append(c["untagged"].get("exact", 0) / c["all"]["n"])
            csets.append(f"{cset} (n {c['all']['n']})")
        assert len(n) == 1, (t, kind, n)
        out[f"{t}|{kind}"] = dict(truth=n.pop(), exact=ex, untagged=un, raw=raw, pon=pon, candidate_sets=csets)
COLS = [(t, k) for k in ("SNV", "INDEL") for t in TUMORS]
rows = [("Truth alleles", lambda t, k: f"{out[f'{t}|{k}']['truth']:,}"),
        ("Exact candidate", lambda t, k: rng(out[f"{t}|{k}"]["exact"])),
        ("Exact candidate, not population-tagged", lambda t, k: rng(out[f"{t}|{k}"]["untagged"])),
        ("Pansoma PASS, before filtering", lambda t, k: rng(out[f"{t}|{k}"]["raw"])),
        ("Pansoma PASS, after filtering", lambda t, k: rng(out[f"{t}|{k}"]["pon"]))]


def comp_cell(c, t, k, a, b):
    x = [comp[f"{c}|{t}_{p}|{k}"] for p in TUMORS[t]]
    u, v = rng([y[a]["tp"] / y[a]["baseline"] for y in x]), rng([y[b]["tp"] / y[b]["baseline"] for y in x])
    return f"{u} (unchanged)" if u == v else f"{u} → {v}"


for c, name in (("deepsomatic", "DeepSomatic-TO"), ("clairs_to", "ClairS-TO")):
    for a, b, lab in (("nopon", "native", "no PoN → native"), ("native", "native_ours", "native → common filter")):
        rows.append((f"{name}, {lab}", lambda t, k, c=c, a=a, b=b: comp_cell(c, t, k, a, b)))
L = ["# Supplementary Table R1.6a on the current pipelines (generated by r16a_current.py)", "",
     "| Recall stage | " + " | ".join(f"{t} {k}" for t, k in COLS) + " |", "|---|" + "---|" * len(COLS)]
L += [f"| {r} | " + " | ".join(fn(t, k) for t, k in COLS) + " |" for r, fn in rows]
L += ["", "## Models in the Pansoma rows", ""] + [f"- {u}" for u in used]
L += ["", "## Not yet available (left out)", ""] + ([f"- {u}" for u in pending] or ["- none"])
L += ["", "## Candidate sets (exact-candidate rows)", ""] + [f"- {t} {k}: " + ", ".join(out[f'{t}|{k}']['candidate_sets']) for t, k in COLS]
(HERE / "r16a_current.json").write_text(json.dumps(dict(cells=out, used=used, pending=pending), indent=1) + "\n")
(HERE / "r16a_current.md").write_text("\n".join(L) + "\n")
print("\n".join(L))

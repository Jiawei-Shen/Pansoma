"""HG008-T: PacBio v6 vs ONT-UL vs Illumina — recall, per-truth patterns, miss classes, labels (writes comparison_tables.md).

Inputs: each run's somatic.recall.tsv and labels.manifest.json, the per-platform miss tables of
tmp/somatic_miss_analysis_{v6,ont,illumina}, and the label ndjson (sites whose representative allele is germline
while another allele is a somatic truth).
"""
import csv, glob, json, sys
from collections import Counter, defaultdict
from pathlib import Path

csv.field_size_limit(sys.maxsize)
D = Path("/scratch/jshen/data/pansoma_v2_tensors")
TMP = Path("/scratch/jshen/Github/Pansoma/tmp")
P = {"PacBio": (D / "Liss_lab_PacBio_Revio_20240125/v6_tensors", TMP / "somatic_miss_analysis_v6"),
     "ONT": (D / "Liss_lab_Northeastern-ONT-UL-20241216/v3_tensors", TMP / "somatic_miss_analysis_ont"),
     "Illumina": (D / "Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors", TMP / "somatic_miss_analysis_illumina")}
NAMES = list(P)
out = []
pr = out.append
pct = lambda a, b: f"{100 * a / b:.1f}%" if b else "-"  # noqa: E731

recall = {n: {r["truth_id"]: r for r in csv.DictReader(open(t / "somatic.recall.tsv"), delimiter="\t")} for n, (t, _) in P.items()}
ids = list(recall["PacBio"])
kind = {t: recall["PacBio"][t]["kind"] for t in ids}
pib = {t: recall["PacBio"][t]["passed"] == "True" and recall["PacBio"][t]["in_bed"] == "True" for t in ids}
grp = lambda s: "tensor" if s.startswith("tensor") else s  # noqa: E731

# 1. recall
pr("## 1. 召回（全部 somatic allele / PASS 且在 somatic BED 内）\n")
pr("| 状态 | 种类 | " + " | ".join(NAMES) + " |\n|---|---|" + "---:|" * 3)
for scope, sel in (("全部", lambda t: True), ("PASS∩BED", lambda t: pib[t])):
    for k in ("SNP", "DEL", "INS"):
        sub = [t for t in ids if kind[t] == k and sel(t)]
        for s, lab in (("tensor_representative", "代表 allele"), ("tensor_non_representative_allele", "非代表 allele"),
                       ("filtered", "filtered"), ("no_candidate", "没有候选")):
            cells = [f"{sum(recall[n][t]['status'] == s for t in sub):,} ({pct(sum(recall[n][t]['status'] == s for t in sub), len(sub))})" for n in NAMES]
            pr(f"| {scope} {lab} | {k} ({len(sub):,}) | " + " | ".join(cells) + " |")

# 2. per-truth pattern across the three platforms
pr("\n## 2. 同一个 truth 在三个平台上的状态（全部 somatic allele）\n")
pr("「T」= 有 tensor（代表或非代表），「F」= filtered，「N」= 没有候选；顺序 PacBio / ONT / Illumina。\n")
pat = defaultdict(Counter)
for t in ids:
    pat[kind[t]]["/".join({"tensor": "T", "filtered": "F", "no_candidate": "N"}.get(grp(recall[n][t]["status"]), "?") for n in NAMES)] += 1
keys = sorted(set().union(*pat.values()), key=lambda x: (-sum(pat[k][x] for k in pat), x))
pr("| PacBio/ONT/Illumina | SNP | DEL | INS |\n|---|---:|---:|---:|")
for x in keys:
    pr(f"| {x} | {pat['SNP'][x]} | {pat['DEL'][x]} | {pat['INS'][x]} |")
tot = {k: sum(v.values()) for k, v in pat.items()}
for label, test in (("三个平台都有 tensor", lambda x: x == "T/T/T"), ("至少一个平台有 tensor", lambda x: "T" in x),
                    ("三个平台都没有候选", lambda x: x == "N/N/N")):
    pr(f"| **{label}** | " + " | ".join(f"**{sum(pat[k][x] for x in keys if test(x))}** ({pct(sum(pat[k][x] for x in keys if test(x)), tot[k])})" for k in ("SNP", "DEL", "INS")) + " |")

# 3. classes of truths missed everywhere (by each platform's own class)
cls = {}
for n, (_, a) in P.items():
    for f in ("indel_no_candidate.tsv", "snv_no_candidate.tsv"):
        for r in csv.DictReader(open(a / f), delimiter="\t"):
            cls[(n, r["chrom"], r["vcf_pos"], r["vcf_ref"], r["vcf_alt"])] = r["final_class"].split("_")[0]
truth = {r["truth_id"]: r for r in csv.DictReader(open(D / "truth/somatic.graph.tsv"), delimiter="\t")}
pr("\n## 3. 三个平台都没有候选的 truth：按 PacBio 的类别，看另两个平台是否同类\n")
agree = Counter(); by = Counter()
for t in ids:
    if all(recall[n][t]["status"] == "no_candidate" for n in NAMES):
        r = truth[t]; c = [cls.get((n, r["chrom"], r["vcf_pos"], r["vcf_ref"], r["vcf_alt"]), "?") for n in NAMES]
        by[(kind[t], c[0])] += 1
        agree[(kind[t], c[0], "三个平台同类" if len(set(c)) == 1 else "类别不同")] += 1
pr("| 种类 | PacBio 类别 | 个数 | 三个平台同类 |\n|---|---|---:|---:|")
for (k, c), n in sorted(by.items(), key=lambda x: (x[0][0], -x[1])):
    pr(f"| {k} | {c} | {n} | {agree[(k, c, '三个平台同类')]} |")

# 4. miss classes side by side
pr("\n## 4. 没有候选的原因分类（每个平台自己的未命中）\n")
cc = {n: Counter() for n in NAMES}
for (n, *_), c in cls.items():
    pass
for n, (_, a) in P.items():
    for f in ("indel_no_candidate.tsv", "snv_no_candidate.tsv"):
        for r in csv.DictReader(open(a / f), delimiter="\t"):
            cc[n][(r["kind"], r["final_class"])] += 1
allc = sorted({k for c in cc.values() for k in c}, key=lambda x: ({"INS": 0, "DEL": 1, "SNP": 2}[x[0]], x[1][:1], int("".join(ch for ch in x[1].split("_")[0][1:] if ch.isdigit()) or 0), x[1]))
pr("| 种类 | 类别 | " + " | ".join(NAMES) + " |\n|---|---|" + "---:|" * 3)
for k, c in allc:
    pr(f"| {k} | {c} | " + " | ".join(str(cc[n][(k, c)]) for n in NAMES) + " |")
for k in ("INS", "DEL", "SNP"):
    pr(f"| {k} | **合计** | " + " | ".join(f"**{sum(v for (kk, _), v in cc[n].items() if kk == k)}**" for n in NAMES) + " |")

# 5. labels
pr("\n## 5. 标签（truth-labels-v2）\n")
ORDER = [("somatic", "representative_allele_is_somatic_truth"), ("germline", "representative_allele_is_germline_truth"),
         ("non", "confident_no_truth_allele"), ("ignore", "outside_confident_region"), ("ignore", "not_on_unique_grch38_node"),
         ("ignore", "near_truth_allele_mismatch"), ("ignore", "truth_matches_non_representative_allele"),
         ("ignore", "germline_truth_filtered"), ("ignore", "somatic_truth_filtered")]
VAL = {"somatic": 1, "germline": 2, "non": 0, "ignore": -1}
for k in ("SNV", "INDEL"):
    man = {n: json.load(open(t / k / "labels.manifest.json")) for n, (t, _) in P.items()}
    pr(f"\n**{k}**\n\n| 标签 | reason | " + " | ".join(NAMES) + " |\n|---:|---|" + "---:|" * 3)
    for name, reason in ORDER:
        pr(f"| {VAL[name]} | {reason} | " + " | ".join(f"{man[n]['reasons'].get(reason, 0):,} ({pct(man[n]['reasons'].get(reason, 0), man[n]['tensors'])})" for n in NAMES) + " |")
    pr("| | 合计 | " + " | ".join(f"{man[n]['tensors']:,}" for n in NAMES) + " |")
    pr("| | 0 : 1 | " + " | ".join(f"{man[n]['totals'].get('non', 0) / max(1, man[n]['totals'].get('somatic', 0)):.0f} : 1" for n in NAMES) + " |")

# 6. label-quality items
pr("\n## 6. 标签质量问题\n")
rows = []
for n, (t, a) in P.items():
    res = list(csv.DictReader(open(a / "residual_edits_labelled_germline.tsv"), delimiter="\t"))
    snv = list(csv.DictReader(open(a / "snv_no_candidate.tsv"), delimiter="\t"))
    both = 0
    for f in glob.glob(f"{t}/*/chr*_labels.ndjson"):
        with open(f) as fh:
            for line in fh:
                if '"somatic": []' in line or '"label": 2,' not in line:
                    continue
                l = json.loads(line)
                if any(g["representative"] for g in l["germline"]) and not any(s["representative"] for s in l["somatic"]):
                    both += 1
    rows.append((n, len({r["residual_candidate"] for r in res}), sum(r["germline_same_allele"] == "True" for r in snv), both))
pr("| | " + " | ".join(NAMES) + " |\n|---|" + "---:|" * 3)
pr("| 残余 edit 的 tensor 被标成 germline（不同候选数） | " + " | ".join(str(r[1]) for r in rows) + " |")
pr("| SNV 未命中里 germline 真值在同一位置、同一 ALT | " + " | ".join(str(r[2]) for r in rows) + " |")
pr("| site 的代表 allele 是 germline、另一个 allele 是 somatic truth → 标 2 | " + " | ".join(str(r[3]) for r in rows) + " |")
Path(__file__).with_name("comparison_tables.md").write_text("\n".join(out) + "\n")
print("\n".join(out))

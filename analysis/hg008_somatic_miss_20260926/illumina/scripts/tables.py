"""Print the summary tables of REPORT.md (markdown) from the v6 recall report, transitions.tsv and the miss tables."""
import csv, json, sys
from collections import Counter, defaultdict
from pathlib import Path
csv.field_size_limit(sys.maxsize)
HERE = Path(__file__).resolve().parent
RUN = Path("/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313")
PB = HERE.parent / "somatic_miss_analysis_v6"
ONT = HERE.parent / "somatic_miss_analysis_ont"
STATUS = ["tensor_representative", "tensor_non_representative_allele", "filtered", "no_candidate", "complex_allele",
          "no_unique_grch38_node_key"]
KINDS = ["SNP", "DEL", "INS"]
pct = lambda a, b: f"{100 * a / b:.1f}%" if b else "-"  # noqa: E731


def recall_table(path, title):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    print(f"\n### {title}\n")
    for scope, sel in (("全部 somatic 等位基因", lambda r: True),
                       ("PASS 且在 somatic BED 内", lambda r: r["passed"] == "True" and r["in_bed"] == "True")):
        c = defaultdict(Counter)
        for r in rows:
            if sel(r):
                c[r["kind"]][r["status"]] += 1
        print(f"\n**{scope}**\n")
        print("| 状态 | " + " | ".join(KINDS) + " |\n|---|" + "---:|" * len(KINDS))
        for s in STATUS:
            if any(c[k][s] for k in KINDS):
                print(f"| {s} | " + " | ".join(f"{c[k][s]:,} ({pct(c[k][s], sum(c[k].values()))})" for k in KINDS) + " |")
        print("| 合计 | " + " | ".join(f"{sum(c[k].values()):,}" for k in KINDS) + " |")


def transitions():
    rows = list(csv.DictReader(open(HERE / "transitions.tsv"), delimiter="\t"))
    short = {"tensor_representative": "rep", "tensor_non_representative_allele": "non-rep", "filtered": "filtered",
             "no_candidate": "no_cand", "complex_allele": "complex", "no_unique_grch38_node_key": "no_key"}
    for source, name in (("ont", "ONT"), ("pacbio", "PacBio v6")):
        for kind in ("SNP", "DEL", "INS"):
            c = Counter((r[f"{source}_status"], r["illumina_status"]) for r in rows if r["kind"] == kind)
            cols = [s for s in STATUS if any(c[(a, s)] for a in STATUS)]
            print(f"\n**{kind}：{name} 状态（行）→ Illumina 状态（列）**\n")
            print(f"| {name} \\ Illumina | " + " | ".join(short[s] for s in cols) + " | 合计 |\n|---|" + "---:|" * (len(cols) + 1))
            for a in STATUS:
                n = sum(c[(a, s)] for s in cols)
                if n:
                    print(f"| {short[a]} | " + " | ".join(str(c[(a, s)]) for s in cols) + f" | {n} |")
    cols = ["tensor_representative", "tensor_non_representative_allele", "filtered", "no_candidate"]
    for source, name in (("ont", "ONT"), ("pacbio", "PacBio v6")):
        print(f"\n**{name} 没有候选的 truth，按 {name} 类别看 Illumina 状态**\n")
        c = Counter((r["kind"], r[f"{source}_class"].split("_")[0], r["illumina_status"]) for r in rows if r[f"{source}_class"])
        classes = sorted({r[f"{source}_class"].split("_")[0] for r in rows if r[f"{source}_class"]},
                         key=lambda x: (x[0], int(x[1:].rstrip("abc") or 0), x))
        print(f"| {name} 类别 | 种类 | " + " | ".join(short[s] for s in cols) + " |\n|---|---|" + "---:|" * 4)
        for kind in ("SNP", "DEL", "INS"):
            for cl in classes:
                n = [c[(kind, cl, s)] for s in cols]
                if sum(n):
                    print(f"| {cl} | {kind} | " + " | ".join(map(str, n)) + " |")


def miss_classes():
    rows = list(csv.DictReader(open(HERE / "indel_no_candidate.tsv"), delimiter="\t"))
    ont = list(csv.DictReader(open(ONT / "indel_no_candidate.tsv"), delimiter="\t"))
    pb = list(csv.DictReader(open(PB / "indel_no_candidate.tsv"), delimiter="\t"))  # PacBio v6
    n, no, np_ = (Counter(r["kind"] for r in x) for x in (rows, ont, pb))
    c = Counter((r["kind"], r["final_class"][:3].rstrip("_")) for r in rows)
    co = Counter((r["kind"], r["final_class"][:3].rstrip("_")) for r in ont)
    cp = Counter((r["kind"], r["final_class"][:3].rstrip("_")) for r in pb)
    classes = sorted({r["final_class"] for r in rows}, key=lambda x: int(x[1:].split("_")[0]))
    print(f"\n### Illumina 没有候选的 INDEL：原因分类（Illumina INS {n['INS']:,}、DEL {n['DEL']:,}；"
          f"ONT INS {no['INS']:,}、DEL {no['DEL']:,}；PacBio v6 INS {np_['INS']:,}、DEL {np_['DEL']:,}）\n")
    print("| 类别 | INS Illumina | INS ONT | INS PacBio | DEL Illumina | DEL ONT | DEL PacBio | ALT 比例中位数 (Illumina) | "
          "按无误 ALT reads 判定 (Illumina) | PASS∩BED (Illumina) |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for cl in classes:
        k = cl[:3].rstrip("_")
        sub = [r for r in rows if r["final_class"] == cl]
        fr = sorted(float(r["alt_fraction"]) for r in sub if r["alt_fraction"])
        print(f"| {cl} | {c[('INS', k)]} ({pct(c[('INS', k)], n['INS'])}) | {co[('INS', k)]} | {cp[('INS', k)]} | "
              f"{c[('DEL', k)]} ({pct(c[('DEL', k)], n['DEL'])}) | {co[('DEL', k)]} | {cp[('DEL', k)]} | "
              f"{fr[len(fr) // 2] if fr else '-'} | {sum(r['from_exact_reads'] == 'True' for r in sub)} | "
              f"{sum(r['passed'] == 'True' and r['in_bed'] == 'True' for r in sub)} |")
    print(f"| 合计 | {n['INS']} | {no['INS']} | {np_['INS']} | {n['DEL']} | {no['DEL']} | {np_['DEL']} | | | |")
    print("\n**更细的原因（reads 的多数表示）**\n")
    c = Counter((r["kind"], r["reason"]) for r in rows)
    print("| reason | INS | DEL |\n|---|---:|---:|")
    for reason in sorted({r["reason"] for r in rows}, key=lambda x: -(c[("INS", x)] + c[("DEL", x)])):
        print(f"| {reason} | {c[('INS', reason)]} | {c[('DEL', reason)]} |")
    # length
    print("\n**长度（Illumina）**\n")
    tr = list(csv.DictReader(open(HERE / "transitions.tsv"), delimiter="\t"))
    truth = {r["truth_id"]: r for r in csv.DictReader(open(RUN / "truth/somatic.graph.tsv"), delimiter="\t")}
    b = lambda n: "1" if n == 1 else "2-5" if n <= 5 else "6-20" if n <= 20 else "21-50" if n <= 50 else ">50"  # noqa: E731
    L = Counter()
    for r in tr:
        t = truth[r["truth_id"]]; n = max(len(t["ref"]), len(t["alt"]))
        if r["kind"] == "SNP":
            continue
        st = r["illumina_status"]
        grp = "有 tensor" if st.startswith("tensor") else "filtered" if st == "filtered" else "没有候选" if st == "no_candidate" else "其它"
        L[(r["kind"], grp, b(n))] += 1
    bs = ["1", "2-5", "6-20", "21-50", ">50"]
    print("| | " + " | ".join(bs) + " |\n|---|" + "---:|" * 5)
    for k in ("INS", "DEL"):
        for grp in ("有 tensor", "filtered", "没有候选"):
            print(f"| {k} {grp} | " + " | ".join(str(L[(k, grp, x)]) for x in bs) + " |")


def residual_labels():
    rows = list(csv.DictReader(open(HERE / "indel_no_candidate.tsv"), delimiter="\t"))
    print("\n### 残余 edits（≥3 条 ALT reads）的 tensor 标签，按位点计（一个位点可有多个标签）\n")
    cats = ["somatic", "germline", "non", "ignore", "no tensor", "(无 ≥3 reads 的残余 edit)"]
    c = defaultdict(Counter)
    for r in rows:
        labs = {x.rsplit("=", 1)[1].split(":")[0].split("(")[0] for x in r["residual_edits(reads)=tensor_label"].split("; ") if x}
        labs = labs or {"(无 ≥3 reads 的残余 edit)"}
        for lab in labs:
            c[(r["kind"], r["final_class"][:3].rstrip("_"))][lab] += 1
    print("| 类别 | " + " | ".join(cats) + " |\n|---|" + "---:|" * len(cats))
    for key in sorted(c):
        print(f"| {key[0]} {key[1]} | " + " | ".join(str(c[key][x]) for x in cats) + " |")


def controls():
    rows = list(csv.DictReader(open(HERE / "controls.tsv"), delimiter="\t"))
    c = Counter((r["kind"], r["reason"]) for r in rows)
    print("\n### 对照（有 tensor 的 truth，同样方法分析）\n")
    print("| reason | SNP | INS | DEL |\n|---|---:|---:|---:|")
    for reason in sorted({r["reason"] for r in rows}, key=lambda x: -(c[("SNP", x)] + c[("INS", x)] + c[("DEL", x)])):
        print(f"| {reason} | {c[('SNP', reason)]} | {c[('INS', reason)]} | {c[('DEL', reason)]} |")


def snv():
    rows = list(csv.DictReader(open(HERE / "snv_no_candidate.tsv"), delimiter="\t"))
    ont = list(csv.DictReader(open(ONT / "snv_no_candidate.tsv"), delimiter="\t"))
    pb = list(csv.DictReader(open(PB / "snv_no_candidate.tsv"), delimiter="\t"))
    c, o, p = (Counter(r["final_class"] for r in x) for x in (rows, ont, pb))
    print(f"\n### SNV 没有候选（Illumina {len(rows)}；ONT {len(ont)}；PacBio v6 {len(pb)}）\n")
    print("| 类别 | Illumina | Illumina 占比 | Illumina ALT 比例中位数 | ONT | PacBio v6 |\n|---|---:|---:|---:|---:|---:|")
    for k in sorted(set(c) | set(o) | set(p)):
        fr = sorted(float(r["alt_fraction"]) for r in rows if r["final_class"] == k and r["alt_fraction"])
        print(f"| {k} | {c[k]} | {pct(c[k], len(rows))} | {fr[len(fr) // 2] if fr else '-'} | {o[k]} | {p[k]} |")


if __name__ == "__main__":
    what = sys.argv[1:] or ["recall", "transitions", "classes", "residual", "controls", "snv"]
    if "recall" in what:
        recall_table(RUN / "v3_tensors/somatic.recall.tsv", "Illumina somatic 召回")
    if "transitions" in what:
        print("\n### ONT / PacBio v6 → Illumina 状态变化")
        transitions()
    if "classes" in what:
        miss_classes()
    if "residual" in what:
        residual_labels()
    if "controls" in what:
        controls()
    if "snv" in what:
        snv()

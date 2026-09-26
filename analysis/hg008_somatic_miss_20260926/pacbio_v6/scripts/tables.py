"""Print the summary tables of REPORT.md (markdown) from the v6 recall report, transitions.tsv and the miss tables."""
import csv, json, sys
from collections import Counter, defaultdict
from pathlib import Path
csv.field_size_limit(sys.maxsize)
HERE = Path(__file__).resolve().parent
RUN = Path("/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125")
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
    for kind in ("DEL", "INS"):
        c = Counter((r["v5_status"], r["v6_status"]) for r in rows if r["kind"] == kind)
        cols = [s for s in STATUS if any(c[(a, s)] for a in STATUS)]
        print(f"\n**{kind}：v5 状态（行）→ v6 状态（列）**\n")
        print("| v5 \\ v6 | " + " | ".join(short[s] for s in cols) + " | 合计 |\n|---|" + "---:|" * (len(cols) + 1))
        for a in STATUS:
            n = sum(c[(a, s)] for s in cols)
            if n:
                print(f"| {short[a]} | " + " | ".join(str(c[(a, s)]) for s in cols) + f" | {n} |")
    print("\n**v5 没有候选的 INDEL，按 v5 类别看 v6 状态**\n")
    c = Counter((r["kind"], r["v5_class"][:2], r["v6_status"]) for r in rows if r["v5_class"])
    cols = ["tensor_representative", "tensor_non_representative_allele", "filtered", "no_candidate"]
    print("| v5 类别 | " + " | ".join(f"{k} {short[s]}" for k in ("DEL", "INS") for s in cols) + " |\n|---|" + "---:|" * 8)
    for cl in ("I1", "I2", "I3", "I4", "I5", "I6", "I7"):
        print(f"| {cl} | " + " | ".join(str(c[(k, cl, s)]) for k in ("DEL", "INS") for s in cols) + " |")
    print("| 合计 | " + " | ".join(str(sum(c[(k, cl, s)] for cl in ("I1", "I2", "I3", "I4", "I5", "I6", "I7")))
                                  for k in ("DEL", "INS") for s in cols) + " |")


def miss_classes():
    rows = list(csv.DictReader(open(HERE / "indel_no_candidate.tsv"), delimiter="\t"))
    v5 = list(csv.DictReader(open(HERE / "v5_reference/v5_indel_no_candidate.tsv"), delimiter="\t"))
    n6 = Counter(r["kind"] for r in rows); n5 = Counter(r["kind"] for r in v5)
    c6 = Counter((r["kind"], r["final_class"]) for r in rows); c5 = Counter((r["kind"], r["final_class"][:2]) for r in v5)
    classes = sorted({r["final_class"] for r in rows}, key=lambda x: int(x[1:].split("_")[0]))
    print(f"\n### v6 没有候选的 INDEL：原因分类（v6 INS {n6['INS']:,}、DEL {n6['DEL']:,}；v5 INS {n5['INS']:,}、DEL {n5['DEL']:,}）\n")
    print("| 类别 | INS v6 | INS v5 | DEL v6 | DEL v5 | ALT 比例中位数 (v6) | 按无误 ALT reads 判定 (v6) | PASS∩BED (v6) |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    for cl in classes:
        sub = [r for r in rows if r["final_class"] == cl]
        fr = sorted(float(r["alt_fraction"]) for r in sub if r["alt_fraction"])
        print(f"| {cl} | {c6[('INS', cl)]} ({pct(c6[('INS', cl)], n6['INS'])}) | {c5[('INS', cl[:2])]} | "
              f"{c6[('DEL', cl)]} ({pct(c6[('DEL', cl)], n6['DEL'])}) | {c5[('DEL', cl[:2])]} | "
              f"{fr[len(fr) // 2] if fr else '-'} | {sum(r['from_exact_reads'] == 'True' for r in sub)} | "
              f"{sum(r['passed'] == 'True' and r['in_bed'] == 'True' for r in sub)} |")
    print(f"| 合计 | {n6['INS']} | {n5['INS']} | {n6['DEL']} | {n5['DEL']} | | | |")
    print("\n**更细的原因（reads 的多数表示）**\n")
    c = Counter((r["kind"], r["reason"]) for r in rows)
    print("| reason | INS | DEL |\n|---|---:|---:|")
    for reason in sorted({r["reason"] for r in rows}, key=lambda x: -(c[("INS", x)] + c[("DEL", x)])):
        print(f"| {reason} | {c[('INS', reason)]} | {c[('DEL', reason)]} |")
    # length
    print("\n**长度（v6）**\n")
    tr = list(csv.DictReader(open(HERE / "transitions.tsv"), delimiter="\t"))
    truth = {r["truth_id"]: r for r in csv.DictReader(open("/scratch/jshen/data/pansoma_v2_tensors/truth/somatic.graph.tsv"), delimiter="\t")}
    b = lambda n: "1" if n == 1 else "2-5" if n <= 5 else "6-20" if n <= 20 else "21-50" if n <= 50 else ">50"  # noqa: E731
    L = Counter()
    for r in tr:
        t = truth[r["truth_id"]]; n = max(len(t["ref"]), len(t["alt"]))
        grp = "有 tensor" if r["v6_status"].startswith("tensor") else "filtered" if r["v6_status"] == "filtered" else "没有候选" if r["v6_status"] == "no_candidate" else "其它"
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
    print("\n### 对照（有 tensor 的 INDEL truth，同样方法分析）\n")
    print("| reason | INS | DEL |\n|---|---:|---:|")
    for reason in sorted({r["reason"] for r in rows}, key=lambda x: -(c[("INS", x)] + c[("DEL", x)])):
        print(f"| {reason} | {c[('INS', reason)]} | {c[('DEL', reason)]} |")


def snv():
    rows = list(csv.DictReader(open(HERE / "snv_no_candidate.tsv"), delimiter="\t"))
    c = Counter(r["final_class"] for r in rows)
    print(f"\n### SNV 没有候选（{len(rows)}；分类沿用 v5 分析，SNV tensors 两版相同）\n")
    print("| 类别 | 个数 | 占比 |\n|---|---:|---:|")
    for k, n in sorted(c.items()):
        print(f"| {k} | {n} | {pct(n, len(rows))} |")


if __name__ == "__main__":
    what = sys.argv[1:] or ["recall", "transitions", "classes", "residual", "controls", "snv"]
    if "recall" in what:
        recall_table(RUN / "v6_tensors/somatic.recall.tsv", "v6 somatic 召回（新标签规则；召回与规则无关）")
    if "transitions" in what:
        print("\n### v5 → v6 状态变化（INDEL）")
        transitions()
    if "classes" in what:
        miss_classes()
    if "residual" in what:
        residual_labels()
    if "controls" in what:
        controls()
    if "snv" in what:
        snv()

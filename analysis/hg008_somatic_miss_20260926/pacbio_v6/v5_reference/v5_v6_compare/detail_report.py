"""Print markdown tables from detail_table.pkl."""
import collections, pickle
import sys
V = sys.argv[1:]; D = {v: pickle.load(open(f"detail_{v}.pkl", "rb")) for v in V}
T = D[V[0]]["truths"]; R = {v: D[v]["result"] for v in V}
CATS = ["somatic (1)", "germline (2)", "non (0)", "−1 near_truth_allele_mismatch", "−1 truth_matches_non_representative_allele",
        "−1 outside_confident_region", "−1 not_on_unique_grch38_node", "−1 other", "only SNV tensor", "no related tensor"]
def best(hits):
    ind = [h for h in hits if h[2] == "INDEL"]
    if not ind: return "only SNV tensor" if hits else "no related tensor"
    labels = {h[0] for h in ind}
    for lab, name in ((1, "somatic (1)"), (2, "germline (2)"), (0, "non (0)")):
        if lab in labels: return name
    reasons = {h[1] for h in ind}
    for r in ("near_truth_allele_mismatch", "truth_matches_non_representative_allele", "outside_confident_region", "not_on_unique_grch38_node"):
        if r in reasons: return "−1 " + r
    return "−1 other"
def table(title, select):
    idx = [i for i, t in enumerate(T) if select(t)]
    kinds = ("DEL", "INS")
    n = {k: sum(T[i]["kind"] == k for i in idx) for k in kinds}
    print(f"\n### {title}  (DEL {n['DEL']}, INS {n['INS']})\n")
    print("| 最好的相关 tensor 标签 | " + " | ".join(f"{k} {v}" for k in kinds for v in V) + " |")
    print("|---|" + "---:|" * (len(kinds) * len(V)))
    c = {(k, v): collections.Counter(best(R[v]["per_truth"][i]) for i in idx if T[i]["kind"] == k) for k in kinds for v in V}
    for cat in CATS:
        print(f"| {cat} | " + " | ".join(str(c[(k, v)][cat]) for k in kinds for v in V) + " |")
    print("| **有相关 INDEL tensor** | " + " | ".join(f"**{n[k] - c[(k, v)]['only SNV tensor'] - c[(k, v)]['no related tensor']}**" for k in kinds for v in V) + " |")
def tensor_table(title, select):
    idx = [i for i, t in enumerate(T) if select(t)]
    print(f"\n### {title}：相关 INDEL tensors（去重）按标签\n")
    print("| 标签 / reason | " + " | ".join(V) + " |"); print("|---|" + "---:|" * len(V))
    cnt = {}
    for v in V:
        seen = {}
        for i in idx:
            for h in R[v]["per_truth"][i]:
                if h[2] == "INDEL": seen[(h[0], h[1], h[3], h[4], i)] = (h[0], h[1])
        cnt[v] = collections.Counter(seen.values())
    for key in sorted(set().union(*cnt.values()), key=lambda k: (-k[0] if k[0] > 0 else 5 + (k[0] == -1), k[1])):
        print(f"| {key[0]} {key[1]} | " + " | ".join(str(cnt[v][key]) for v in V) + " |")
def global_table():
    print("\n### 全部 tensors：0 和 −1 按 reason\n")
    for tk in ("SNV", "INDEL"):
        print(f"\n**{tk}**\n")
        print("| 标签 | reason | " + " | ".join(V) + " |"); print("|---|---|" + "---:|" * len(V))
        agg = {v: collections.Counter() for v in V}
        for v in V:
            for (k, lab, reason, sub), n in R[v]["glob"].items():
                if k == tk: agg[v][(lab, reason)] += n
        for key in sorted(set().union(*agg.values()), key=lambda k: ({0: 0, -1: 1, 1: 2, 2: 3}[k[0]], k[1])):
            print(f"| {key[0]} | {key[1]} | " + " | ".join(f"{agg[v][key]:,}" for v in V) + " |")
        print("| | **合计** | " + " | ".join(f"**{sum(agg[v].values()):,}**" for v in V) + " |")
    print("\n**INDEL 0 / −1 按类型和长度**\n")
    subs = ("DEL1", "DEL>1", "INS1", "INS>1")
    print("| 标签 | " + " | ".join(f"{s} {v}" for s in subs for v in V) + " |"); print("|---|" + "---:|" * (len(subs) * len(V)))
    for lab in (0, -1):
        cells = []
        for s in subs:
            for v in V:
                cells.append(f"{sum(n for (k, l, r, ss), n in R[v]['glob'].items() if k == 'INDEL' and l == lab and ss == s):,}")
        print(f"| {lab} | " + " | ".join(cells) + " |")
table("I2 部分在图里 + 残余 edit", lambda t: t["final_class"].startswith("I2"))
table("I3 走 GRCh38，edits 写法不同", lambda t: t["final_class"].startswith("I3"))
table("I1 完全在图里（没有 edits）", lambda t: t["final_class"].startswith("I1"))
table("I4–I7", lambda t: t["final_class"][:2] in ("I4", "I5", "I6", "I7"))
table("全部 v5 no-candidate INDEL", lambda t: True)
tensor_table("I2 + I3", lambda t: t["final_class"][:2] in ("I2", "I3"))
global_table()

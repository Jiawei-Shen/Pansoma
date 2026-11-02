#!/usr/bin/env python3
"""Check which in-region rule rtg vcfeval -e applied in eval_bed_nograph_all: a truth record is kept iff the first base
of its anchor-trimmed allele (SNV: the base; DEL: first deleted base; INS: the base after the anchor) lies inside
eval_chr1.bed. Also lists the graph-carried SNV that stayed in the eval-BED denominator."""
import bisect, gzip, json
from pathlib import Path
R = Path("/scratch/jshen/data/pansoma_net_v2_runs"); X = Path("/scratch/jshen/data/reviewer_experiments_20261013/E23_nograph")
G = Path("/scratch/jshen/data/reviewer_experiments_20261013/E23_graph_alleles.mgi_non_HG005.d9.atomized.tsv.gz")
RUNS = {"SNV": "HapMap_Illumina_SNV_full", "INDEL": "HapMap_Illumina_INDEL"}
bed = [tuple(map(int, l.split("\t")[1:3])) for l in open(X / "eval_chr1.bed")]; bs = [b[0] for b in bed]
def inb(x):
    i = bisect.bisect_right(bs, x) - 1
    return i >= 0 and bed[i][0] <= x < bed[i][1]
def recs(d):
    out = {}
    for n in ("tp-baseline.vcf.gz", "fn.vcf.gz"):
        for line in gzip.open(d / n, "rt"):
            if line[0] == "#": continue
            c = line.rstrip("\n").split("\t"); v = dict(zip(c[8].split(":"), c[9].split(":")))
            af = v.get("AF", "."); out[(c[0], int(c[1]), c[3], c[4])] = None if af in (".", "") else float(af)
    return out
graph = {tuple([l.split("\t")[0], int(l.split("\t")[1])] + l.rstrip("\n").split("\t")[2:]) for l in gzip.open(G, "rt") if l.startswith("chr1\t")}
masked = {tuple([l.split("\t")[0], int(l.split("\t")[1])] + l.rstrip("\n").split("\t")[2:]) for l in gzip.open(X / "graph_carried_truth.tsv.gz", "rt") if l.startswith("chr1\t")}
res = {}
for kind, run in RUNS.items():
    F = {k for k, a in recs(R / run / "vcf_chr1/eval_bed/vcfeval/raw_calls").items() if a is not None and 0.0099 <= a <= 0.1651}
    E = {k for k, a in recs(R / run / "vcf_chr1/eval_bed_nograph_all/vcfeval/raw_calls").items() if a is not None and 0.0099 <= a <= 0.1651}
    agree = dis = 0; bad = []
    for k in F:
        _, p, r, a = k
        n = 0
        while n < min(len(r), len(a)) and r[n] == a[n]: n += 1  # shared leading bases (the anchor)
        first = p - 1 + n  # 0-based first base after the shared prefix (INS: the base after the insertion point)
        pred = inb(first)
        if pred == (k in E): agree += 1
        else: dis += 1; bad.append(k)
    stay = sorted(k for k in E if k in graph)
    res[kind] = dict(F=len(F), E=len(E), rule_agrees=agree, rule_disagrees=dis, disagree_examples=[list(b) for b in bad[:10]],
                     graph_carried_kept=[list(k) + [k in masked] for k in stay])
    print(kind, json.dumps(res[kind]), flush=True)
Path(__file__).with_suffix(".json").write_text(json.dumps(res, indent=1) + "\n")

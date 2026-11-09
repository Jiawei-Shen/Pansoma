# usage: python count.py SET   -- SNV candidates filtered by the AF floor (0.08) only, by AF, per chromosome, and the
# truth alleles (somatic / germline truth keys of truth/*.graph.tsv) among them: what an AF floor of 0.05 would add.
# Prefilter records carry af_upper_bound (ALT support upper bound / coverage); counted records carry af.
import bisect, json, sys
from collections import Counter, defaultdict
S = sys.argv[1]; D = f"/scratch/jshen/data/pansoma_v2_tensors/{S}"
R = [l.rstrip("\n").split("\t") for l in open("/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/hprc-v2.1-mc-grch38.d46.chr_node_ranges.tsv")][1:]
R = sorted((int(a), int(b), c) for c, a, b, *_ in R); starts = [r[0] for r in R]
def chrom(node):
    i = bisect.bisect_right(starts, node) - 1
    return R[i][2] if i >= 0 and node <= R[i][1] else "other"
truth = {}
for kind in ("somatic", "germline"):
    f = open(f"{D}/truth/{kind}.graph.tsv"); h = f.readline().rstrip("\n").split("\t")
    for l in f:
        x = dict(zip(h, l.rstrip("\n").split("\t")))
        if x["kind"] != "SNP" or x["passed"] != "True": continue
        for k in x["keys"].split(","):
            if k: truth.setdefault(k, (kind, x["in_bed"] == "True", x["chrom"]))
BINS = [0.0, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08]
lab = lambda af: next(f"[{BINS[i]:.2f}, {BINS[i+1]:.2f})" for i in range(len(BINS) - 1) if af < BINS[i + 1])
n = Counter(); sites = defaultdict(set); tr = Counter(); kinds = Counter()
for line in open(f"{D}/tensors/SNV/filtered_candidates.ndjson", "rb"):
    if b'"reasons": ["min_af"]' not in line: continue
    r = json.loads(line); af = r.get("af", r.get("af_upper_bound")); kinds["af" if "af" in r else "af_upper_bound"] += 1
    b = lab(af); c = chrom(r["node_id"]); g = "chr1" if c == "chr1" else "other chromosomes"
    n[(g, b)] += 1; sites[(g, b)].add(r["site_id"])
    t = truth.get(r["candidate_id"])
    if t: tr[(t[0], t[2] == "chr1", t[1], b)] += 1
print(S, "min_af-only records:", dict(kinds))
print("\n| AF (upper bound) | chr1 alleles | chr1 sites | other chromosomes alleles | sites |\n|---|---|---|---|---|")
for i in range(len(BINS) - 1):
    b = f"[{BINS[i]:.2f}, {BINS[i+1]:.2f})"
    print(f"| {b} | {n[('chr1', b)]:,} | {len(sites[('chr1', b)]):,} | {n[('other chromosomes', b)]:,} | {len(sites[('other chromosomes', b)]):,} |")
print("\ntruth alleles among them (kind, chr1, in BED, bin): count")
for k, v in sorted(tr.items()): print(" ", k, v)
cur = Counter()
for c in [f"chr{i}" for i in range(1, 23)]:
    try: cur["chr1" if c == "chr1" else "other"] += sum(1 for _ in open(f"{D}/tensors/SNV/{c}_variant_summary.ndjson", "rb"))
    except FileNotFoundError: pass
print("\ncurrent SNV tensors (sites): chr1", cur["chr1"], "other chromosomes", cur["other"])

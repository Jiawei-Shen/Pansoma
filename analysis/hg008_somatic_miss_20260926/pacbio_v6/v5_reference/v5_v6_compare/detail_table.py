"""Truth-level and tensor-level comparison of v5 no-candidate somatic INDELs across label sets.

Related tensors of a truth = (a) GRCh38-placed tensors within [span_lo-10, span_hi+10], where the span is the
truth's equivalence span (right-shifted through the repeat), plus (b) tensors on the truth's residual-edit nodes
(from the read analysis) within ±15 bp of the residual edit offset, including tensors on non-GRCh38 nodes.
Usage: python3 detail_table.py v5 v6 v6new
"""
import bisect, collections, csv, glob, json, pickle, re, sys
csv.field_size_limit(sys.maxsize)
ROOT = "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125"
DIRS = dict(v5=f"{ROOT}/v5_tensors", v6=f"{ROOT}/v6_tensors",
            v6new="/scratch/jshen/Github/Pansoma/tmp/v5_v6_compare/v6_relabel")
MISS = "/scratch/jshen/Github/Pansoma/tmp/somatic_miss_analysis_20260923/indel_no_candidate.tsv"
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
RESID = re.compile(r"(\d+):(\d+):(SNP|INS|DEL):")

fai = {l.split("\t")[0]: tuple(map(int, l.split("\t")[1:5])) for l in open(FASTA + ".fai")}
fa = open(FASTA, "rb")
def seq(chrom, start, end):
    length, offset, bases, width = fai[chrom]
    end = min(end, length); out = []
    fa.seek(offset + start // bases * width + start % bases)
    raw = fa.read((end - start) + (end - start) // bases * (width - bases) + width)
    return raw.replace(b"\n", b"").replace(b"\r", b"")[:end - start].decode().upper()

def span(r):
    """0-based [lo, hi) of the truth's equivalence span (VCF is left-aligned; shift right through the repeat)."""
    p0, ref, alt = int(r["vcf_pos"]) - 1, r["vcf_ref"].upper(), r["vcf_alt"].upper()
    chrom = r["chrom"]; ctx = seq(chrom, p0, p0 + len(ref) + 5000)
    if len(ref) > len(alt):                      # deletion of ref[1:]
        d = ref[len(alt):]; s = 1; e = s + len(d)
        while e < len(ctx) and ctx[e] == ctx[s]: s += 1; e += 1
        return p0 + 1, p0 + e
    ins = alt[len(ref):]; i = 0                  # insertion after ref
    while len(ref) + i < len(ctx) and ctx[len(ref) + i] == ins[i % len(ins)]: i += 1
    return p0 + len(ref), p0 + len(ref) + i

def parse_cand(cid):
    node, off, kind, allele = cid.split(":", 3)
    ref, alt = allele.split("@")[0].split(">")
    return int(node), int(off), kind, len(alt) - len(ref)

truths = [r for r in csv.DictReader(open(MISS), delimiter="\t")]
want = collections.defaultdict(list)             # node -> [(truth index, residual offset)]
for ti, r in enumerate(truths):
    r["span"] = span(r)
    for m in RESID.finditer(r["residual_edits(reads)=tensor_label"]):
        want[int(m.group(1))].append((ti, int(m.group(2))))

def load(version):
    bypos, bynode, glob_counts = collections.defaultdict(list), collections.defaultdict(list), collections.Counter()
    for tk in ("SNV", "INDEL"):
        for f in glob.glob(f"{DIRS[version]}/{tk}/chr*_labels.ndjson"):
            for line in open(f):
                l = json.loads(line); node, off, kind, dlen = parse_cand(l["candidate_id"])
                key = (l["label"], l["reason"])
                sub = "SNV" if tk == "SNV" else f"{kind}{'1' if abs(dlen) == 1 else '>1'}"
                glob_counts[(tk, l["label"], l["reason"], sub)] += 1
                g = l.get("grch38")
                item = (l["candidate_id"], tk, l["label"], l["reason"], dlen)
                if g and g.get("pos0") is not None: bypos[g["chrom"]].append((g["pos0"],) + item)
                if node in want: bynode[node].append((off,) + item + (g is not None,))
    for c in bypos: bypos[c].sort()
    return bypos, bynode, glob_counts


result = {}
for version in sys.argv[1:]:
    bypos, bynode, glob_counts = load(version)
    keys = {c: [x[0] for x in v] for c, v in bypos.items()}
    per_truth = []
    for ti, r in enumerate(truths):
        lo, hi = r["span"]; ks = keys.get(r["chrom"], [])
        hits = {h[1]: (h[3], h[4], h[2], h[5], False) for h in bypos[r["chrom"]][bisect.bisect_left(ks, lo - 10):bisect.bisect_right(ks, hi + 10)]}
        for m in RESID.finditer(r["residual_edits(reads)=tensor_label"]):
            node, off = int(m.group(1)), int(m.group(2))
            for h in bynode.get(node, []):
                if abs(h[0] - off) <= 15:
                    hits.setdefault(h[1], (h[3], h[4], h[2], h[5], not h[6]))
        per_truth.append(list(hits.values()))
    result[version] = dict(per_truth=per_truth, glob=glob_counts)
    print(version, "loaded", flush=True)
    del bypos, bynode, keys
    pickle.dump(dict(truths=[{k: r[k] for k in ("final_class", "chrom", "vcf_pos", "vcf_ref", "vcf_alt", "kind", "reason", "span")}
                             for r in truths], result=result[version]), open(f"detail_{version}.pkl", "wb"))

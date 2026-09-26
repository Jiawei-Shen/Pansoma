"""Overlap of the somatic (v0.2) and germline (dipcall) truth sets, as the labelling code sees them."""
import csv, sys, json, glob, collections, bisect
csv.field_size_limit(sys.maxsize)
T = "/scratch/jshen/data/pansoma_v2_tensors/truth"
som = list(csv.DictReader(open(f"{T}/somatic.graph.tsv"), delimiter="\t"))
skey = collections.defaultdict(list)
for r in som:
    for k in (r["keys"].split(",") if r["keys"] else []):
        skey[k].append(r["truth_id"])
spos = {(r["chrom"], r["pos0"]) for r in som}
sall = {(r["chrom"], r["pos0"], r["ref"], r["alt"]): r for r in som}
svcf = {(r["chrom"], r["vcf_pos"], r["vcf_ref"], r["vcf_alt"]) for r in som}
same_allele, same_pos, key_shared, vcf_same = {}, set(), collections.defaultdict(set), set()
germ_pos = collections.defaultdict(list)
n = 0
with open(f"{T}/germline.graph.tsv") as f:
    cols = next(f).rstrip("\n").split("\t")
    for line in f:
        x = dict(zip(cols, line.rstrip("\n").split("\t")))
        n += 1
        a = (x["chrom"], x["pos0"], x["ref"], x["alt"])
        if a in sall: same_allele[sall[a]["truth_id"]] = x
        if (x["chrom"], x["pos0"]) in spos: same_pos.add((x["chrom"], x["pos0"]))
        if (x["chrom"], x["vcf_pos"], x["vcf_ref"], x["vcf_alt"]) in svcf: vcf_same.add((x["chrom"], x["vcf_pos"]))
        for k in (x["keys"].split(",") if x["keys"] else []):
            if k in skey:
                for t in skey[k]: key_shared[t].add((x["truth_id"], x["passed"], x["in_bed"], x["gt"]))
        germ_pos[x["chrom"]].append(int(x["pos0"]))
for c in germ_pos: germ_pos[c].sort()
kind = {r["truth_id"]: r["kind"] for r in som}
pib = {r["truth_id"]: r["passed"] == "True" and r["in_bed"] == "True" for r in som}
def by_kind(ids): 
    c = collections.Counter(kind[t] for t in ids); return dict(c), sum(c.values())
print("somatic alleles", len(som), collections.Counter(r["kind"] for r in som), "germline alleles", n)
print("same allele (normalized chrom,pos0,ref,alt):", by_kind(same_allele), "of which somatic PASS&BED:", sum(pib[t] for t in same_allele))
print("  germline side of those: PASS&BED", sum(v["passed"] == "True" and v["in_bed"] == "True" for v in same_allele.values()), collections.Counter(v["gt"] for v in same_allele.values()).most_common(5))
print("same raw VCF record (chrom,pos,ref,alt):", len(vcf_same))
print("somatic whose position also has a germline allele:", sum((r["chrom"], r["pos0"]) in same_pos for r in som))
print("somatic sharing >=1 node key with a germline allele:", by_kind(key_shared))
good = [t for t, v in key_shared.items() if any(p == "True" and b == "True" for _, p, b, _ in v)]
print("   ... with a PASS & in-BED germline allele:", by_kind(good))
near = 0
for r in som:
    a = germ_pos.get(r["chrom"], []); p = int(r["pos0"])
    near += bisect.bisect_right(a, p + 10) > bisect.bisect_left(a, p - 10)
print("somatic with any germline allele within 10 bp:", near)
json.dump(sorted(key_shared), open("/tmp/claude-10008/-scratch-jshen-Github-Pansoma/c9020fc2-7f52-41b9-b09f-9e3636ed53ce/scratchpad/key_shared.json", "w"))
# tensors carrying both a somatic and a germline match
for name, d in (("PacBio v6", "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/v6_tensors"),
                ("ONT", "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_Northeastern-ONT-UL-20241216/v3_tensors")):
    c = collections.Counter()
    for f in glob.glob(f"{d}/*/chr*_labels.ndjson"):
        for line in open(f):
            if '"somatic": []' in line or '"germline": []' in line: continue
            l = json.loads(line)
            rs = any(s["representative"] for s in l["somatic"]); rg = any(g["representative"] for g in l["germline"])
            c[(f.split("/")[-2], l["label_name"], "rep_somatic" if rs else "nonrep_somatic", "rep_germline" if rg else "nonrep_germline")] += 1
    print(name, "tensors matching both a somatic and a germline truth:", sum(c.values()))
    for k, v in sorted(c.items()): print("   ", k, v)

"""Per-locus outcome of altlen_tumor_asm_22.tsv, the R3.7 rule: each tumor haplotype is 'truth' if any hit equals
GRCh38 + truth, else 'normal' if any hit equals a COLO829BL haplotype sequence, else 'novel' (a second hit is a paralog
copy); then truth on any tumor hap -> somatic, all haps normal -> germline-like, else unresolved. Checked against R3.7's
summary for the 18 shared loci. usage: python summarize.py (writes altlen_tumor_asm_22_summary.tsv)"""
import csv, collections
rows = list(csv.DictReader(open("altlen_tumor_asm_22.tsv"), delimiter="\t"))
by = collections.defaultdict(lambda: collections.defaultdict(list))
first = {}
for r in rows:
    first.setdefault(r["truth_id"], r)
    if r["tumor_hap"] != "none" and r["len_T"] != "":
        by[r["truth_id"]][r["tumor_hap"]].append(r)
def hap_state(hits):
    if any(h["T_eq_A"] == "True" for h in hits): return "truth"
    if any(h["T_eq_hapX"] == "True" or h["T_eq_hapY"] == "True" for h in hits): return "normal"
    return "novel"
def outcome(t):
    st = [hap_state(v) for v in by[t].values()]
    if not st: return "no hit"
    if "truth" in st: return "tumor carries the truth allele (not in normal)"
    if all(s == "normal" for s in st): return "tumor carries only normal haplotype sequences"
    if "normal" in st: return "one tumor hap = normal hap, other hap novel non-truth sequence"
    return "novel non-truth sequence only"
old = {r["truth_id"]: r["outcome"] for r in csv.DictReader(open("/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/misc/r37_tumor_asm_35_summary.tsv"), delimiter="\t")}
c, agree = collections.Counter(), collections.Counter()
with open("altlen_tumor_asm_22_summary.tsv", "w") as w:
    w.write("truth_id\tchrom\tvcf_pos\tvcf_ref\tvcf_alt\talt_len_haps\tpf_class\toutcome\tr37_outcome\n")
    for t in [x.strip() for x in open("ids_22.txt") if x.strip()]:
        o, h = outcome(t), first[t]
        c[o] += 1
        if t in old:
            agree[o == old[t]] += 1
        w.write("\t".join([t, h["chrom"], h["vcf_pos"], h["vcf_ref"], h["vcf_alt"], h["alt_len_haps"], h["pf_class"], o, old.get(t, "")]) + "\n")
print(dict(c)); print("R3.7 loci reproduced:", dict(agree))

"""Pick example loci per miss class (PASS, in BED, most error-free ALT reads) and write show_loci.txt."""
import csv, sys
from collections import defaultdict
csv.field_size_limit(sys.maxsize)
rows = [r for r in csv.DictReader(open("indel_no_candidate.tsv"), delimiter="\t") if r["passed"] == "True" and r["in_bed"] == "True"]
by = defaultdict(list)
for r in rows:
    by[(r["final_class"], r["kind"])].append(r)
with open("show_loci.txt", "w") as out, open("examples_pick.tsv", "w") as tab:
    tab.write("final_class\tkind\tchrom\tvcf_pos\tvcf_ref\tvcf_alt\talt_exact\talt_like\tspanning\treason\tresidual\n")
    for (cl, kind), items in sorted(by.items()):
        items.sort(key=lambda r: (-int(r["alt_exact"]), -int(r["alt_like"])))
        for r in items[:2]:
            out.write(f"{r['chrom']} {r['vcf_pos']}\n")
            tab.write("\t".join([cl, kind, r["chrom"], r["vcf_pos"], r["vcf_ref"][:30], r["vcf_alt"][:30], r["alt_exact"],
                                 r["alt_like"], r["spanning"], r["reason"], r["residual_edits(reads)=tensor_label"][:200]]) + "\n")
print(open("examples_pick.tsv").read())

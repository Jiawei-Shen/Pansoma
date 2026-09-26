"""(ONT v3) Somatic SNVs with vs without a candidate: GRCh38 node length at the site, germline variants at / near the site."""
import csv, sys, collections
import numpy as np
sys.path.insert(0, "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_Northeastern-ONT-UL-20241216/v3_run/source")  # the run's frozen package (indexed_gam_pipeline_v2 was removed from the repo)
from indexed_gam_pipeline_v3.tensor_postprocessing.reference_path import ReferencePath
from indexed_gam_pipeline_v3.tensor_postprocessing.truth_labels import Locator
csv.field_size_limit(sys.maxsize)
R = "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_Northeastern-ONT-UL-20241216/v3_tensors"
T = "/scratch/jshen/data/pansoma_v2_tensors/truth"
status = {r["truth_id"]: r["status"] for r in csv.DictReader(open(f"{R}/somatic.recall.tsv"), delimiter="\t")}
som = [r for r in csv.DictReader(open(f"{T}/somatic.graph.tsv"), delimiter="\t") if r["kind"] == "SNP"]
germ = collections.defaultdict(list)
germ_alt = {}
with open(f"{T}/germline.graph.tsv") as f:
    next(f)
    for line in f:
        x = line.split("\t", 10)
        chrom, pos0, ref, alt, kind = x[1], int(x[5]), x[6], x[7], x[8]
        germ[chrom].append(pos0)
        if kind == "SNP":
            germ_alt[(chrom, pos0)] = germ_alt.get((chrom, pos0), "") + alt
germ = {c: np.array(sorted(v)) for c, v in germ.items()}
loc = Locator(ReferencePath("/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path"))
tab = collections.defaultdict(collections.Counter)
for r in sorted(som, key=lambda r: (r["chrom"], int(r["pos0"]))):
    st = status[r["truth_id"]]
    st = "hit" if st == "tensor_representative" else "no_candidate" if st == "no_candidate" else None
    if st is None:
        continue
    p = int(r["pos0"]); c = r["chrom"]
    hit = loc.at(c, p)
    L = hit[2] if hit else None
    tab[st]["node 1bp" if L == 1 else "node 2-32bp" if L and L <= 32 else "node >32bp" if L else "no unique node"] += 1
    a = germ.get(c, np.array([]))
    n50 = int(np.searchsorted(a, p + 50, "right") - np.searchsorted(a, p - 50, "left"))
    tab[st]["germline at same pos"] += (c, p) in germ_alt
    tab[st]["germline same pos, same ALT"] += r["alt"] in germ_alt.get((c, p), "")
    tab[st]["germline within 50bp: 0" if n50 == 0 else "germline within 50bp: 1" if n50 == 1 else "germline within 50bp: >=2"] += 1
    tab[st]["total"] += 1
for st, c in tab.items():
    tot = c["total"]
    print(st, tot)
    for k in sorted(c):
        if k != "total":
            print(f"   {k:32s} {c[k]:5d}  {100 * c[k] / tot:5.1f}%")

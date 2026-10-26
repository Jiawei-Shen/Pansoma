#!/usr/bin/env python3
"""chr1 HapMap truth alleles (benchmark BED) by the number of minor donors that carry them (INFO/<donor>_GT), and the
share that are paths of the graph (graph_alleles.sh output): supports 'the in-graph positives are mostly shared alleles'.
usage: python donors_in_graph.py GRAPH_ALLELES.tsv.gz [GRAPH_ALLELES_ATOMIZED.tsv.gz]
"""
import gzip
import subprocess
import sys
from collections import Counter

TRUTH = "/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/HapMap_somatic_snv_indel.vcf.gz"
BED = "/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/SMHTHAPMAP6_GRCh38_v1.6_somatic_benchmark_autosomes.bed"
DONORS = ["HG002", "HG00438", "HG02257", "HG02486", "HG02622"]


def graph(path):
    out = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("chr1\t"):
                c, pos, ref, alt = line.rstrip("\n").split("\t")
                out.add((int(pos), ref, alt))
    return out


graphs = [graph(p) for p in sys.argv[1:]]
fmt = "%POS\t%REF\t%ALT\t" + "\t".join(f"%INFO/{d}_GT" for d in DONORS) + "\n"
q = subprocess.run(f"bcftools view -r chr1 -T {BED} {TRUTH} -Ou | bcftools query -f '{fmt}'", shell=True,
                   capture_output=True, text=True, check=True).stdout
n = Counter()
hit = [Counter() for _ in graphs]
for line in q.splitlines():
    pos, ref, alt, *gts = line.split("\t")
    kind = "SNV" if len(ref) == 1 and len(alt) == 1 else "INDEL"
    carriers = "unknown" if any("*" in g or g == "." for g in gts) else str(sum("1" in g for g in gts))
    key = (int(pos), ref, alt)
    n[kind, carriers] += 1
    for i, g in enumerate(graphs):
        hit[i][kind, carriers] += key in g
print("kind\tdonors carrying\ttruth\t" + "\t".join(f"in graph ({p.split('/')[-1]})" for p in sys.argv[1:]))
for k in sorted(n):
    print(f"{k[0]}\t{k[1]}\t{n[k]}\t" + "\t".join(f"{h[k] / n[k]:.3f}" for h in hit))

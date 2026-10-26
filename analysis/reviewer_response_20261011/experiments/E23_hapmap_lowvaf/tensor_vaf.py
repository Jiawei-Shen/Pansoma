#!/usr/bin/env python3
"""E23 tensor level: what happens to each HapMap truth allele, by its expected VAF.

Inputs (HapMap_Illumina_nonHG005, labels of label_run / label_parallel.py):
  <set>/truth/somatic.graph.tsv     truth_id -> chrom, vcf_pos, vcf_ref, vcf_alt, kind, in_bed
  <set>/tensors/somatic.recall.tsv  truth_id -> status (tensor_representative, tensor_non_representative_allele,
                                    filtered + reasons, no_candidate, complex_allele, no_unique_grch38_node_key)
  the somatic truth VCF             FORMAT AF = the allele's expected VAF in the mixture (missing for some records)
  --graph-alleles (optional)        CHROM POS REF ALT of the graph's alleles (graph_alleles.sh): a truth allele found
                                    there is a path of the graph (absorbed: reads that carry it align without an edit)
Output (--output DIR): tensor_vaf.tsv (kind, scope, VAF bin, in_graph, status, n) and tables.md (shares per bin).
Scopes: all autosomes and chr1 (the test chromosome), both inside the benchmark BED only.
usage: python tensor_vaf.py --set /scratch/jshen/data/pansoma_v2_tensors/HapMap_Illumina_nonHG005 --vcf TRUTH.vcf.gz
       [--graph-alleles GRAPH.tsv.gz] --output DIR
"""
import argparse
import gzip
from collections import Counter, defaultdict
from pathlib import Path

# the truth VCF's AF values are multiples of 0.0025 (0.0025-0.165); bins keep the low end fine
BINS = [(0, 0.004, "0.0025"), (0.004, 0.0075, "0.005"), (0.0075, 0.011, "0.01"), (0.011, 0.016, "0.0125-0.015"),
        (0.016, 0.026, "0.02-0.025"), (0.026, 0.049, "0.03-0.045"), (0.049, 0.079, "0.05-0.075"),
        (0.079, 0.099, "0.08-0.095"), (0.099, 1.0, "0.10-0.165")]
STATES = ["tensor_representative", "tensor_non_representative_allele", "filtered:min_af", "filtered:min_variants",
          "filtered:other", "no_candidate", "other"]


def vaf_bin(af):
    if af is None:
        return "missing"
    for lo, hi, name in BINS:
        if lo <= af < hi:
            return name
    return "other"


def read_vcf_af(path):
    """(chrom, pos, ref, alt) -> AF (float, or None for '.')."""
    out = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t", 10)
            keys, values = c[8].split(":"), c[9].split(":")
            af = dict(zip(keys, values)).get("AF", ".")
            for alt in c[4].split(","):
                out[(c[0], int(c[1]), c[3], alt)] = None if af in (".", "") else float(af)
    return out


def state_of(status, detail):
    if status.startswith("tensor_"):
        return status
    if status == "filtered":
        reasons = set(detail.split(","))
        if "min_af" in reasons:
            return "filtered:min_af"
        if "min_variants" in reasons:
            return "filtered:min_variants"
        return "filtered:other"
    return status if status in STATES else "other"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--set", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--graph-alleles")
    p.add_argument("--output", required=True)
    a = p.parse_args()
    root, out = Path(a.set), Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    af = read_vcf_af(a.vcf)
    graph = None
    if a.graph_alleles:
        graph = set()
        with gzip.open(a.graph_alleles, "rt") as f:
            for line in f:
                c, pos, ref, alt = line.rstrip("\n").split("\t")
                graph.add((c, int(pos), ref, alt))
    table = {}  # truth_id -> (chrom, vcf_pos, vcf_ref, vcf_alt, kind, in_bed); the keys column can be huge: split
    with open(root / "truth" / "somatic.graph.tsv") as f:
        head = f.readline().rstrip("\n").split("\t")
        col = {k: i for i, k in enumerate(head)}
        last = max(col[k] for k in ("truth_id", "chrom", "vcf_pos", "vcf_ref", "vcf_alt", "kind", "in_bed")) + 1
        for line in f:
            r = line.split("\t", last)
            table[r[col["truth_id"]]] = (r[col["chrom"]], int(r[col["vcf_pos"]]), r[col["vcf_ref"]], r[col["vcf_alt"]],
                                         r[col["kind"]], r[col["in_bed"]])
    counts = defaultdict(Counter)  # (kind, scope, bin, in_graph) -> state counts
    with open(root / "tensors" / "somatic.recall.tsv") as f:
        head = f.readline().rstrip("\n").split("\t")
        col = {k: i for i, k in enumerate(head)}
        for line in f:
            r = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, kind, in_bed = table[r[col["truth_id"]]]
            if in_bed != "True" or r[col["passed"]] != "True":
                continue
            kind = "SNV" if kind == "SNP" else "INDEL" if kind in ("INS", "DEL") else kind
            b = vaf_bin(af.get((chrom, pos, ref, alt)))
            s = state_of(r[col["status"]], r[col["detail"]])
            g = "na" if graph is None else ("yes" if (chrom, pos, ref, alt) in graph else "no")
            for scope in ("autosomes", "chr1") if chrom == "chr1" else ("autosomes",):
                counts[kind, scope, b, g][s] += 1
    with open(out / "tensor_vaf.tsv", "w") as f:
        f.write("kind\tscope\tvaf_bin\tin_graph\tstate\tn\n")
        for (kind, scope, b, g), c in sorted(counts.items()):
            for s, n in c.items():
                f.write(f"{kind}\t{scope}\t{b}\t{g}\t{s}\t{n}\n")
    merged = defaultdict(Counter)  # over in_graph, for the main table
    for (kind, scope, b, g), c in counts.items():
        merged[kind, scope, b].update(c)
    order = [name for _, _, name in BINS] + ["missing", "other"]
    lines = ["# E23 tensor level: HapMap truth alleles by expected VAF (benchmark BED, PASS truth)", ""]
    for kind in ("SNV", "INDEL"):
        for scope in ("autosomes", "chr1"):
            rows = [(b, merged[kind, scope, b]) for b in order if merged.get((kind, scope, b))]
            if not rows:
                continue
            lines += [f"## {kind}, {scope}", "",
                      "| VAF | truth | with a tensor | filtered: AF floor | filtered: < 3 ALT reads | other filter | "
                      "no candidate |", "|---|---|---|---|---|---|---|"]
            for b, c in rows:
                n = sum(c.values())
                t = c["tensor_representative"] + c["tensor_non_representative_allele"]
                lines.append(f"| {b} | {n:,} | {t / n:.3f} | {c['filtered:min_af'] / n:.3f} | "
                             f"{c['filtered:min_variants'] / n:.3f} | {c['filtered:other'] / n:.3f} | "
                             f"{c['no_candidate'] / n:.3f} |")
            lines.append("")
            if graph is not None:
                lines += [f"### {kind}, {scope}: by whether the allele is a path of the graph", "",
                          "| VAF | in graph | truth | with a tensor | filtered | no candidate |", "|---|---|---|---|---|---|"]
                for b in order:
                    for g in ("yes", "no"):
                        c = counts.get((kind, scope, b, g))
                        if not c:
                            continue
                        n = sum(c.values())
                        t = c["tensor_representative"] + c["tensor_non_representative_allele"]
                        fl = sum(v for k, v in c.items() if k.startswith("filtered"))
                        lines.append(f"| {b} | {g} | {n:,} | {t / n:.3f} | {fl / n:.3f} | {c['no_candidate'] / n:.3f} |")
                lines.append("")
    (out / "tables.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

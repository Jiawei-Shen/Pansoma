#!/usr/bin/env python3
"""E23 without graph-carried positives (user 2026-10-29): the HapMap truth alleles that are paths of the graph
(the atomized graph alleles of graph_alleles.sh, the definition E23 reports) are taken out of training and evaluation.

Writes into OUT:
  graph_carried_truth.tsv.gz        chrom, vcf_pos, vcf_ref, vcf_alt of every graph-carried truth allele (autosomes)
  exclude.<KIND>.tsv                candidate ids whose somatic truth (exact or partial match) is graph-carried, for
                                    train / predict --exclude-candidates (left out of training, validation and the test)
  eval_chr1.bed                     the benchmark BED on chr1 minus every graph-carried truth allele's span (SNVs and
                                    INDELs, +-1 bp): vcfeval inside it leaves these sites out of recall and precision
  summary.tsv                       per kind and chromosome: tensors with a somatic truth, excluded
usage: python exclude_graph_carried.py --set SET_DIR --graph-alleles G.tsv.gz --out OUT [--processes N]
"""
import argparse
import gzip
import json
from multiprocessing import Pool
from pathlib import Path

GRAPH = set()  # (chrom, vcf_pos, vcf_ref, vcf_alt) of the graph-carried truth alleles, shared with the workers by fork


def scan(path):
    """(kind, chrom, tensors with a somatic truth, excluded ids) of one <chrom>_labels.ndjson."""
    kind, chrom = path.parent.name, path.name.split("_")[0]
    n, out = 0, []
    with open(path) as f:
        for line in f:
            if '"truth_id"' not in line:
                continue
            r = json.loads(line)
            truths = list(r.get("somatic") or [])
            if r.get("partial_truth"):
                truths.append(r["partial_truth"])
            if not truths:
                continue
            n += 1
            if any((t["chrom"], int(t["vcf_pos"]), t["vcf_ref"], t["vcf_alt"]) in GRAPH for t in truths):
                out.append(r["candidate_id"])
    return kind, chrom, n, out


def subtract(regions, holes):
    """Sorted, merged [start, end) regions minus sorted holes."""
    out, j = [], 0
    for s, e in regions:
        while j < len(holes) and holes[j][1] <= s:
            j += 1
        k, cur = j, s
        while k < len(holes) and holes[k][0] < e:
            hs, he = holes[k]
            if hs > cur:
                out.append((cur, hs))
            cur = max(cur, he)
            k += 1
        if cur < e:
            out.append((cur, e))
    return out


def merged(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--set", required=True)
    p.add_argument("--graph-alleles", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--processes", type=int, default=16)
    a = p.parse_args()
    root, out = Path(a.set), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    truth = set()
    with open(root / "truth" / "somatic.graph.tsv") as f:
        col = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
        last = max(col[k] for k in ("chrom", "vcf_pos", "vcf_ref", "vcf_alt")) + 1
        for line in f:
            r = line.split("\t", last)
            truth.add((r[col["chrom"]], int(r[col["vcf_pos"]]), r[col["vcf_ref"]], r[col["vcf_alt"]]))
    with gzip.open(a.graph_alleles, "rt") as f:
        for line in f:
            c, pos, ref, alt = line.rstrip("\n").split("\t")
            key = (c, int(pos), ref, alt)
            if key in truth:
                GRAPH.add(key)
    del truth
    with gzip.open(out / "graph_carried_truth.tsv.gz", "wt") as o:
        o.write("".join(f"{c}\t{pos}\t{ref}\t{alt}\n" for c, pos, ref, alt in sorted(GRAPH)))
    print(f"graph-carried truth alleles: {len(GRAPH):,}", flush=True)
    files = sorted((root / "tensors").glob("*/*_labels.ndjson"))
    summary = []
    with Pool(a.processes) as pool:
        by_kind = {}
        for kind, chrom, n, ids in pool.imap_unordered(scan, files):
            by_kind.setdefault(kind, []).extend((i, chrom) for i in ids)
            summary.append((kind, chrom, n, len(ids)))
            print(f"{kind} {chrom}: {n:,} tensors with a somatic truth, {len(ids):,} graph-carried", flush=True)
    for kind, rows in by_kind.items():
        with open(out / f"exclude.{kind}.tsv", "w") as o:
            o.write("candidate_id\tchrom\n" + "".join(f"{i}\t{c}\n" for i, c in sorted(rows)))
    with open(out / "summary.tsv", "w") as o:
        o.write("kind\tchrom\tsomatic_truth_tensors\texcluded\n")
        o.write("".join(f"{k}\t{c}\t{n}\t{x}\n" for k, c, n, x in sorted(summary)))
    labels = json.load(open(root / "tensors" / "SNV" / "labels.manifest.json"))
    bed = merged((int(s), int(e)) for c, s, e in (l.split("\t")[:3] for l in open(labels["truth"]["somatic"]["bed"]))
                 if c == "chr1")
    holes = merged((pos - 2, pos - 1 + len(ref) + 1) for c, pos, ref, alt in GRAPH if c == "chr1")
    kept = subtract(bed, holes)
    with open(out / "eval_chr1.bed", "w") as o:
        o.write("".join(f"chr1\t{s}\t{e}\n" for s, e in kept))
    size = lambda iv: sum(e - s for s, e in iv)
    print(f"chr1 BED {size(bed):,} bp -> {size(kept):,} bp after removing {len(holes):,} graph-carried spans", flush=True)


if __name__ == "__main__":
    main()

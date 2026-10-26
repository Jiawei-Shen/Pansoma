#!/usr/bin/env python3
"""E23 call level: chr1 recall per truth VAF bin and precision per call AF bin, from a run's rtg vcfeval outputs.

<run>/vcf_chr1/eval_bed/vcfeval/{raw_calls,pon_calls}/ (written by machine_learning/pansoma_net_v2/jobs/vcf.sh):
  tp-baseline + fn   truth records in the benchmark BED, FORMAT AF = the expected VAF in the mixture
  tp + fp            PASS calls, FORMAT AF = the allele fraction in the reads
raw_calls = before population filtering (the detection limit E23 reports), pon_calls = after it.
--graph-alleles (graph_alleles.sh): truth recall is also split by whether the allele is a path of the graph.
usage: python call_vaf.py --run /scratch/jshen/data/pansoma_net_v2_runs/<RUN> --kind SNV|INDEL [--graph-alleles G.tsv.gz]
       --output DIR
"""
import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

from tensor_vaf import BINS, vaf_bin  # the same truth VAF bins as the tensor-level table

CALL_BINS = [(0, 0.02, "< 0.02"), (0.02, 0.03, "0.02-0.03"), (0.03, 0.05, "0.03-0.05"), (0.05, 0.08, "0.05-0.08"),
             (0.08, 0.2, "0.08-0.20"), (0.2, 1.01, "≥ 0.20")]


def records_af(path, with_key=False):
    """FORMAT AF of every record (first sample), None when missing; with_key: ((chrom, pos, ref, alt), AF)."""
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            values = dict(zip(c[8].split(":"), c[9].split(":")))
            af = values.get("AF", ".")
            af = None if af in (".", "") else float(af)
            out.append(((c[0], int(c[1]), c[3], c[4]), af) if with_key else af)
    return out


def call_bin(af):
    for lo, hi, name in CALL_BINS:
        if af is not None and lo <= af < hi:
            return name
    return "missing"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--kind", required=True)
    p.add_argument("--graph-alleles")
    p.add_argument("--output", required=True)
    a = p.parse_args()
    run, out = Path(a.run), Path(a.output)
    graph = None
    if a.graph_alleles:
        graph = set()
        with gzip.open(a.graph_alleles, "rt") as f:
            for line in f:
                c, pos, ref, alt = line.rstrip("\n").split("\t")
                graph.add((c, int(pos), ref, alt))
    out.mkdir(parents=True, exist_ok=True)
    result = {}
    lines = [f"# E23 call level: {run.name}, {a.kind}, chr1, benchmark BED", ""]
    for arm, title in (("raw_calls", "before population filtering"), ("pon_calls", "after population filtering")):
        d = run / "vcf_chr1" / "eval_bed" / "vcfeval" / arm
        tp, fn = Counter(map(vaf_bin, records_af(d / "tp-baseline.vcf.gz"))), Counter(map(vaf_bin, records_af(d / "fn.vcf.gz")))
        ctp, cfp = Counter(map(call_bin, records_af(d / "tp.vcf.gz"))), Counter(map(call_bin, records_af(d / "fp.vcf.gz")))
        result[arm] = dict(truth_tp=dict(tp), truth_fn=dict(fn), call_tp=dict(ctp), call_fp=dict(cfp))
        lines += [f"## {title}", "", "| truth VAF | truth | TP | recall |", "|---|---|---|---|"]
        for b in [name for _, _, name in BINS] + ["missing"]:
            n = tp[b] + fn[b]
            if n:
                lines.append(f"| {b} | {n:,} | {tp[b]:,} | {tp[b] / n:.3f} |")
        n, t = sum(tp.values()) + sum(fn.values()), sum(tp.values())
        lines += [f"| all | {n:,} | {t:,} | {t / n:.3f} |", "", "| call AF | calls | TP | precision |", "|---|---|---|---|"]
        for b in [name for _, _, name in CALL_BINS] + ["missing"]:
            n = ctp[b] + cfp[b]
            if n:
                lines.append(f"| {b} | {n:,} | {ctp[b]:,} | {ctp[b] / n:.3f} |")
        n, t = sum(ctp.values()) + sum(cfp.values()), sum(ctp.values())
        lines += [f"| all | {n:,} | {t:,} | {t / n:.3f} |", ""]
        if graph is not None:
            split = Counter()  # (bin, in graph, TP or FN)
            for name, path in (("TP", "tp-baseline.vcf.gz"), ("FN", "fn.vcf.gz")):
                for key, af in records_af(d / path, with_key=True):
                    split[vaf_bin(af), key in graph, name] += 1
            result[arm]["truth_by_graph"] = {f"{b}|{g}|{k}": v for (b, g, k), v in split.items()}
            lines += ["| truth VAF | allele in graph | truth | recall |", "|---|---|---|---|"]
            for b in [name for _, _, name in BINS] + ["missing"]:
                for g in (True, False):
                    n = split[b, g, "TP"] + split[b, g, "FN"]
                    if n:
                        lines.append(f"| {b} | {'yes' if g else 'no'} | {n:,} | {split[b, g, 'TP'] / n:.3f} |")
            lines.append("")
    (out / f"call_vaf.{run.name}.json").write_text(json.dumps(result, indent=1) + "\n")
    (out / f"call_vaf.{run.name}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

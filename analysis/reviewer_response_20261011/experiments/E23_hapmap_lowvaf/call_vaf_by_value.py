#!/usr/bin/env python3
"""E23 call level at every truth VAF value from 0.01 to 0.165 (the HapMap truth VAFs are discrete: sums of donor
fractions), chr1, benchmark BED, without population filtering (raw_calls of vcf.sh): recall of the truth not carried
by the graph and of all truth per value, and precision per 1 % step of the calls' own allele fraction.
usage: python call_vaf_by_value.py (--run RUN_DIR | --eval-dir VCFEVAL_DIR) --kind SNV|INDEL [--graph-alleles G.tsv.gz]
       [--title T] --output OUT.md   (--eval-dir: any rtg vcfeval output, e.g. Mutect2 calls)
"""
import argparse
import gzip
from collections import Counter
from pathlib import Path

from call_vaf import records_af


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run")
    p.add_argument("--eval-dir", help="an rtg vcfeval output directory instead of <run>/vcf_chr1/eval_bed/vcfeval/raw_calls")
    p.add_argument("--title")
    p.add_argument("--kind", required=True)
    p.add_argument("--graph-alleles", help="without it every truth allele counts as not graph-carried")
    p.add_argument("--output", required=True)
    a = p.parse_args()
    d = Path(a.eval_dir) if a.eval_dir else Path(a.run) / "vcf_chr1" / "eval_bed" / "vcfeval" / "raw_calls"
    graph = set()
    with (gzip.open(a.graph_alleles, "rt") if a.graph_alleles else open("/dev/null")) as f:
        for line in f:
            if line.startswith("chr1\t"):
                c, pos, ref, alt = line.rstrip("\n").split("\t")
                graph.add((c, int(pos), ref, alt))
    n, tp = Counter(), Counter()  # (VAF value, in graph) -> truth, called
    for name, called in (("tp-baseline.vcf.gz", True), ("fn.vcf.gz", False)):
        for key, af in records_af(d / name, with_key=True):
            if af is None or not 0.0099 <= af <= 0.1651:
                continue
            k = (round(af, 4), key in graph)
            n[k] += 1
            tp[k] += called
    name = a.title or (Path(a.run).name if a.run else str(d))
    lines = [f"# E23 at every truth VAF value: {name}, {a.kind}, chr1, BED, no population filtering", "",
             "| truth VAF | not graph-carried | recall | all truth | recall |", "|---|---|---|---|---|"]
    for v in sorted({v for v, _ in n}):
        no, al = n[v, False], n[v, False] + n[v, True]
        rn = f"{tp[v, False] / no:.3f}" if no else "-"
        lines.append(f"| {v:.4f} | {no:,} | {rn} | {al:,} | {(tp[v, False] + tp[v, True]) / al:.3f} |")
    ctp, cfp = Counter(), Counter()
    step = lambda af: min(int(af * 100), 16)  # noqa: E731  1 % steps, the last one 16-16.5 %+
    for name, counter in (("tp.vcf.gz", ctp), ("fp.vcf.gz", cfp)):
        for af in records_af(d / name):
            if af is not None and 0.01 <= af <= 0.1651:
                counter[step(af)] += 1
    lines += ["", "| call AF | calls | precision |", "|---|---|---|"]
    for s in range(1, 17):
        c = ctp[s] + cfp[s]
        if c:
            lines.append(f"| {s}-{s + 1 if s < 16 else 16.5} % | {c:,} | {ctp[s] / c:.3f} |")
    ranges = [(0.0099, 0.0199, "1-2 %"), (0.0199, 0.0499, "2-5 %"), (0.0499, 0.0999, "5-10 %"), (0.0999, 0.1651, "10-16.5 %"),
              (0.0099, 0.1651, "1-16.5 %")]
    rtp, rcal = Counter(), Counter()  # precision by call AF in the same ranges
    for name, counter in (("tp.vcf.gz", rtp), ("fp.vcf.gz", None)):
        for af in records_af(d / name):
            for lo, hi, r in ranges:
                if af is not None and lo <= af < hi:
                    rcal[r] += 1
                    if counter is not None:
                        counter[r] += 1
    lines += ["", "| VAF | recall (not graph-carried) | recall (all truth) | precision (call AF) |", "|---|---|---|---|"]
    for lo, hi, r in ranges:
        keys = [k for k in n if lo <= k[0] < hi]
        no = sum(n[k] for k in keys if not k[1]); tno = sum(tp[k] for k in keys if not k[1])
        al = sum(n[k] for k in keys); tal = sum(tp[k] for k in keys)
        lines.append(f"| {r} | {tno / no:.3f} ({tno:,} / {no:,}) | {tal / al:.3f} ({tal:,} / {al:,}) | "
                     f"{rtp[r] / rcal[r]:.3f} ({rtp[r]:,} / {rcal[r]:,}) |")
    Path(a.output).write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-7:]))


if __name__ == "__main__":
    main()

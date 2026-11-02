#!/usr/bin/env python3
"""Why the R1.5 conditional HapMap denominators (215,094 SNV / 92,079 INDEL) are smaller than truth minus graph-carried
(216,794 / 98,561): classify every truth record (chr1, BED, truth VAF 1-16.5 %) that is NOT graph-carried but is absent
from the eval_chr1.bed evaluation (eval_bed_nograph_all), by its relation to the masked graph-carried spans.

Inputs (all existing):
  full-BED eval   <run>/vcf_chr1/eval_bed/vcfeval/raw_calls/{tp-baseline,fn}.vcf.gz
  eval-BED eval   <run>/vcf_chr1/eval_bed_nograph_all/vcfeval/raw_calls/{tp-baseline,fn}.vcf.gz
  graph alleles   E23_graph_alleles.mgi_non_HG005.d9.atomized.tsv.gz (the definition of call_vaf_by_value.py)
  masked alleles  E23_nograph/graph_carried_truth.tsv.gz (the holes of exclude_graph_carried.py: [pos-2, pos+len(ref)) 0-based)
Also checks Mutect2's eval-BED evaluation has the same truth records.
"""
import bisect
import gzip
import json
from collections import Counter
from pathlib import Path

R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
X = Path("/scratch/jshen/data/reviewer_experiments_20261013/E23_nograph")
G = Path("/scratch/jshen/data/reviewer_experiments_20261013/E23_graph_alleles.mgi_non_HG005.d9.atomized.tsv.gz")
M2 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E23_mutect2")
TRUTH_VCF = Path("/scratch/jshen/data/pansoma_v2_tensors/HapMap_truth/HapMap_somatic_snv_indel.vcf.gz")
RUNS = {"SNV": "HapMap_Illumina_SNV_full", "INDEL": "HapMap_Illumina_INDEL"}


def records(path):
    out = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            v = dict(zip(c[8].split(":"), c[9].split(":")))
            af = v.get("AF", ".")
            out[(c[0], int(c[1]), c[3], c[4])] = None if af in (".", "") else float(af)
    return out


def truth_af():
    """VAF (FORMAT AF) of every chr1 truth record, both kinds (for the masking alleles' VAF)."""
    out = {}
    with gzip.open(TRUTH_VCF, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.startswith("chr1\t"):
                continue
            c = line.rstrip("\n").split("\t")
            v = dict(zip(c[8].split(":"), c[9].split(":")))
            out[("chr1", int(c[1]), c[3], c[4])] = float(v["AF"]) if v.get("AF", ".") not in (".", "") else None
    return out


def in_range(af):
    return af is not None and 0.0099 <= af <= 0.1651


def main():
    graph = set()
    with gzip.open(G, "rt") as f:
        for line in f:
            if line.startswith("chr1\t"):
                c, p, r, a = line.rstrip("\n").split("\t")
                graph.add((c, int(p), r, a))
    masked = []
    with gzip.open(X / "graph_carried_truth.tsv.gz", "rt") as f:
        for line in f:
            if line.startswith("chr1\t"):
                c, p, r, a = line.rstrip("\n").split("\t")
                masked.append((int(p), r, a))
    masked.sort()
    taf = truth_af()
    mpos = {}
    for p, r, a in masked:
        mpos.setdefault(p, []).append((r, a))
    # masked allele REF spans, 0-based [s, e), and the +-1 bp holes
    spans = sorted((p - 1, p - 1 + len(r), p, r, a) for p, r, a in masked)
    starts = [s[0] for s in spans]
    maxlen = max(e - s for s, e, *_ in spans)
    bed = [tuple(map(int, l.split("\t")[1:3])) for l in open(X / "eval_chr1.bed")]
    bstarts = [b[0] for b in bed]

    def in_eval(x):  # 0-based base x inside eval_chr1.bed
        i = bisect.bisect_right(bstarts, x) - 1
        return i >= 0 and bed[i][0] <= x < bed[i][1]

    def near(s, e, pad):
        """masked alleles whose span [ms, me) padded by pad overlaps [s, e)."""
        i = bisect.bisect_left(starts, s - maxlen - pad - 1)
        out = []
        while i < len(spans) and spans[i][0] < e + pad:
            ms, me = spans[i][0], spans[i][1]
            if ms - pad < e and s < me + pad:
                out.append(spans[i])
            i += 1
        return out

    result = {}
    for kind, run in RUNS.items():
        full = {}
        for n in ("tp-baseline.vcf.gz", "fn.vcf.gz"):
            full.update(records(R / run / "vcf_chr1/eval_bed/vcfeval/raw_calls" / n))
        tp_full = set(records(R / run / "vcf_chr1/eval_bed/vcfeval/raw_calls/tp-baseline.vcf.gz"))
        ev = {}
        for n in ("tp-baseline.vcf.gz", "fn.vcf.gz"):
            ev.update(records(R / run / "vcf_chr1/eval_bed_nograph_all/vcfeval/raw_calls" / n))
        m2 = {}
        for n in ("tp-baseline.vcf.gz", "fn.vcf.gz"):
            m2.update(records(M2 / f"eval.pass.{kind}" / n))
        F = {k for k, af in full.items() if in_range(af)}
        E = {k for k, af in ev.items() if in_range(af)}
        M = {k for k, af in m2.items() if in_range(af)}
        GC = F & graph
        nonG = F - graph
        D = nonG - E
        cat, why_kind, why_vaf, det = Counter(), Counter(), Counter(), Counter()
        examples = {}
        for k in sorted(D):
            _, p, r, a = k
            s, e = p - 1, p - 1 + len(r)
            same = [m for m in mpos.get(p, []) if m != (r, a)]
            ov = near(s, e, 0)
            adj = near(s, e, 1)
            if same:
                c = "same VCF position as a graph-carried allele (another allele of the same site)"
                hit = [(p, mr, ma) for mr, ma in same]
            elif ov:
                c = "REF span overlaps a graph-carried allele's REF span (different position)"
                hit = [(x[2], x[3], x[4]) for x in ov]
            elif adj:
                c = "within 1 bp of a graph-carried allele (the +-1 bp margin of the mask)"
                hit = [(x[2], x[3], x[4]) for x in adj]
            else:
                c = "no graph-carried allele within 1 bp"
                hit = []
            # the record's own first base and whole span vs the evaluation BED
            first_in = in_eval(s)
            whole_in = all(in_eval(x) for x in range(s, e))
            det[(c, "first base in eval BED" if first_in else "first base masked",
                 "whole REF span in eval BED" if whole_in else "REF span partly masked")] += 1
            cat[c] += 1
            hk = sorted({"SNV" if len(hr) == len(ha) else "INDEL" for _, hr, ha in hit})
            why_kind[c + " | masking allele kind: " + "+".join(hk or ["-"])] += 1
            hv = [taf.get(("chr1", hp, hr, ha)) for hp, hr, ha in hit]
            why_vaf[c + " | a masking allele at VAF 1-16.5 %: " + str(any(in_range(v) for v in hv))] += 1
            examples.setdefault(c, []).append((k, hit[:2]))
        # retained records that overlap a mask (tests rtg -e semantics)
        kept_touch = sum(1 for k in E if near(k[1] - 1, k[1] - 1 + len(k[2]), 1) and k not in graph)
        res = dict(
            truth_F=len(F), graph_carried=len(GC), nonG=len(nonG), eval_E=len(E), E_minus_F=len(E - F),
            E_and_graph=len(E & graph), extra_D=len(D), mutect2_eval_equals_pansoma_eval=(M == E), mutect2_eval_n=len(M),
            D_tp_in_full_eval=len(D & tp_full), nonG_tp_full=len(nonG & tp_full),
            categories=dict(cat), by_mask_kind=dict(why_kind), by_mask_vaf=dict(why_vaf),
            detail={" / ".join(k): v for k, v in det.items()}, kept_records_within_1bp_of_mask=kept_touch,
            examples={c: [[list(k), [list(h) for h in hit]] for k, hit in v[:5]] for c, v in examples.items()},
        )
        result[kind] = res
        print(kind, json.dumps({k: v for k, v in res.items() if k != "examples"}, indent=1), flush=True)
    out = Path(__file__).with_suffix(".json")
    out.write_text(json.dumps(result, indent=1) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()

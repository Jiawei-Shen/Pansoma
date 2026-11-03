#!/usr/bin/env python3
"""Why Pansoma's HCC1395T chr1 SNV F1 is not the highest: truth SNVs found by Pansoma (six-tumor / four-line models,
`--norm dataset` re-predictions `vcf_chr1_HCC1395T_<p>_d46_normfix`) and by DeepSomatic-TO / ClairS-TO (E8 scores,
`eval_bed/pass_pon`), all after the same population filter in the SEQC2 high-confidence regions, and the candidate
status of each truth SNV in Pansoma's tensor set (`tensors/somatic.recall.tsv`: tensor, or filtered at the builder's
floors: allele fraction 0.08 / three ALT reads, or no candidate). Writes tables.md and summary.json next to itself.
"""
import gzip
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
T = Path("/scratch/jshen/data/pansoma_v2_tensors")


def recs(vcf):
    out = set()
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c, p, _, r, a = line.split("\t")[:5]
            if c == "chr1" and len(r) == 1 and all(len(x) == 1 for x in a.split(",")):
                out.add((int(p), r, a))
    return out


def tp_fn(d):
    return recs(d / "tp-baseline.vcf.gz"), recs(d / "fn.vcf.gz")


def status(p):
    s = {}
    with open(T / f"HCC1395T_{p}_d46/tensors/somatic.recall.tsv") as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            x = dict(zip(h, line.rstrip("\n").split("\t")))
            if x["chrom"] == "chr1" and x["kind"] == "SNP":
                st = x["status"]
                if st == "filtered":
                    st = "filtered:" + (x["detail"] or "?")
                s[int(x["vcf_pos"])] = st
    return s


def main():
    res, L = {}, ["# HCC1395T chr1 SNVs: where Pansoma's recall is lost relative to the comparators", ""]
    for p in ("Illumina", "PacBio", "ONT"):
        st = status(p)
        calls = {}
        for m in ("multicancer6", "multicancer4"):
            calls[m] = tp_fn(R / f"{m}_{p}_SNV/vcf_chr1_HCC1395T_{p}_d46_normfix/eval_bed/vcfeval/pon_calls")
        for c in ("deepsomatic", "clairs_to"):
            calls[c] = tp_fn(E8 / c / f"HCC1395T_{p}/SNV/eval_bed/pass_pon")
        truth = calls["multicancer6"][0] | calls["multicancer6"][1]
        for k, (tp, fn) in calls.items():
            assert tp | fn == truth, (p, k, len(tp | fn), len(truth))
        stc = Counter(st.get(x[0], "not in table") for x in truth)
        row = dict(truth=len(truth), status=dict(stc))
        P = calls["multicancer6"][0]
        for c in ("deepsomatic", "clairs_to"):
            C = calls[c][0]
            only_c = C - P
            row[c] = dict(recall=len(C) / len(truth), pansoma_recall=len(P) / len(truth), both=len(C & P),
                          comparator_only=len(only_c), pansoma_only=len(P - C),
                          comparator_only_status=dict(Counter(st.get(x[0], "not in table") for x in only_c)))
        res[p] = row
    (HERE / "summary.json").write_text(json.dumps(res, indent=1))
    order = ["tensor_representative", "filtered:min_af", "filtered:min_variants", "no_candidate"]
    name = {"tensor_representative": "tensor (scored by the model)", "filtered:min_af": "candidate below allele fraction 0.08",
            "filtered:min_variants": "candidate with fewer than 3 ALT reads", "no_candidate": "no candidate"}
    L += ["| chr1 SNVs in the BED | " + " | ".join(res) + " |", "|---|" + "---|" * len(res)]
    L.append("| truth SNVs | " + " | ".join(str(res[p]["truth"]) for p in res) + " |")
    for k in order:
        L.append(f"| {name[k]} | " + " | ".join(f"{res[p]['status'].get(k, 0)} ({res[p]['status'].get(k, 0) / res[p]['truth']:.1%})" for p in res) + " |")
    for c, cn in (("deepsomatic", "DeepSomatic-TO"), ("clairs_to", "ClairS-TO")):
        L.append(f"| recall: six-tumor Pansoma / {cn} | " + " | ".join(f"{res[p][c]['pansoma_recall']:.3f} / {res[p][c]['recall']:.3f}" for p in res) + " |")
        L.append(f"| found by {cn} only / by Pansoma only | " + " | ".join(f"{res[p][c]['comparator_only']} / {res[p][c]['pansoma_only']}" for p in res) + " |")
        for k in order:
            L.append(f"| — {cn} only: {name[k]} | " + " | ".join(str(res[p][c]["comparator_only_status"].get(k, 0)) for p in res) + " |")
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

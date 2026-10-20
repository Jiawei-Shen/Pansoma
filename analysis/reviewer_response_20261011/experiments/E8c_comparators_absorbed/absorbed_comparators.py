"""E8c (R1.4): do the linear tumor-only callers keep the somatic alleles that the HPRC v1.1 graph absorbs?

Absorbed alleles: the "perfect bypass" truths of the graph-absorbed analyses (ALT reads align without edit to an
existing HPRC v1.1 (d9) path, so Pansoma gets no candidate): HG008T `analysis/graph_absorbed_somatic_20261001/
per_variant.tsv` (<platform>_perfect), COLO829T `analysis/graph_absorbed_somatic_colo829t_20261005/per_variant.tsv`
(<platform>_perfect). Per dataset, the alleles absorbed on that platform.
Comparators: the user's existing tumor-only runs re-scored by E8 (DeepSomatic 1.9.0; ClairS-TO 0.4.2 on HG008T, 0.4.1
on COLO829T; chr1 only), rtg vcfeval in the benchmark BED at four stages (`score/<caller>/<set>/<KIND>/eval_bed/`):
all records (any FILTER), all records after our PoN, PASS, PASS after our PoN. A truth allele is kept at a stage when
it is in that stage's tp-baseline; it is in the BED when it is in tp-baseline or fn of the `all` stage.
With --graph d46 (2026-10-20): the alleles the HPRC v2.1 (d46) graph absorbs, from E5b (`../E5b_absorbed_d46/
lost_alleles.tsv.gz`): no candidate on d46 on that platform (`<platform>_d46` == N) and the ALT spelled by a d46 path
(`d46_spelled`, the strict test of E5b s4, a lower bound); also HCC1395T (comparators run on tequila 2026-10-19/20).
Writes tables.md and summary.json (d9), tables_d46.md and summary_d46.json (d46) next to this script. Standard library
+ gzip only.
usage: python3 absorbed_comparators.py [--graph d9|d46]
"""
import argparse
import csv
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SCORE = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
ABSORBED = {"HG008T": REPO / "analysis/graph_absorbed_somatic_20261001/per_variant.tsv",
            "COLO829T": REPO / "analysis/graph_absorbed_somatic_colo829t_20261005/per_variant.tsv"}
DATASETS = [("HG008T", "Illumina", "HG008T_Illumina", "Illumina"), ("HG008T", "PacBio HiFi", "HG008T_PacBio", "PacBio"),
            ("HG008T", "ONT", "HG008T_ONT", "ONT"), ("COLO829T", "Illumina", "COLO829T_Illumina", "Illumina"),
            ("COLO829T", "Fiber-seq", "COLO829T_fiberseq", "fiberseq"), ("COLO829T", "ONT", "COLO829T_ONT", "ONT")]
DATASETS_D46 = DATASETS + [("HCC1395T", "Illumina", "HCC1395T_Illumina", "Illumina"),
                           ("HCC1395T", "PacBio HiFi", "HCC1395T_PacBio", "PacBio"), ("HCC1395T", "ONT", "HCC1395T_ONT", "ONT")]
E5B = HERE.parent / "E5b_absorbed_d46" / "lost_alleles.tsv.gz"
CALLERS = [("deepsomatic", "DeepSomatic-TO"), ("clairs_to", "ClairS-TO")]
STAGES = [("all", "any record"), ("pass", "PASS"), ("pass_pon", "PASS + our PoN")]


def keys(vcf):
    out = set()
    if not vcf.exists():
        return None
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] != "#":
                c = line.split("\t", 5)
                for alt in c[4].split(","):
                    out.add((c[0], int(c[1]), c[3], alt))
    return out


def absorbed(sample, col):
    out = {"SNV": set(), "INDEL": set()}
    with open(ABSORBED[sample]) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["chrom"] == "chr1" and r.get(f"{col}_perfect") == "True":
                kind = r.get("kind") or r.get("kind2")
                out["SNV" if kind == "SNV" else "INDEL"].add((r["chrom"], int(r["vcf_pos"]), r["vcf_ref"], r["vcf_alt"]))
    return out


def absorbed_d46(sample, col):
    out = {"SNV": set(), "INDEL": set()}
    with gzip.open(E5B, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["sample"] == sample and r["chrom"] == "chr1" and r.get(f"{col}_d46") == "N" and r["d46_spelled"] == "1":
                out["SNV" if r["kind"] == "SNP" else "INDEL"].add((r["chrom"], int(r["vcf_pos"]), r["vcf_ref"], r["vcf_alt"]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", choices=["d9", "d46"], default="d9")
    graph = ap.parse_args().graph
    rows = []
    for sample, plat, ds, col in (DATASETS if graph == "d9" else DATASETS_D46):
        ab = absorbed(sample, col) if graph == "d9" else absorbed_d46(sample, col)
        for kind in ("SNV", "INDEL"):
            for caller, label in CALLERS:
                d = SCORE / caller / ds / kind / "eval_bed"
                tp = {st: keys(d / st / "tp-baseline.vcf.gz") for st, _ in STAGES}
                fn_all = keys(d / "all" / "fn.vcf.gz")
                if tp["all"] is None or fn_all is None:
                    continue
                truth = tp["all"] | fn_all
                a = ab[kind] & truth
                row = dict(sample=sample, platform=plat, kind=kind, caller=label, truth=len(truth), absorbed=len(a))
                for st, _ in STAGES:
                    t = tp[st] or set()
                    row[f"{st}_absorbed"] = len(a & t)
                    row[f"{st}_all"] = len(truth & t)
                rows.append(row)
    pct = lambda n, d: f"{n} ({100 * n / d:.0f}%)" if d else "-"  # noqa: E731
    name = "HPRC v1.1 (d9)" if graph == "d9" else "HPRC v2.1 (d46)"
    rule = ("ALT reads align without edit to an existing HPRC v1.1 (d9) path on that platform (no Pansoma candidate)."
            if graph == "d9" else
            "no candidate on HPRC v2.1 (d46) on that platform and the ALT spelled by a d46 path (E5b strict test, a lower bound).")
    L = [f"# E8c: linear tumor-only callers at the alleles the {name} graph absorbs (chr1, benchmark BED)", "",
         "Absorbed = " + rule,
         "Kept = the truth allele is a true positive of the caller at that stage (rtg tp-baseline). In brackets: share of the",
         "absorbed alleles, and for comparison the share of all chr1 truth alleles of that kind in the BED.", "",
         "| tumor | platform | kind | caller | absorbed / all truth | any record: absorbed (all truth) | PASS | PASS + our PoN |",
         "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        cells = [f"{pct(r[f'{st}_absorbed'], r['absorbed'])} ({100 * r[f'{st}_all'] / r['truth']:.0f}%)" for st, _ in STAGES]
        L.append(f"| {r['sample']} | {r['platform']} | {r['kind']} | {r['caller']} | {r['absorbed']} / {r['truth']} | " + " | ".join(cells) + " |")
    suf = "" if graph == "d9" else "_d46"
    (HERE / f"tables{suf}.md").write_text("\n".join(L) + "\n")
    (HERE / f"summary{suf}.json").write_text(json.dumps(rows, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

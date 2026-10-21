"""E9b x E9c (R3.4 point 5): which graph-coordinate PASS SNV calls a linear reference can represent, by what they are.

Rows: the unplaced PASS SNV calls (linear_vcf records without a GRCh38 position) of each in-sample graph run that has
both an E9b projection (`<scratch>/runs_e9b/<run>/vcf_chr1/projected.tsv`, chm13_class; `<scratch>/runs/` for runs
projected after the CHM13 path was built, 2026-10-21) and E9c categories
(`../E9c_offref_filter_validation/per_call.tsv.gz`, group unplaced_pass; COLO829T ONT d46: E9c's random 3,000).
Linear reference: CHM13 = the event's node is on the CHM13 path (chm13_class direct); GRCh38 = projection class i
(a GRCh38 allele from the flanking GRCh38 nodes, E9a); neither = the rest (specific to the graph).
E9c category (exclusive, in this order): population (graph PoN: an HPRC allele), germline (in the normal's assembly),
tumor_only / not_in_normal (assembly support for somatic), unsupported, unverifiable.
Writes results/chm13_by_category.tsv and prints per run and pooled shares.
usage: python3 chm13_by_category.py
"""
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
S = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection")
E9C = HERE.parent / "E9c_offref_filter_validation" / "per_call.tsv.gz"
CATS = ("population", "germline", "tumor_only", "not_in_normal", "unsupported", "unverifiable", "no_tensor")
REFS = ("CHM13", "GRCh38 only", "neither")


def main():
    cat = defaultdict(dict)
    with gzip.open(E9C, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["run"].startswith("rcpu__") and r["run"].endswith("__vcf_chr1") and r["group"] == "unplaced_pass" \
                    and "__norescue" not in r["run"]:
                cat[r["run"][len("rcpu__"):-len("__vcf_chr1")]][r["candidate_id"]] = r["category"]
    out, pooled = [], defaultdict(Counter)
    for run in sorted(cat):
        p = next((q for q in (S / d / run / "vcf_chr1" / "projected.tsv" for d in ("runs_e9b", "runs")) if q.exists()), None)
        if p is None:
            print(f"{run}: no projection, skipped")
            continue
        c = defaultdict(Counter)
        with open(p) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["chm13_class"] == "not_built":
                    raise SystemExit(f"{p}: CHM13 classes not built")
                if r["filter"] != "PASS" or r["candidate_id"] not in cat[run]:
                    continue
                ref = "CHM13" if r["chm13_class"] == "direct" else "GRCh38 only" if r["proj_class"] == "i" else "neither"
                k = cat[run][r["candidate_id"]]
                c[ref][k] += 1
                pooled[ref][k] += 1
        for ref in REFS:
            n = sum(c[ref].values())
            out.append([run, ref, n] + [c[ref][k] for k in CATS])
        print(run, " | ".join(f"{ref}: {sum(c[ref].values())} " +
                                ", ".join(f"{k} {c[ref][k]}" for k in CATS if c[ref][k]) for ref in REFS))
    for ref in REFS:
        n = sum(pooled[ref].values())
        out.append(["pooled", ref, n] + [pooled[ref][k] for k in CATS])
        print("pooled", ref, n, ", ".join(f"{k} {pooled[ref][k]} ({100 * pooled[ref][k] / n:.1f}%)"
                                         for k in CATS if pooled[ref][k]))
    (HERE / "results").mkdir(exist_ok=True)
    with open(HERE / "results" / "chm13_by_category.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "linear reference", "calls"] + list(CATS))
        w.writerows(out)


if __name__ == "__main__":
    main()

"""E9b: the CHM13 VCFs (pansoma_net_v2 linear_vcf on the graph's CHM13 reference path) of the runs next to their GRCh38
VCFs: which PASS calls have a position on GRCh38, on CHM13, on both, on neither. A PASS call is a PASS record of the graph
VCF (one candidate ID); it is placed on a reference when its ID is in that reference's VCF (in any record: merging can
put it into another record's ID list, with that record's FILTER).

    python collect_chm13_vcf.py     (writes results/chm13_vcf_by_run.tsv)

Inputs: <runs>/<RUN>/vcf_chr1/<set>.<KIND>.linear{,.unplaced}.vcf.gz and
/scratch/jshen/data/reviewer_experiments_20261013/E9b_chm13_vcf/<RUN>/vcf_chr1/<set>.<KIND>.chm13{,.unplaced}.vcf.gz
(tmp/reviewer_cpu_20261019/E9b_chm13/chm13_vcf.sbatch).
"""
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
CHM13 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9b_chm13_vcf")


def ids(path, pass_only=False):
    out = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                x = line.split("\t", 8)
                if not pass_only or x[6] == "PASS":
                    out.update(x[2].split(";"))
    return out


def main():
    rows = []
    for d in sorted(CHM13.glob("*/vcf_chr1")):
        run = d.parent.name
        (chm13,) = d.glob("*.chm13.vcf.gz")
        name = chm13.name[: -len(".chm13.vcf.gz")]
        stats = json.loads((d / f"{name}.chm13.stats.json").read_text())
        calls = ids(RUNS / run / "vcf_chr1" / f"{name}.graph.vcf.gz", pass_only=True)
        g, c = ids(RUNS / run / "vcf_chr1" / f"{name}.linear.vcf.gz"), ids(chm13)
        for path in (RUNS / run / "vcf_chr1" / f"{name}.linear.unplaced.vcf.gz", d / f"{name}.chm13.unplaced.vcf.gz"):
            if ids(path) & (g if "linear" in path.name else c):
                raise SystemExit(f"{path}: a candidate is both placed and unplaced")
        rows.append(dict(run=run, kind=name.rsplit(".", 1)[1], graph="d46" if "_d46_" in run else "d9",
                         graph_records=stats["graph_records"], chm13_records=stats["records"],
                         chm13_unplaced=stats["unplaced"], pass_calls=len(calls), pass_both=len(calls & g & c),
                         pass_grch38_only=len((calls & g) - c), pass_chm13_only=len((calls & c) - g),
                         pass_neither=len(calls - g - c)))
    out = HERE / "results" / "chm13_vcf_by_run.tsv"
    cols = list(rows[0])
    out.write_text("\t".join(cols) + "\n" + "".join("\t".join(str(r[k]) for k in cols) + "\n" for r in rows))
    for r in rows:
        print("\t".join(str(r[k]) for k in cols))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

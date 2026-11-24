"""ver4 offref (2): Response Table 3.11, CHM13 representation of the graph-coordinate SNV PASS calls of the nine
six-tumor call sets (HPRC v2.1, chr1, default output).

CHM13 allele = E9b `chm13_class == direct` of the E9a / E9b projection (experiments/E9a_unplaced_projection/
project_unplaced.py, unmodified; job tmp/ver4_20261123/offref/r34_task.sh): the event's node(s) have one visit on the
embedded T2T-CHM13 path, so the call can be written as a CHM13 record. Definition of tmp/r33_hcc1395_20261104/r34/
r34_results.py (the ver3 HCC1395T rows). Independent count, asserted equal for every call set: PASS calls of the graph
VCF that the CHM13 VCF (linear_vcf on the d46 CHM13 path) places and the GRCh38 VCF does not.

usage: python table_3_11.py           -> table_3.11.md, table_3.11.json
       python table_3_11.py --repro   -> the ver3 HCC1395T rows from the 2026-11-05 outputs (prints; writes nothing)
"""
import csv
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
X = Path("/scratch/jshen/data/reviewer_experiments_20261013")
PLAT = {"Illumina": "Illumina", "PacBio_HiFi": "PacBio HiFi", "Fiber-seq": "Fiber-seq", "ONT": "ONT"}
csv.field_size_limit(sys.maxsize)


def ids(path, pass_only=False):
    out = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                x = line.split("\t", 8)
                if not pass_only or x[6] == "PASS":
                    out.update(x[2].split(";"))
    return out


def one(name, vdir, o):
    vdir = Path(vdir)
    pj = json.loads((o / "e9a" / name / "projection.json").read_text())
    proj = [r for r in csv.DictReader(open(o / "e9a" / name / "projected.tsv"), delimiter="\t") if r["filter"] == "PASS"]
    u = len(proj)
    rep = json.loads((vdir / "eval_nobed/report.json").read_text())
    assert u == pj["unplaced_pass"] == rep["unplaced_pass"], (name, u, pj["unplaced_pass"], rep["unplaced_pass"])
    assert Path(pj["vcf_dir"]).resolve() == vdir.resolve(), (name, pj["vcf_dir"])
    cls = {k: sum(r["chm13_class"] == k for r in proj) for k in ("direct", "anchor_interval", "none")}
    assert sum(cls.values()) == u and cls == {k: pj["chm13_pass"].get(k, 0) for k in cls}, (name, cls, pj["chm13_pass"])
    (gv,) = vdir.glob("*.graph.vcf.gz")
    stem = gv.name[: -len(".graph.vcf.gz")]
    calls = ids(gv, pass_only=True)
    g = ids(vdir / f"{stem}.linear.vcf.gz")
    c = ids(o / "chm13_vcf" / name / f"{stem}.chm13.vcf.gz")
    chm13_only = len((calls & c) - g)
    assert chm13_only == cls["direct"], (name, chm13_only, cls["direct"])
    return dict(name=name, calls=u, chm13_allele=cls["direct"], chm13_anchor_interval=cls["anchor_interval"],
                chm13_none=cls["none"], chm13_vcf_only=chm13_only, grch38_class_i=sum(r["proj_class"] == "i" for r in proj),
                grch38_or_chm13=sum(r["proj_class"] == "i" or r["chm13_class"] == "direct" for r in proj))


def cell(a, n):
    p = 100 * a / n
    return f"{a} ({'100' if a == n else f'{p:.1f}'}%)"


def main():
    if "--repro" in sys.argv:
        o = X / "r34_hcc1395t_20261105"
        for line in open("/scratch/jshen/Github/Pansoma/tmp/r33_hcc1395_20261104/r34/tasks.txt"):
            name, vdir = line.split()
            if name.endswith("_off"):
                r = one(name, vdir, o)
                print(name, r["calls"], cell(r["chm13_allele"], r["calls"]), "CHM13 VCF check", r["chm13_vcf_only"])
        return
    o = X / "ver4_20261123/offref/r34"
    rows = []
    for line in open(HERE / "runs.txt"):
        name, vdir, tumor, plat = line.split()
        rows.append(dict(tumor=tumor, platform=PLAT[plat], **one(name, vdir, o)))
    md = ["| **Tumor** | **Platform** | **Calls** | **CHM13 allele** |", "|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['tumor']} | {r['platform']} | {r['calls']} | {cell(r['chm13_allele'], r['calls'])} |")
    tot = {}
    for label, sel in (("all", None), ("excluding HCC1395T", ("HG008T", "COLO829T")), ("HG008T", ("HG008T",)),
                       ("COLO829T", ("COLO829T",)), ("HCC1395T", ("HCC1395T",))):
        s = [r for r in rows if sel is None or r["tumor"] in sel]
        tot[label] = dict(calls=sum(r["calls"] for r in s), chm13_allele=sum(r["chm13_allele"] for r in s),
                          chm13_anchor_interval=sum(r["chm13_anchor_interval"] for r in s),
                          chm13_none=sum(r["chm13_none"] for r in s), grch38_class_i=sum(r["grch38_class_i"] for r in s),
                          grch38_or_chm13=sum(r["grch38_or_chm13"] for r in s))
    a, e = tot["all"], tot["excluding HCC1395T"]
    md.append(f"| **Total** | | **{a['calls']}** | **{cell(a['chm13_allele'], a['calls'])}** |")
    md.append(f"| Total excluding HCC1395T | | {e['calls']} | {cell(e['chm13_allele'], e['calls'])} |")
    table = "\n".join(md) + "\n"
    notes = ["", "Per tumor: " + "; ".join(f"{t} {cell(tot[t]['chm13_allele'], tot[t]['calls'])} of {tot[t]['calls']}"
                                           for t in ("HG008T", "COLO829T", "HCC1395T")) + ".",
             f"Other calls (all nine): CHM13 anchor interval only {a['chm13_anchor_interval']}, no CHM13 position {a['chm13_none']}; "
             f"GRCh38 allele through the bubble's flanking nodes (E9a class i) {a['grch38_class_i']}; "
             f"GRCh38 (class i) or CHM13 allele {cell(a['grch38_or_chm13'], a['calls'])}.",
             "CHM13 VCF check (PASS calls placed on CHM13 but not on GRCh38 = CHM13 allele by projection): " +
             ", ".join(f"{r['tumor']} {r['platform']} {r['chm13_vcf_only']} / {r['chm13_allele']}" for r in rows) + "."]
    print(table + "\n".join(notes))
    (HERE / "table_3.11.md").write_text(
        "# Response Table 3.11 (ver4): CHM13 representation of graph-coordinate SNV PASS calls\n\n"
        "Six-tumor current SNV models (snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, seed 0), HPRC v2.1, chr1, "
        "default output. Definitions in table_3_11.py.\n\n" + table + "\n".join(notes) + "\n")
    (HERE / "table_3.11.json").write_text(json.dumps(dict(rows=rows, totals=tot), indent=1) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Letter Supplementary Tables R1.6b (SNVs) and R1.6c (INDELs): why truth alleles get no candidate, HPRC v2.1 (d46),
chr1-22, all truth alleles (inside and outside the benchmark regions). User 2026-11-10: add the HCC1395T columns.

Inputs (read-only)
  /scratch/jshen/data/pansoma_v2_tensors/<SET>/tensors/somatic.recall.tsv   truth / tensor (T + A) / filtered / no_candidate
  <SET>/snv_no_candidate.tsv, <SET>/indel_no_candidate.tsv (this directory)  final_class of every no_candidate allele
                                                                             (scripts/report.py, classes S1-S4c, I1-I9)
The no_candidate count of somatic.recall.tsv must equal the number of classified alleles. The HG008T / COLO829T cells
are checked against the published tables (EXPECTED below). Sets whose classification is not finished are skipped.
Second SNV table on the current SNV pipeline (user 2026-11-10/11: "我们现在用的是AF 0.05的"): PacBio HiFi / Fiber-seq /
Illumina sets of COLO829T and HCC1395T are built at SNV AF 0.05 (<SET>_snvaf05); HG008T (not rebuilt, user 11-09) and
ONT (original recipe) stay at 0.08. Only the tensor / floor rows change: the no_candidate truth ids of each _snvaf05 set
are asserted identical to those of its 0.08 set, so the S1-S4 classes carry over.
Writes r16bc_tables.md here.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

csv.field_size_limit(sys.maxsize)
HERE = Path(__file__).resolve().parent
T = Path("/scratch/jshen/data/pansoma_v2_tensors")
SETS = [("HG008T_Illumina_d46", "HG008T Illumina"), ("HG008T_PacBio_d46", "HG008T PacBio HiFi"), ("HG008T_ONT_d46", "HG008T ONT"),
        ("COLO829T_Illumina_d46", "COLO829T Illumina"), ("COLO829T_fiberseq_d46", "COLO829T Fiber-seq"), ("COLO829T_ONT_d46", "COLO829T ONT"),
        ("HCC1395T_Illumina_d46", "HCC1395T Illumina"), ("HCC1395T_PacBio_d46", "HCC1395T PacBio HiFi"), ("HCC1395T_ONT_d46", "HCC1395T ONT")]
TENSOR = {"tensor_representative", "tensor_non_representative_allele"}
# published cells (letter ver3 Suppl. Tables R1.6b / R1.6c), HG008T x 3 then COLO829T x 3
EXPECTED = {
    "SNV": dict(truth=[8690] * 3 + [42035] * 3, tensor=[8311, 8252, 8293, 37862, 37736, 37515], floor=[84, 52, 68, 3452, 3673, 3997],
                none=[295, 386, 329, 721, 626, 523], S1=[51, 63, 1, 244, 265, 15], S2=[30, 25, 25, 278, 118, 367], S3=[60, 51, 36, 55, 50, 44],
                S4a=[39, 65, 58, 63, 133, 20], S4b=[72, 133, 160, 7, 11, 12], S4c=[43, 49, 49, 74, 49, 65]),
    "INDEL": dict(truth=[8496] * 3 + [1912] * 3, tensor=[3924, 4623, 4608, 1169, 1300, 1342], floor=[811, 704, 1088, 100, 149, 180],
                  none=[3761, 3169, 2800, 643, 463, 390], I1=[1639, 1259, 1022, 365, 240, 199], I2=[1597, 1680, 1568, 158, 145, 133],
                  I3=[108, 109, 125, 57, 36, 37], I4=[300, 53, 44, 35, 15, 8], I5=[101, 30, 2, 20, 18, 3], other=[16, 38, 39, 8, 9, 10]),
}


AF05 = {"COLO829T_Illumina_d46", "COLO829T_fiberseq_d46", "HCC1395T_Illumina_d46", "HCC1395T_PacBio_d46"}


def recall_rows(s, kind):
    return [r for r in csv.DictReader(open(T / s / "tensors/somatic.recall.tsv"), delimiter="\t")
            if (r["kind"] == "SNP") == (kind == "SNV")]


def af05(s, x):
    """SNV cells of set s on its AF 0.05 build (x = the 0.08 cells)."""
    a, b = recall_rows(s, "SNV"), recall_rows(s + "_snvaf05", "SNV")
    nc = lambda rec: sorted(r["truth_id"] for r in rec if r["status"] == "no_candidate")
    assert len(a) == len(b) and nc(a) == nc(b), s
    st = Counter(r["status"] for r in b)
    return dict(x, tensor=sum(st[k] for k in TENSOR), floor=st["filtered"])


def load(s, kind):
    rec = recall_rows(s, kind)
    st = Counter(r["status"] for r in rec)
    other = set(st) - TENSOR - {"filtered", "no_candidate"}
    assert not other, (s, kind, other)
    f = HERE / s / ("snv_no_candidate.tsv" if kind == "SNV" else "indel_no_candidate.tsv")
    if not f.exists():
        return None
    cls = Counter(r["final_class"].split("_")[0] for r in csv.DictReader(open(f), delimiter="\t"))
    assert sum(cls.values()) == st["no_candidate"], (s, kind, sum(cls.values()), st["no_candidate"])
    x = dict(truth=len(rec), tensor=sum(st[k] for k in TENSOR), floor=st["filtered"], none=st["no_candidate"])
    x.update(cls)
    x["other"] = sum(v for k, v in cls.items() if k not in ("I1", "I2", "I3", "I4", "I5"))
    return x


data = {kind: {s: load(s, kind) for s, _ in SETS} for kind in ("SNV", "INDEL")}
for kind, exp in EXPECTED.items():
    for i, (s, _) in enumerate(SETS[:6]):
        for k, v in exp.items():
            got = data[kind][s].get(k, 0)
            assert got == v[i], f"{kind} {s} {k}: {got} != published {v[i]}"
done = [(s, lab) for s, lab in SETS if data["SNV"][s] is not None and data["INDEL"][s] is not None]
f = "{:,}".format


def pct(n, d):
    return f"{f(n)} ({100 * n / d:.1f}%)"


hdr = "| chr1–22, all truth alleles | " + " | ".join(lab for _, lab in done) + " |"
sep = "|---|" + "---|" * len(done)
S = {s: data["SNV"][s] for s, _ in done}
I = {s: data["INDEL"][s] for s, _ in done}
L = ["# Supplementary Tables R1.6b / R1.6c with HCC1395T (generated by r16bc_tables.py)", "",
     "HG008T / COLO829T columns reproduce the published tables exactly (asserted).", "", "## R1.6b SNVs", "", hdr, sep]
L.append("| truth SNVs | " + " | ".join(f(S[s]["truth"]) for s, _ in done) + " |")
L.append("| tensor | " + " | ".join(f(S[s]["tensor"]) for s, _ in done) + " |")
L.append("| candidate below the builder's floors (allele fraction 0.08, three ALT reads) | " + " | ".join(f(S[s]["floor"]) for s, _ in done) + " |")
L.append("| no candidate | " + " | ".join(f(S[s]["none"]) for s, _ in done) + " |")
L.append("| S1 no read with MAPQ > 10 | " + " | ".join(f(S[s].get("S1", 0)) for s, _ in done) + " |")
L.append("| S2 fewer than 3 ALT reads | " + " | ".join(f(S[s].get("S2", 0)) for s, _ in done) + " |")
L.append("| S3 the ALT is a graph SNP allele | " + " | ".join(f(S[s].get("S3", 0)) for s, _ in done) + " |")
L.append("| S4 repeat (a / b / c) | " + " | ".join(
    f"{f(sum(S[s].get(k, 0) for k in ('S4a', 'S4b', 'S4c')))} ({S[s].get('S4a', 0)} / {S[s].get('S4b', 0)} / {S[s].get('S4c', 0)})" for s, _ in done) + " |")
S5 = {s: af05(s, S[s]) if s in AF05 else S[s] for s, _ in done}
L += ["", "## R1.6b SNVs on the current SNV pipeline (AF 0.05 for COLO829T / HCC1395T PacBio HiFi, Fiber-seq, Illumina)", "", hdr, sep]
L.append("| truth SNVs | " + " | ".join(f(S5[s]["truth"]) for s, _ in done) + " |")
L.append("| tensor | " + " | ".join(f(S5[s]["tensor"]) for s, _ in done) + " |")
L.append("| candidate below the builder's floors (allele fraction 0.05, or 0.08 for HG008T and ONT; three ALT reads) | "
         + " | ".join(f(S5[s]["floor"]) for s, _ in done) + " |")
L.append("| no candidate (S1-S4 as above) | " + " | ".join(f(S5[s]["none"]) for s, _ in done) + " |")
L += ["", "## R1.6c INDELs", "", hdr, sep]
L.append("| truth INDELs | " + " | ".join(f(I[s]["truth"]) for s, _ in done) + " |")
L.append("| tensor | " + " | ".join(f(I[s]["tensor"]) for s, _ in done) + " |")
L.append("| candidate below the builder's floors | " + " | ".join(f(I[s]["floor"]) for s, _ in done) + " |")
L.append("| no candidate | " + " | ".join(f(I[s]["none"]) for s, _ in done) + " |")
L.append("| I1 allele fully in the graph (complete absorption) | " + " | ".join(pct(I[s].get("I1", 0), I[s]["none"]) for s, _ in done) + " |")
L.append("| I2 allele partly in the graph (partial absorption) | " + " | ".join(pct(I[s].get("I2", 0), I[s]["none"]) for s, _ in done) + " |")
L.append("| I3 reads on GRCh38, edits spell another allele | " + " | ".join(f(I[s].get("I3", 0)) for s, _ in done) + " |")
L.append("| I4 fewer than 3 ALT reads | " + " | ".join(f(I[s].get("I4", 0)) for s, _ in done) + " |")
L.append("| I5 no read with MAPQ > 10 | " + " | ".join(f(I[s].get("I5", 0)) for s, _ in done) + " |")
L.append("| I6 longer than 50 bp, I7 unclear, I9 other | " + " | ".join(f(I[s]["other"]) for s, _ in done) + " |")
missing = [lab for s, lab in SETS if (s, lab) not in done]
if missing:
    L += ["", "Not finished (skipped): " + ", ".join(missing) + "."]
(HERE / "r16bc_tables.md").write_text("\n".join(L) + "\n")
print("\n".join(L))

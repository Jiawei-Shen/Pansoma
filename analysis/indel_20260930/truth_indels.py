"""chr1 somatic INDEL truths (PASS, in BED): candidate status x repeat context x allele in the d9 graph x composite
(GRCh38 allele != the HG008-N-assembly allele) x PoN tag. Read-only analysis."""
import collections
import csv
import gzip
import json
from pathlib import Path

import pysam

T = "/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors"
TRUTH = "/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz"
FASTA = pysam.FastaFile("/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
D9 = pysam.TabixFile("/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz")
PON = "/scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_INDEL_base/vcf_chr1/eval_bed/truth.INDEL.pon.vcf.gz"
OUT = Path(__file__).with_name("truth_indels.tsv")

status = {}
with open(f"{T}/somatic.recall.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["chrom"] == "chr1" and r["kind"] in ("INS", "DEL"):
            status[(int(r["vcf_pos"]), r["kind"])] = (r["status"], r["detail"], r["passed"] == "True", r["in_bed"] == "True")

pon = {}
for line in gzip.open(PON, "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    info = dict(x.partition("=")[::2] for x in f[7].split(";"))
    pon[(int(f[1]), f[3], f[4])] = info.get("PANSOMA_PON", "")


def apply(seq, start, pos0, ref, alt):
    """seq = reference from `start`; replace ref at pos0 by alt."""
    i = pos0 - start
    assert seq[i:i + len(ref)].upper() == ref.upper(), (pos0, ref, seq[i:i + len(ref)])
    return seq[:i] + alt + seq[i + len(ref):]


def period(s):
    for p in range(1, len(s) + 1):
        if len(s) % p == 0 and s[:p] * (len(s) // p) == s:
            return s[:p]
    return s


def repeat_context(pos0, ref, alt):
    """(unit, copies in the reference after the padding base): the tandem run the event sits in (left-aligned)."""
    ev = alt[1:] if len(alt) > len(ref) else ref[1:]
    unit = period(ev.upper())
    after = FASTA.fetch("chr1", pos0 + 1, pos0 + 1 + 400).upper()
    n = 0
    while after[n * len(unit):(n + 1) * len(unit)] == unit:
        n += 1
    return unit, n


def in_graph(pos0, ref, alt):
    lo, hi = pos0 - 60, pos0 + len(ref) + 60
    for rec in D9.fetch("chr1", max(0, lo), hi):
        g = rec.split("\t")
        gpos0, gref = int(g[1]) - 1, g[3]
        if len(gref) > 5000:
            continue
        s, e = min(pos0, gpos0) - 30, max(pos0 + len(ref), gpos0 + len(gref)) + 30
        seq = FASTA.fetch("chr1", s, e).upper()
        want = apply(seq, s, pos0, ref.upper(), alt.upper())
        for galt in g[4].split(","):
            if galt.startswith("<") or galt == "*":
                continue
            if apply(seq, s, gpos0, gref.upper(), galt.upper()) == want:
                return True
    return False


rows = []
for line in gzip.open(TRUTH, "rt"):
    if line[0] == "#":
        continue
    f = line.rstrip("\n").split("\t")
    if f[0] != "chr1" or f[6] not in ("PASS", "."):
        continue
    pos, ref, alt = int(f[1]), f[3], f[4]
    if "," in alt or len(ref) == len(alt):
        continue
    kind = "INS" if len(alt) > len(ref) else "DEL"
    st = status.get((pos, kind))
    if st is None or not st[3]:
        continue
    info = dict(x.partition("=")[::2] for x in f[7].split(";"))
    n = info.get("HG008Nv62SOMATICVARIANT", "")
    nref, nalt = n.split(":")[1].split("-")[1:3] if n else ("", "")
    unit, copies = repeat_context(pos - 1, ref, alt)
    rows.append(dict(pos=pos, kind=kind, len=len(alt) - len(ref), status=st[0], detail=st[1], unit=unit, copies=copies,
                     hp=len(unit) == 1 and copies >= 5, str_=len(unit) > 1 and copies >= 3,
                     composite=(len(nalt) - len(nref)) != (len(alt) - len(ref)), in_graph=in_graph(pos - 1, ref, alt),
                     pon=pon.get((pos, ref, alt), "?")))

with open(OUT, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
C = collections.Counter
print("chr1 INDEL truths PASS in BED:", len(rows))
print("status", C(r["status"] for r in rows))
print("context", C("homopolymer>=5" if r["hp"] else "STR>=3" if r["str_"] else "other" for r in rows))
print("size", C("1bp" if abs(r["len"]) == 1 else "2-5bp" if abs(r["len"]) <= 5 else ">5bp" for r in rows))
print("in d9 graph", C(r["in_graph"] for r in rows), " composite", C(r["composite"] for r in rows))
print("PoN", C(("CoLoRSdb" in r["pon"], "gnomAD" in r["pon"]) for r in rows))
print("status x in_graph", sorted(C((r["status"], r["in_graph"]) for r in rows).items()))
print("status x context", sorted(C((r["status"], "hp" if r["hp"] else "str" if r["str_"] else "other") for r in rows).items()))
print("status x composite", sorted(C((r["status"], r["composite"]) for r in rows).items()))
print("in_graph x composite", sorted(C((r["in_graph"], r["composite"]) for r in rows).items()))

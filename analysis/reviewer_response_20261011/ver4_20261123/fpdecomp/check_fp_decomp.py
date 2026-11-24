#!/usr/bin/env python3
"""Checks of fp_decomp.json (ver4 false-positive decomposition). Writes checks.txt next to this script.

1. Every cell / arm / stage: FP allele keys == FP records (so shares of alleles are shares of FP records); class sums.
2. Truth (TP + FN) equal in the two arms at both stages (also asserted by fp_decomp.py against report.json and the
   rtg tp-baseline / fn files).
3. Before-filter FP counts and F1 against section 4 of latest_results_20261123.md (18 cells x GRCh38-only / v2.1);
   after-filter P / R / F1 against sections 2a / 2b.
4. Independent recount of one new cell, HCC1395T Illumina SNV (before and after the filter), without calls_lib: plain
   (pos, ref, alt) of single-base records in fp.vcf.gz, the germline VCF (alleles a genotype carries), the carrier
   table key chr1:pos:ref:alt, and the single-base tensors of HCC1395T_Illumina_d46_snvaf05."""
import gzip
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
LATEST = Path("/wanglab/jshen/scratch_jshen/Pansoma_papers/ver1/reviewer_response/latest_results_20261123.md")
V21 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition/v2.1/hprc_alleles/chr1.tsv.gz")
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
out = []
bad = 0


def say(ok, text):
    global bad
    bad += not ok
    out.append(("OK   " if ok else "FAIL ") + text)


J = json.loads((HERE / "fp_decomp.json").read_text())
assert len(J) == 18

# 1 + 2
n1 = n2 = 0
for e in J:
    for st, s in e["stages"].items():
        ok = s["lin_fp_records"] == s["lin_fp_keys"] and s["d46_fp_records"] == s["d46_fp_keys"]
        ok &= s["germline_in_graph"]["n"] + s["germline_not_in_graph"]["n"] + s["non_germline"]["n"] == s["all"]["n"] == s["lin_fp_keys"]
        ok &= s["germline_in_graph"]["n"] + s["germline_not_in_graph"]["n"] == s["germline_all"]["n"]
        n1 += ok
        r = e["report"]
        t = {r[f"{arm}_{st}"]["tp"] + r[f"{arm}_{st}"]["fn"] for arm in ("lin", "d46")} | {r[f"{arm}_{st}"]["truth"] for arm in ("lin", "d46")}
        n2 += len(t) == 1 and r[f"lin_{st}"]["fp"] == s["lin_fp_records"] and r[f"d46_{st}"]["fp"] == s["d46_fp_records"]
say(n1 == 36, f"FP allele keys == FP records in both arms and class sums == total: {n1} of 36 cell x stage")
say(n2 == 36, f"truth (TP + FN) equal in the two arms and equal to report.json baseline; FP == report.json fp: {n2} of 36 cell x stage")
truths = {(e["tumor"], e["kind"]): e["report"]["lin_raw"]["truth"] for e in J}
say(len({(e["tumor"], e["kind"], e["report"][f"{a}_{s}"]["truth"]) for e in J for a in ("lin", "d46") for s in ("raw", "pon")}) == 6,
    f"one truth count per tumor x kind over platforms, arms and stages: {truths}")

# 3
text = LATEST.read_text().split("\n## ")
sec = {s.split("\n", 1)[0]: s for s in text}
s4 = next(v for k, v in sec.items() if k.startswith("4. Six-tumor models, chr1 PASS before"))
s2 = next(v for k, v in sec.items() if k.startswith("2. Six-tumor models on the three references, chr1 after"))
by = {(e["kind"], e["dataset"]): e for e in J}
f1 = lambda x: 2 * x["tp"] / (2 * x["tp"] + x["fp"] + x["fn"])
n = 0
for line in s4.split("\n"):
    c = [x.strip() for x in line.strip("|").split("|")]
    if len(c) != 6 or c[0] not in ("SNV", "INDEL"):
        continue
    e = by[(c[0], f"{c[1]} {c[2]}")]
    for col, arm in ((3, "lin"), (5, "d46")):
        m = re.fullmatch(r"([\d.]+) \(([\d,]+) FP\)", c[col])
        r = e["report"][f"{arm}_raw"]
        ok = int(m.group(2).replace(",", "")) == r["fp"] == e["stages"]["raw"][f"{arm}_fp_records"] and abs(float(m.group(1)) - f1(r)) < 0.00101  # the file prints report.json's rounded f1
        n += ok
        if not ok:
            say(False, f"latest_results section 4 {c[0]} {c[1]} {c[2]} {arm}: {c[col]} vs FP {r['fp']} F1 {f1(r):.4f}")
say(n == 36, f"before-filter FP count (exact) and F1 (within 0.001) of both arms equal to latest_results_20261123.md section 4: {n} of 36")
n = 0
for kind, sub in (("SNV", "### 2a. SNV"), ("INDEL", "### 2b. INDEL")):
    block = s2.split(sub)[1].split("\n### ")[0]
    for line in block.split("\n"):
        c = [x.strip() for x in line.strip("|").split("|")]
        if len(c) < 6 or (kind, f"{c[0]} {c[1]}") not in by:
            continue
        e = by[(kind, f"{c[0]} {c[1]}")]
        ok = int(c[2].replace(",", "")) == e["report"]["lin_pon"]["truth"]
        for col, arm in ((3, "lin"), (5, "d46")):
            p, r, f = (float(x.replace("*", "")) for x in c[col].split(" / "))
            a = e["report"][f"{arm}_pon"]
            ok &= abs(p - a["tp"] / (a["tp"] + a["fp"])) < 0.00101 and abs(r - a["tp"] / a["truth"]) < 0.00101 and abs(f - f1(a)) < 0.00101
        n += ok
        if not ok:
            say(False, f"latest_results section 2 {kind} {c[0]} {c[1]}: {c[2:6]} vs {e['report']['lin_pon']} {e['report']['d46_pon']}")
say(n == 18, f"after-filter truth count and P / R / F1 of both arms equal to latest_results_20261123.md sections 2a / 2b (within rounding): {n} of 18 rows")

# 4
e = by[("SNV", "HCC1395T Illumina")]


def simple(path):
    keys, n_rec, n_other = set(), 0, 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 6)
            if c[0] != "chr1":
                continue
            n_rec += 1
            alts = c[4].split(",")
            if len(c[3]) == 1 and all(len(a) == 1 and a in "ACGT" for a in alts):
                keys |= {(int(c[1]), c[3].upper(), a.upper()) for a in alts}
            else:
                n_other += 1
    return keys, n_rec, n_other


germ = set()
seen = False
with gzip.open(e["germline_vcf"], "rt") as f:
    for line in f:
        if line[0] == "#":
            continue
        c = line.rstrip("\n").split("\t")
        if c[0] != "chr1":
            if seen:
                break
            continue
        seen = True
        gt = dict(zip(c[8].split(":"), c[9].split(":"))).get("GT", ".")
        alts = c[4].split(",")
        for g in re.split(r"[|/]", gt):
            if g.isdigit() and int(g) > 0:
                a = alts[int(g) - 1].upper()
                r = c[3].upper()
                if len(r) == len(a):  # SNV or MNP: base by base
                    germ |= {(int(c[1]) + i, r[i], a[i]) for i in range(len(r)) if r[i] != a[i]}
fps = {}
for st, which in (("raw", "raw_calls"), ("pon", "pon_calls")):
    fps[st] = simple(Path(e["lin_dir"]) / "eval_bed" / "vcfeval" / which / "fp.vcf.gz")
need = fps["raw"][0] | fps["pon"][0]
car = {}
with gzip.open(V21, "rt") as f:
    f.readline()
    for line in f:
        c = line.split("\t", 7)
        k = c[0].split(":")
        if len(k) == 4 and len(k[2]) == 1 and len(k[3]) == 1:
            key = (int(k[1]), k[2], k[3])
            if key in need:
                car[key] = int(c[6])
tens = set()
for line in open(TENSORS / e["d46_set"] / "tensors" / "SNV" / "chr1_labels.ndjson"):
    g = json.loads(line).get("grch38")
    if g and g["chrom"] == "chr1" and len(g["ref"]) == 1 and len(g["alt"]) == 1:
        k = (g["pos0"] + 1, g["ref"].upper(), g["alt"].upper())
        if k in need:
            tens.add(k)
for st in ("raw", "pon"):
    keys, n_rec, n_other = fps[st]
    a = [k for k in keys if k in germ and car.get(k, 0) >= 46]
    b = [k for k in keys if k in germ and car.get(k, 0) < 46]
    o = [k for k in keys if k not in germ]
    mine = dict(records=n_rec, in_graph=(len(a), sum(k in tens for k in a)), not_in_graph=(len(b), sum(k in tens for k in b)),
                non_germline=(len(o), sum(k in tens for k in o)))
    s = e["stages"][st]
    theirs = dict(records=s["lin_fp_records"], in_graph=(s["germline_in_graph"]["n"], s["germline_in_graph"]["d46_tensor"]),
                  not_in_graph=(s["germline_not_in_graph"]["n"], s["germline_not_in_graph"]["d46_tensor"]),
                  non_germline=(s["non_germline"]["n"], s["non_germline"]["d46_tensor"]))
    say(mine == theirs, f"independent recount HCC1395T Illumina SNV, {st} (records; (n, candidate on v2.1) per class; {n_other} non-single-base FP records): "
                        f"recount {mine} vs fp_decomp.json {theirs}")

out.append(f"{'ALL CHECKS PASSED' if not bad else str(bad) + ' CHECK(S) FAILED'}")
(HERE / "checks.txt").write_text("\n".join(out) + "\n")
print("\n".join(out))

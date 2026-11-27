#!/usr/bin/env python3
"""strata.md: Response Table 2.7 (HG008T chr1 SNV by the number of nearby germline variants) with the HPRC rows on the
six-tumor models of this batch, so that no model of the table takes the population site features as an input
(user 2026-11-24).

The stratification is E19's (../reviewer_response_20261011/experiments/E19_germline_density/strata_latest.py, the
source of the table): the rtg records of the PASS calls after the population filter in the benchmark BED
(<vcf_dir>/eval_bed/vcfeval/pon_calls; tp-baseline + fn = truth, tp + fp = calls), each counted by the germline
records of HG008-N (dipcall, any FILTER) that intersect its normalized interval widened by 50 bp, the records that
intersect the variant itself not counted; strata 0 / 1 / >=2; recall by the truth variant's stratum, precision by the
call's. Same functions (calls_lib: Germline.count_near, normalize, interval, prf); only the call sets differ.

  HPRC v1.1 / v2.1  $PANSOMA_RUNS/m6_nositef_graphs_20261118/multicancer6_<P>_SNV_{d9,d46}_nositef/vcf_chr1_<SET>
                    (stratified here)
  GRCh38-only       m4m6_graphs_20261110/multicancer6_<P>_SNV_linear_nositef: its call sets are on the old cluster,
                    so its rows are those of strata_latest.md section 7. For the F1 difference "before rounding" the
                    counts of a row are recovered from its printed values (truth x recall = TP; TP calls = TP, which
                    must reproduce the printed precision and F1 to three decimals; asserted).

Checks (asserted): every call set's totals equal its report.json; the truth strata equal those of strata_latest.md
(the same truth after the population filter in every row); the counts recovered from strata_latest.md's rows with the
site features reproduce its printed F1 differences (the worst deviation is printed).

The germline VCF and the FASTA of labels.manifest.json are old-cluster paths; their storage2 copies are read (the
VCF's sha256 is checked against the manifest). calls_lib.chrom_seq needs pysam, replaced here by a .fai read.
usage: /home/jiawei.s/miniconda3/envs/convnextv2/bin/python strata.py   (login node; needs numpy)
"""
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
E = HERE.parent / "reviewer_response_20261011/experiments"
sys.path.insert(0, str(E / "call_level_common"))
sys.path.insert(0, str(HERE))
import calls_lib as L  # noqa: E402
import collect as C  # noqa: E402

DATA = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset")
TENSORS = DATA / "pansoma_v2_tensors"
FASTA = DATA / "REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
OLD = "/scratch/jshen/data/"  # manifest prefix of the truth files (ris/submit.sh mounts the storage2 copies there)
LATEST = E / "E19_germline_density/strata_latest.md"
TUMOR, W = "HG008T", 50
PARTS = ("tp-baseline", "fn", "tp", "fp")
SL = ("0", "1", "2+")
_SEQ = {}


def chrom_seq(fasta, chrom=L.CHROM):
    key = (str(fasta), chrom)
    if key not in _SEQ:
        for line in open(str(fasta) + ".fai"):
            name, length, offset, bases, width = line.split("\t")[:5]
            if name == chrom:
                break
        else:
            sys.exit(f"{chrom} not in {fasta}.fai")
        length, bases, width = int(length), int(bases), int(width)
        with open(fasta, "rb") as f:
            f.seek(int(offset))
            _SEQ[key] = f.read(length + (length // bases + 1) * (width - bases)).replace(b"\n", b"").decode()[:length].upper()
        assert len(_SEQ[key]) == length, (fasta, chrom)
    return _SEQ[key]


L.chrom_seq = chrom_seq


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def stratify(vcf_dir, germ, seq):
    """{stratum: {part: records}} of one call set, and its report.json."""
    n = {s: dict.fromkeys(PARTS, 0) for s in SL + ("all",)}
    for part in PARTS:
        for rec in L.read_records(vcf_dir / "eval_bed" / "vcfeval" / "pon_calls" / f"{part}.vcf.gz"):
            key = next((k for k in (L.normalize(seq, rec.pos, rec.ref, alt) for alt in rec.alts) if k), None)
            if key is None:
                continue
            a, b = L.interval(*key)
            k = germ.count_near(a, b, W)
            n["0" if k == 0 else "1" if k == 1 else "2+"][part] += 1
            n["all"][part] += 1
    return n, json.loads((vcf_dir / "eval_bed" / "report.json").read_text())


def row(x):
    """truth, fp, P, R, F1 of one stratum's counts."""
    p, r, f = L.prf(x["tp-baseline"], x["fn"], x["tp"], x["fp"])
    return dict(truth=x["tp-baseline"] + x["fn"], fp=x["fp"], p=p, r=r, f=f)


def latest():
    """Section 7 of strata_latest.md: truth strata and {(platform, reference): {stratum: printed cells}}."""
    lines = LATEST.read_text().splitlines()
    i = next(i for i, x in enumerate(lines) if x.startswith(f"## 7. chr1 SNV, {TUMOR}, six-tumor model"))
    truth = [int(v) for v in re.search(r"(\d+) / (\d+) / (\d+)\)", lines[i]).groups()]
    out, plat = {}, None
    for line in lines[i + 4:]:
        if not line.startswith("|"):
            break
        c = [x.strip() for x in line.strip("|").split("|")]
        plat = c[0] or plat
        tri = [[v.replace("–", "").strip() for v in x.split(" / ")] for x in c[2:7]]
        out[plat, c[1]] = {s: dict(truth=truth[j], fp=int(tri[0][j]), p=float(tri[1][j]), r=float(tri[2][j]), f=float(tri[3][j]),
                                   d=float(tri[4][j]) if len(tri[4]) == 3 and tri[4][j] else None) for j, s in enumerate(SL)}
    return dict(zip(SL, truth)), out


def recover(c):
    """Exact F1 of a printed row of strata_latest.md (TP = truth x recall, TP calls = TP)."""
    tp = round(c["r"] * c["truth"])
    p, r, f = L.prf(tp, c["truth"] - tp, tp, c["fp"])
    assert (round(p, 3), round(r, 3), round(f, 3)) == (c["p"], c["r"], c["f"]), c
    return f


def f3(x, sign=False):
    return "–" if x is None else f"{x:+.3f}" if sign else f"{x:.3f}"


def tri(vals, sign=False):
    return " / ".join(str(v) if isinstance(v, int) else f3(v, sign) for v in vals)


def main():
    truth, old = latest()
    base = {k: {s: recover(c[s]) for s in SL} for k, c in old.items()}
    dev = max(abs(base[k][s] - base[k[0], "GRCh38-only"][s] - c[s]["d"]) for k, c in old.items() for s in SL if c[s]["d"] is not None)
    assert dev < 0.0011, dev  # the printed differences are rounded once more
    res, germ_vcf = {}, None
    for p in C.PLATFORMS:
        for g, gl in C.GRAPHS:
            st = C.test_set(TUMOR, p, g)
            m = json.loads((TENSORS / st / "tensors" / "SNV" / "labels.manifest.json").read_text())["truth"]["germline"]
            assert m["vcf"].startswith(OLD), m["vcf"]
            vcf = DATA / m["vcf"][len(OLD):]
            if germ_vcf is None:
                assert sha256(vcf) == m["vcf_sha256"], vcf
            assert germ_vcf in (None, vcf)
            germ_vcf = vcf
            n, rep = stratify(C.RUNS / f"multicancer6_{p}_SNV_{g}_nositef" / f"vcf_chr1_{st}", L.germline(str(vcf), FASTA), chrom_seq(FASTA))
            a, t = rep["pon"]["at_threshold"], n["all"]
            assert (t["tp-baseline"], t["fn"], t["tp"], t["fp"]) == (a["tp_baseline"], a["fn"], a["tp_call"], a["fp"]), st
            res[C.plat(TUMOR, p), gl] = {s: row(n[s]) for s in SL + ("all",)}
            assert all(res[C.plat(TUMOR, p), gl][s]["truth"] == truth[s] for s in SL), st
    out = [f"# Response Table 2.7 without the population site features in any model ({TUMOR} chr1 SNV)", "",
           "Generated by `strata.py`. chr1, benchmark regions, PASS calls after the population filter. Stratum = germline "
           "variants of HG008-N within 50 bp of the variant (0 / 1 / ≥2); recall by the truth variant's stratum, precision "
           "by the call's stratum. Six-tumor models, seed 0, none with the population site features as a model input: "
           "GRCh38-only = `multicancer6_<P>_SNV_linear_nositef` (rows of `E19_germline_density/strata_latest.md` section 7, "
           "unchanged), HPRC v1.1 / v2.1 = `m6_nositef_graphs_20261118/multicancer6_<P>_SNV_{d9,d46}_nositef` (stratified "
           "here). ΔF1 = HPRC − GRCh38-only of the same platform, from the unrounded F1. Truth SNVs 0 / 1 / ≥2: "
           f"{tri(list(truth.values()))}.", "",
           "| Platform | Reference | False positives (0 / 1 / ≥2) | Precision (0 / 1 / ≥2) | Recall (0 / 1 / ≥2) | F1 (0 / 1 / ≥2) | ΔF1, HPRC − GRCh38-only (0 / 1 / ≥2) |",
           "|---|---|---|---|---|---|---|"]
    delta = {}
    for p in C.PLATFORMS:
        pl = C.plat(TUMOR, p)
        o = old[pl, "GRCh38-only"]
        out.append(f"| {pl} | GRCh38-only | {tri([o[s]['fp'] for s in SL])} | {tri([o[s]['p'] for s in SL])} | "
                   f"{tri([o[s]['r'] for s in SL])} | {tri([o[s]['f'] for s in SL])} | – |")
        for _, gl in C.GRAPHS:
            S = res[pl, gl]
            delta[pl, gl] = [S[s]["f"] - base[pl, "GRCh38-only"][s] for s in SL]
            out.append(f"|  | {gl} | {tri([S[s]['fp'] for s in SL])} | {tri([S[s]['p'] for s in SL])} | "
                       f"{tri([S[s]['r'] for s in SL])} | {tri([S[s]['f'] for s in SL])} | {tri(delta[pl, gl], True)} |")
    out += ["", "## Against the HPRC models with the site features (the table until 2026-11-24)", "",
            "| Platform | Reference | ΔF1 with the site features (0 / 1 / ≥2) | ΔF1 without (0 / 1 / ≥2) | whole chr1 P / R / F1 without |",
            "|---|---|---|---|---|"]
    for (pl, gl), d in delta.items():
        a = res[pl, gl]["all"]
        out.append(f"| {pl} | {gl} | {tri([old[pl, gl][s]['d'] for s in SL], True)} | {tri(d, True)} | {tri([a['p'], a['r'], a['f']])} |")
    lo = lambda j: (min(d[j] for d in delta.values()), max(d[j] for d in delta.values()))  # noqa: E731
    near = [x for d in delta.values() for x in d[1:]]
    out += ["", "## Ranges of ΔF1 over the six HPRC rows", "",
            f"- no nearby germline variant: {f3(lo(0)[0], True)} to {f3(lo(0)[1], True)}",
            f"- one: {f3(lo(1)[0], True)} to {f3(lo(1)[1], True)}",
            f"- at least two: {f3(lo(2)[0], True)} to {f3(lo(2)[1], True)}",
            f"- one or at least two: {f3(min(near), True)} to {f3(max(near), True)}",
            f"- rows whose ΔF1 with one nearby variant exceeds that with none: {sum(d[1] > d[0] for d in delta.values())} of 6; "
            f"with at least two: {sum(d[2] > d[0] for d in delta.values())} of 6", "",
            "## Checks", "",
            "- totals of every call set equal its `eval_bed/report.json` (asserted)",
            f"- truth strata of every call set equal strata_latest.md's ({tri(list(truth.values()))}; asserted)",
            f"- germline VCF: `{germ_vcf}` (sha256 equal to the manifests')",
            f"- GRCh38-only F1 recovered from the printed rows: strata_latest.md's own ΔF1 column reproduced within {dev:.4f}"]
    (HERE / "strata.md").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()

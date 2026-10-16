#!/usr/bin/env python3
"""E18 (R3.5): SBS96 spectra of the chr1 truth SNVs and of Pansoma's PASS-after-PoN true and false positives.

    /wanglab/jshen/anaconda3/bin/python e18_sbs96.py

Input: the SNV units of ../E16_bootstrap/callsets.py (rtg eval_bed/vcfeval/pon_calls: tp = call-side TPs, fp, fn,
tp-baseline), comparators.tsv rows of kind SNV (slot, none yet), and per set the labels manifest (truth, BED, germline
truth, FASTA).
Spectra (single-base substitutions on chr1 inside the set's benchmark BED; trinucleotide context from the labels'
GRCh38 FASTA, pyrimidine-centred: a G or A reference base is reverse-complemented with its context):
  truth     tp-baseline + fn of the unit (= every truth SNV in the BED; identical across the units of one truth / BED)
  TP / FP   the unit's PASS calls after the PoN that rtg matched / did not match
  FN        truth SNVs missed after the PoN (recall per substitution class = TPb / (TPb + FN))
  germline  the set's germline truth (dipcall of the normal, or the CASTLE / SEQC2 germline VCF) on chr1 in the same BED,
            FILTER PASS / . / GAP1 / GAP2 (the labels' germline rule), single-base ALTs carried by the genotype
Signature references: COSMIC v3.4 SBS GRCh38 from the installed SigProfilerAssignment package (local, not downloaded).
Hallmarks: UV = C>T at a 5' pyrimidine (T[C>T]N + C[C>T]N), share of all SNVs; CpG = N[C>T]G share (SBS1-like);
CC>TT doublets = adjacent SNV pairs (POS, POS+1) with REF CC -> TT or GG -> AA (not in SBS96; the truth VCFs on chr1
have no MNP records).
Writes sbs96.tsv (all spectra, COSMIC channel order), summary.json and tables.md next to this script; the germline
spectra are cached in tmp/reviewer_cpu_20261013/E18_sbs96/germline_cache.json.
"""
import csv
import glob
import gzip
import json
import math
import sys
from bisect import bisect_right
from collections import Counter
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "E16_bootstrap"))
import callsets as C  # noqa: E402

CACHE = Path("/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/E18_sbs96/germline_cache.json")
COSMIC_GLOB = "/wanglab/jshen/anaconda3/lib/python3*/site-packages/SigProfilerAssignment/data/Reference_Signatures/GRCh38/COSMIC_v3.4_SBS_GRCh38.txt"
NAMED = ("SBS1", "SBS5", "SBS7a", "SBS7b", "SBS18", "SBS40a")
COMP = str.maketrans("ACGT", "TGCA")
CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
CHROM = "chr1"


def load_cosmic():
    paths = sorted(glob.glob(COSMIC_GLOB))
    if not paths:
        return None, None, None
    with open(paths[0]) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    channels = [r[0] for r in rows]
    sigs = {name: {r[0]: float(r[i]) for r in rows} for i, name in enumerate(header) if i > 0}
    return str(Path(paths[0]).resolve()), channels, sigs


def default_channels():
    return [f"{a}[{s}]{b}" for s in CLASSES for a in "ACGT" for b in "ACGT"]


class Genome:
    def __init__(self, fasta):
        import pysam
        self.seq = pysam.FastaFile(fasta).fetch(CHROM).upper()

    def key(self, pos, ref, alt):
        """SBS96 channel of a single-base substitution at 1-based pos, or None."""
        ctx = self.seq[pos - 2:pos + 1]
        if len(ctx) != 3 or "N" in ctx or ctx[1] != ref or alt not in "ACGT" or alt == ref:
            return None
        if ref in "GA":
            ctx, ref, alt = ctx.translate(COMP)[::-1], ref.translate(COMP), alt.translate(COMP)
        return f"{ctx[0]}[{ref}>{alt}]{ctx[2]}"


def snvs(records):
    """(pos, ref, alt) single-base substitutions on chr1."""
    return [(r[1], r[2], r[3]) for r in records if r[0] == CHROM and len(r[2]) == 1 and len(r[3]) == 1]


def spectrum(genome, items):
    c, skipped = Counter(), 0
    for pos, ref, alt in items:
        k = genome.key(pos, ref, alt)
        if k:
            c[k] += 1
        else:
            skipped += 1
    return c, skipped


def doublets(items):
    """(adjacent SNV pairs, CC>TT / GG>AA pairs) among (pos, ref, alt)."""
    by = {p: (r, a) for p, r, a in items}
    pairs = cctt = 0
    for p, (r, a) in by.items():
        if p + 1 in by:
            pairs += 1
            r2, a2 = by[p + 1]
            if (r + r2, a + a2) in (("CC", "TT"), ("GG", "AA")):
                cctt += 1
    return pairs, cctt


def cosine(a, b):
    keys = set(a) | set(b)
    dot = sum(a.get(k, 0) * b.get(k, 0) for k in keys)
    na, nb = math.sqrt(sum(v * v for v in a.values())), math.sqrt(sum(v * v for v in b.values()))
    return dot / (na * nb) if na and nb else None


def mix_truth_germline(c, truth, germline, channels):
    """Non-negative least squares of a spectrum (proportions) on the truth and germline spectra (proportions):
    (truth-like share of the fit, cosine of the fit); None when a spectrum is empty."""
    import numpy as np
    from scipy.optimize import nnls
    vec = lambda x: np.array([x.get(k, 0) for k in channels], dtype=float)  # noqa: E731
    y, t, g = vec(c), vec(truth), vec(germline)
    if not y.sum() or not t.sum() or not g.sum():
        return None
    y, t, g = y / y.sum(), t / t.sum(), g / g.sum()
    w, _ = nnls(np.column_stack([t, g]), y)
    fit = w[0] * t + w[1] * g
    cos = float(fit @ y / (np.linalg.norm(fit) * np.linalg.norm(y))) if fit.any() else None
    return (float(w[0] / w.sum()) if w.sum() else None, cos)


def hallmarks(c):
    n = sum(c.values())
    if not n:
        return dict(n=0)
    cls = Counter()
    for k, v in c.items():
        cls[k[2:5]] += v
    uv = sum(v for k, v in c.items() if k[2:5] == "C>T" and k[0] in "CT")
    cpg = sum(v for k, v in c.items() if k[2:5] == "C>T" and k[6] == "G")
    return dict(n=n, uv=uv / n, cpg=cpg / n, **{f"frac_{x}": cls[x] / n for x in CLASSES})


def read_bed(bed):
    iv = []
    with open(bed) as f:
        for line in f:
            c = line.split("\t")
            if c[0] == CHROM:
                iv.append((int(c[1]), int(c[2])))
    iv.sort()
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [m[0] for m in merged], [m[1] for m in merged]


def in_bed(bed, pos):
    starts, ends = bed
    i = bisect_right(starts, pos - 1) - 1
    return i >= 0 and pos - 1 < ends[i]


def germline_items(vcf, bed):
    """Germline single-base SNVs (pos, ref, alt) on chr1 in the BED: FILTER PASS / . / GAP1 / GAP2, ALT in the GT."""
    import pysam
    b = read_bed(bed)
    out = []

    def take(line):
        c = line.rstrip("\n").split("\t")
        if c[0] != CHROM:
            return
        if not (c[6] in ("PASS", ".") or set(c[6].split(";")) <= {"GAP1", "GAP2"}):
            return
        gt = c[9].split(":")[0].replace("|", "/").split("/") if len(c) > 9 else ["1"]
        alts = c[4].split(",")
        pos, ref = int(c[1]), c[3].upper()
        for i, alt in enumerate(alts, 1):
            if str(i) in gt and len(ref) == 1 and len(alt) == 1 and in_bed(b, pos):
                out.append((pos, ref, alt.upper()))
    if Path(str(vcf) + ".tbi").exists():
        for line in pysam.TabixFile(str(vcf)).fetch(CHROM):
            take(line)
    else:  # sorted VCF without an index: stream until chr1 ends
        seen = False
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                on = line.startswith(CHROM + "\t")
                if on:
                    seen = True
                    take(line)
                elif seen:
                    break
    return out


def germline_spectrum(genome, vcf, bed, cache):
    st = Path(vcf).stat()
    key = f"{vcf}|{st.st_size}|{int(st.st_mtime)}|{bed}"
    if key not in cache:
        items = germline_items(vcf, bed)
        c, skipped = spectrum(genome, items)
        cache[key] = dict(counts=dict(c), skipped=skipped, doublets=doublets(items))
    e = cache[key]
    return Counter(e["counts"]), e["skipped"], tuple(e["doublets"])


def main():
    units = [u for u in C.discover(kinds=("SNV",))]
    if not units:
        raise SystemExit("no finished SNV call sets")
    cosmic_path, channels, sigs = load_cosmic()
    channels = channels or default_channels()
    genome = Genome(C.fasta_path(units))
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    spectra = []          # rows for sbs96.tsv
    truth_by = {}         # (truth, bed) -> dict
    unit_rows = []
    for u in units:
        w = u.rtg["bed"]
        rec = {k: snvs(C.read_records(w / f"{k}.vcf.gz")) for k in C.RTG_FILES}
        tkey = (u.truth, u.bed)
        if tkey not in truth_by:
            items = rec["tp-baseline"] + rec["fn"]
            c, sk = spectrum(genome, items)
            m = C.manifest(u.set, "SNV")
            g = None
            if m and m["truth"].get("germline", {}).get("vcf"):
                gc, gsk, gd = germline_spectrum(genome, m["truth"]["germline"]["vcf"], u.bed, cache)
                g = dict(counts=gc, skipped=gsk, doublets=gd, vcf=m["truth"]["germline"]["vcf"])
            truth_by[tkey] = dict(sample=u.sample, set=u.set, counts=c, skipped=sk, n_records=len(items),
                                  doublets=doublets(items), germline=g, units=[])
            spectra.append((f"truth {u.sample}", u.sample, "-", "truth", c))
            if g:
                spectra.append((f"germline {u.sample}", u.sample, "-", "germline", g["counts"]))
        t = truth_by[tkey]
        t["units"].append(u.uid)
        n_truth_unit = len(rec["tp-baseline"]) + len(rec["fn"])
        sp = {}
        for name, items in (("TP", rec["tp"]), ("FP", rec["fp"]), ("FN", rec["fn"]), ("TPb", rec["tp-baseline"])):
            sp[name], _ = spectrum(genome, items)
        for name in ("TP", "FP", "FN"):
            spectra.append((f"{name} {u.label}", u.sample, u.uid, name, sp[name]))
        allc = sp["TP"] + sp["FP"]
        row = dict(uid=u.uid, label=u.label, role=u.role, variant=u.variant, sample=u.sample, set=u.set,
                   n_truth=n_truth_unit, truth_consistent=n_truth_unit == t["n_records"],
                   n_tp=sum(sp["TP"].values()), n_fp=sum(sp["FP"].values()), n_fn=sum(sp["FN"].values()),
                   cos_tp_truth=cosine(sp["TP"], t["counts"]), cos_fp_truth=cosine(sp["FP"], t["counts"]),
                   cos_all_truth=cosine(allc, t["counts"]), cos_fn_truth=cosine(sp["FN"], t["counts"]),
                   cos_fp_germline=cosine(sp["FP"], t["germline"]["counts"]) if t["germline"] else None,
                   cos_tp_germline=cosine(sp["TP"], t["germline"]["counts"]) if t["germline"] else None,
                   tp=hallmarks(sp["TP"]), fp=hallmarks(sp["FP"]), fn=hallmarks(sp["FN"]),
                   tp_doublets=doublets(rec["tp"]), fp_doublets=doublets(rec["fp"]))
        rec_cls = {}
        tagged = pon_tagged(u)
        for x in CLASSES:
            tpb = sum(v for k, v in sp["TPb"].items() if k[2:5] == x)
            fn = sum(v for k, v in sp["FN"].items() if k[2:5] == x)
            n_tag = None
            if tagged is not None:
                n_tag = sum(1 for pos, ref, alt in rec["tp-baseline"] + rec["fn"]
                            if (pos, ref, alt) in tagged and (genome.key(pos, ref, alt) or "     ")[2:5] == x)
            rec_cls[x] = (tpb, fn, tpb / (tpb + fn) if tpb + fn else None, n_tag)
        row["recall_by_class"] = rec_cls
        if t["germline"]:
            row["mix_fp"] = mix_truth_germline(sp["FP"], t["counts"], t["germline"]["counts"], channels)
            row["mix_fn"] = mix_truth_germline(sp["FN"], t["counts"], t["germline"]["counts"], channels)
        uv_tpb = sum(v for k, v in sp["TPb"].items() if k[2:5] == "C>T" and k[0] in "CT")
        uv_fn = sum(v for k, v in sp["FN"].items() if k[2:5] == "C>T" and k[0] in "CT")
        row["recall_uv"] = (uv_tpb, uv_fn, uv_tpb / (uv_tpb + uv_fn) if uv_tpb + uv_fn else None)
        if sigs:
            row["cos_named"] = {s: dict(TP=cosine(sp["TP"], sigs[s]), FP=cosine(sp["FP"], sigs[s])) for s in NAMED
                                if s in sigs}
            row["best_fp"] = best_sigs(sp["FP"], sigs)
            row["best_tp"] = best_sigs(sp["TP"], sigs)
        unit_rows.append(row)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache) + "\n")
    truth_rows = []
    for (truth, bed), t in truth_by.items():
        r = dict(sample=t["sample"], truth=truth, bed=bed, n_records=t["n_records"], skipped=t["skipped"],
                 doublets=t["doublets"], hallmarks=hallmarks(t["counts"]), units=t["units"])
        if t["germline"]:
            g = t["germline"]
            r.update(germline_vcf=g["vcf"], germline_hallmarks=hallmarks(g["counts"]), germline_doublets=g["doublets"],
                     cos_truth_germline=cosine(t["counts"], g["counts"]))
        if sigs:
            r["cos_named"] = {s: cosine(t["counts"], sigs[s]) for s in NAMED if s in sigs}
            r["best"] = best_sigs(t["counts"], sigs)
            if t["germline"]:
                r["germline_best"] = best_sigs(t["germline"]["counts"], sigs)
        truth_rows.append(r)
    with open(HERE / "sbs96.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["spectrum", "sample", "unit", "kind", "n", *channels])
        for name, sample, uid, kind, c in spectra:
            w.writerow([name, sample, uid, kind, sum(c.values()), *[c.get(k, 0) for k in channels]])
    summary = dict(generated=datetime.now().isoformat(timespec="seconds"), cosmic=cosmic_path, named=NAMED,
                   truth=truth_rows, units=unit_rows)
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    (HERE / "tables.md").write_text(tables(truth_rows, unit_rows, cosmic_path) + "\n")
    print(f"{len(unit_rows)} SNV call sets, {len(truth_rows)} truth sets; COSMIC {cosmic_path}")


def pon_tagged(u):
    """{(pos, ref, alt)} of the chr1 truth SNVs the PoN tags (eval_bed/truth.SNV.pon.vcf.gz), or None."""
    if not u.truth_pon:
        return None
    return {(r[1], r[2], r[3]) for r in C.read_records(u.truth_pon)
            if r[0] == CHROM and "PanelOfNormals" in r[4].split(";")}


def best_sigs(c, sigs, k=3):
    if not sum(c.values()):
        return []
    scored = sorted(((cosine(c, s), name) for name, s in sigs.items()), reverse=True)
    return [(name, round(v, 3)) for v, name in scored[:k]]


def mix(m):
    return "-" if not m or m[0] is None else f"{100 * m[0]:.0f} % ({m[1]:.3f})"


def pct(x):
    return "-" if x is None else f"{100 * x:.1f} %"


def tables(truth_rows, unit_rows, cosmic_path):
    f = C.fmt
    lines = ["# E18: SBS96 spectra (chr1, benchmark BED; Pansoma PASS after the PoN)", "",
             f"Generated by `e18_sbs96.py` ({datetime.now():%Y-%m-%d %H:%M}). Trinucleotide context from the labels' "
             "GRCh38 FASTA, pyrimidine-centred. COSMIC: " + (f"`{cosmic_path}` (v3.4 SBS GRCh38, local)" if cosmic_path
                                                            else "no local file; hallmarks only") + ". UV = C>T at TpC "
             "or CpC (share of all SNVs); CpG = C>T at CpG; CC>TT = adjacent SNV pairs CC>TT / GG>AA (DBS, not in "
             "SBS96).", "",
             "## 1. Truth and germline spectra per sample (chr1, in the BED)", "",
             "| sample | truth SNVs | C>T | UV (TpC/CpC C>T) | CpG C>T | CC>TT doublets (adjacent pairs) | best COSMIC "
             "(cosine) | " + " | ".join(NAMED) + " | germline SNVs | germline UV / CpG | cos(truth, germline) | germline "
             "best COSMIC |", "|---|---:|---:|---:|---:|---|---|" + "---:|" * len(NAMED) + "---:|---|---:|---|"]
    for r in truth_rows:
        h = r["hallmarks"]
        gh = r.get("germline_hallmarks") or {}
        best = ", ".join(f"{n} {v:.3f}" for n, v in r.get("best", []))
        gbest = ", ".join(f"{n} {v:.3f}" for n, v in r.get("germline_best", []))
        named = " | ".join(f(r.get("cos_named", {}).get(s)) for s in NAMED)
        lines.append(f"| {r['sample']} | {h['n']} | {pct(h.get('frac_C>T'))} | {pct(h.get('uv'))} | {pct(h.get('cpg'))} | "
                     f"{r['doublets'][1]} ({r['doublets'][0]}) | {best} | {named} | {gh.get('n', '-')} | "
                     f"{pct(gh.get('uv'))} / {pct(gh.get('cpg'))} | {f(r.get('cos_truth_germline'))} | {gbest} |")
    lines += ["", "## 2. Pansoma PASS-after-PoN calls vs the truth spectrum", "",
              "Cosine to the same sample's truth spectrum (TP, FP, TP + FP, FN) and to its germline spectrum (FP). "
              "TP is a subset of the truth, so its cosine is high by construction; the informative columns are FP "
              "(do false calls look like the tumour's signature or like germline / artefact?) and FN (are misses "
              "context-biased?).", "",
              "'FP truth-like share': non-negative least squares of the FP spectrum (proportions) on the sample's truth "
              "and germline spectra; the truth weight's share (fit cosine in brackets). An estimate of how much of "
              "the FP spectrum looks like the tumour's own mutations rather than like germline variants; separable "
              "only where the two spectra differ (cos(truth, germline) in section 1).", "",
              "| call set | role | TP | FP | cos TP | cos FP | cos TP+FP | cos FN | cos FP-germline | FP truth-like share "
              "(fit cos) | UV share TP / FP | CpG share TP / FP | CC>TT pairs TP / FP | best COSMIC of FP |",
              "|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|---|"]
    for r in unit_rows:
        best = ", ".join(f"{n} {v:.3f}" for n, v in r.get("best_fp", []))
        lines.append(f"| {r['label']} | {r['role']} | {r['n_tp']} | {r['n_fp']} | {f(r['cos_tp_truth'])} | "
                     f"{f(r['cos_fp_truth'])} | {f(r['cos_all_truth'])} | {f(r['cos_fn_truth'])} | "
                     f"{f(r['cos_fp_germline'])} | {mix(r.get('mix_fp'))} | "
                     f"{pct(r['tp'].get('uv'))} / {pct(r['fp'].get('uv'))} | "
                     f"{pct(r['tp'].get('cpg'))} / {pct(r['fp'].get('cpg'))} | {r['tp_doublets'][1]} / "
                     f"{r['fp_doublets'][1]} | {best} |")
    lines += ["", "## 3. Recall after the PoN by substitution class (truth side: TPb / (TPb + FN))", "",
              "Cell: recall (truth SNVs of the class in the BED; of them tagged by the PoN, i.e. lost to any caller "
              "after the PoN).", "",
              "| call set | " + " | ".join(CLASSES) + " | UV C>T (TpC/CpC) |", "|---|" + "---:|" * (len(CLASSES) + 1)]
    for r in unit_rows:
        cells = [f"{f(v[2])} ({v[0] + v[1]}; {'-' if v[3] is None else v[3]})"
                 for v in (r["recall_by_class"][x] for x in CLASSES)]
        uv = r["recall_uv"]
        lines.append(f"| {r['label']} | " + " | ".join(cells) + f" | {f(uv[2])} ({uv[0] + uv[1]}) |")
    lines += ["", "## 4. Cosine of TP / FP spectra to named COSMIC signatures", "",
              "| call set | " + " | ".join(f"{s} TP / FP" for s in NAMED) + " |", "|---|" + "---|" * len(NAMED)]
    for r in unit_rows:
        cn = r.get("cos_named", {})
        lines.append(f"| {r['label']} | " + " | ".join(
            f"{f(cn.get(s, {}).get('TP'))} / {f(cn.get(s, {}).get('FP'))}" for s in NAMED) + " |")
    comp = [r for r in unit_rows if r["role"] == "comparator"]
    lines += ["", "## 5. Comparators", "",
              (f"{len(comp)} comparator call sets (rows above with role 'comparator')." if comp else
               "None yet: `../E16_bootstrap/comparators.tsv` has no active rows. Add the re-scored ClairS-TO / DeepSomatic "
               "rtg directories there (score_comparator.py `eval_bed/pass_pon`) and rerun `collect.sh`; they enter "
               "sections 2-4 as role 'comparator'."), "",
              "## 6. Checks", "",
              f"- Every unit's truth (TPb + FN) count equals its truth set's: "
              f"{sum(r['truth_consistent'] for r in unit_rows)} of {len(unit_rows)}.",
              "- Truth SNVs without a valid context (N, or REF differs from the FASTA): " +
              ", ".join(f"{r['sample']} {r['skipped']}" for r in truth_rows) + "."]
    return "\n".join(lines)


if __name__ == "__main__":
    main()

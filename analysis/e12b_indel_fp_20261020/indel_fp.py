"""Why the comparators' INDEL F1 is far above the multi-cancer models' on HCC1395T (E12b): where the false calls of the
multicancer4 / multicancer6 INDEL models are. Standard library only (RIS login node):

    python3 indel_fp.py [RUNS]     # RUNS: the copy of the tequila runs (default below)

Input: the E12b INDEL runs copied from tequila on 2026-10-20 (`<run>/vcf_<chrom>_<set>_normfix/` and
`<run>/test_<chrom>_normfix/`, no checkpoints), the E12b summary.json (comparator P / R), the RIS copy of the labels' FASTA.
Per run and test set (in the BED, after the PoN, rtg vcfeval of the pipeline):
  1. the PASS threshold, PASS calls, TP / FP / FN, the recall ceiling, the best point of the curve, and the false calls
     on the curve at the comparators' TP (TP = recall x truth records of E12b summary.json);
  2. TP, FP and FN by GRCh38 repeat class (RepeatContext of indexed_gam_pipeline_v4 truth_labels, the class R13 uses:
     HP>=7, HP4-6, STR, VNTR, none, complex) and the label reason of each false call's tensor(s) (the VCF ID lists the
     candidate ids; reason from test_<chrom>_normfix predictions);
  3. the tensors at or above the PASS threshold by label reason (repeat_indel_without_truth = R13's -1);
  4. per label group of the test chromosome's tensors (somatic, non, germline, R13 -1, other -1): count, median
     p_somatic and the share at or above the PASS threshold.
Writes tables.md next to this script.
"""
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = Path(sys.argv[1] if len(sys.argv) > 1 else
            "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs/from_tequila_20261020/runs")
FASTA = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset/REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
E12B = HERE.parent / "reviewer_response_20261011" / "experiments" / "E12b_multicancer_models" / "summary.json"
PLAT = {"Illumina": ("Illumina", "Illumina"), "PacBio": ("PacBio", "fiberseq"), "ONT": ("ONT", "ONT")}
CLASSES = ("HP>=7", "HP4-6", "STR", "VNTR", "none", "complex")


class Fasta:
    """fetch(chrom, start, end) from a FASTA with its .fai (0-based, end exclusive)."""

    def __init__(self, path):
        self.f = open(path, "rb")
        self.index = {}
        for line in open(str(path) + ".fai"):
            name, length, offset, bases, width = line.split("\t")[:5]
            self.index[name] = (int(length), int(offset), int(bases), int(width))

    def fetch(self, chrom, start, end):
        length, offset, bases, width = self.index[chrom]
        start, end = max(0, start), min(length, end)
        if end <= start:
            return ""
        a = offset + (start // bases) * width + start % bases
        b = offset + (end // bases) * width + end % bases
        self.f.seek(a)
        return self.f.read(b - a).replace(b"\n", b"").replace(b"\r", b"").decode()


# minimal_period, tract and classify: as indexed_gam_pipeline_v4/tensor_postprocessing/truth_labels.py (RepeatContext)
def minimal_period(s):
    for p in range(1, len(s) + 1):
        if len(s) % p == 0 and s[:p] * (len(s) // p) == s:
            return p
    return len(s)


def tract(seq, lo, hi, period):
    best = (lo, lo)
    i = max(0, lo - 200)
    while i < min(len(seq) - period, hi + 200):
        if seq[i] != seq[i + period]:
            i += 1
            continue
        j = i
        while j < len(seq) - period and seq[j] == seq[j + period]:
            j += 1
        a, b = i, j + period
        if a <= hi and b >= lo and b - a > best[1] - best[0] and b - a >= 2 * period:
            best = (a, b)
        i = j + 1
    return best


FLANK = 300


def repeat_class(genome, chrom, pos0, ref, alt):
    if ref and alt:
        return "complex"
    start = max(0, pos0 - FLANK)
    seq = genome.fetch(chrom, start, pos0 + len(ref) + FLANK).upper()
    event = ref or alt
    unit, lo = minimal_period(event), pos0 - start
    hi = lo + len(ref)
    if unit <= 6:
        a, b = tract(seq, lo if ref else lo - 1, hi if ref else hi + 1, unit)
    else:
        n = 0
        while seq[hi + n * unit:hi + (n + 1) * unit] == event:
            n += 1
        a, b = lo, hi + n * unit
    length = b - a
    if unit == 1:
        return "HP>=7" if length >= 7 else "HP4-6" if length >= 4 else "none"
    if unit <= 6:
        return "STR" if length >= 3 * unit else "none"
    return "VNTR" if length >= 2 * unit else "none"


def vcf_records(path):
    """(chrom, pos0, ref, alt, ids, af) of a VCF's records, anchor base removed (pos0: first deleted base / the base
    after the insertion point)."""
    out = []
    if not path.exists():
        return out
    for line in gzip.open(path, "rt"):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3].upper(), f[4].upper()
        af = None
        if len(f) > 9:
            keys, vals = f[8].split(":"), f[9].split(":")
            if "AF" in keys and keys.index("AF") < len(vals) and vals[keys.index("AF")] not in (".", ""):
                af = float(vals[keys.index("AF")].split(",")[0])
        for a in alt.split(","):
            r, x = ref, a
            k = 0  # common prefix (the anchor)
            while k < min(len(r), len(x)) and r[k] == x[k]:
                k += 1
            r, x = r[k:], x[k:]
            ids = [] if vid == "." else vid.replace(";", ",").split(",")
            out.append((chrom, pos - 1 + k, r, x, ids, af))
    return out


def roc(path):
    rows = []
    for line in gzip.open(path, "rt"):
        if line[0] == "#":
            continue
        f = line.split("\t")
        rows.append((float(f[0]), float(f[1]), float(f[2])))
    return rows


def fp_at_tp(rows, tp):
    for score, t, fp in rows:  # rtg: descending score, cumulative counts
        if t >= tp:
            return int(round(t)), int(round(fp)), score
    return None


def main():
    genome = Fasta(FASTA)
    e12b = {(r["target"], r["chrom"], r["platform"]): r for r in json.load(open(E12B)) if r["kind"] == "INDEL"}
    L = ["# Where the multi-cancer INDEL false calls are (E12b, copied from tequila 2026-10-20)\n",
         "(`indel_fp.py`; in the BED, after the PoN, rtg vcfeval of the pipeline. Repeat class: GRCh38 RepeatContext of "
         "truth_labels (the class R13 uses). DS / CT = DeepSomatic-TO / ClairS-TO, their TP = recall x truth of E12b "
         "summary.json, FP = TP (1 - P) / P.)\n"]
    head = ["\n## 1. PASS, best point and false calls at the comparators' TP\n",
            "| run | test | threshold | PASS TP / FP (FN) | ceiling | best point: p_somatic >=, TP / FP | at DS TP: TP / FP (DS FP) "
            "| at CT TP: TP / FP (CT FP) |", "|---|---|---:|---|---:|---|---|---|"]
    cls = ["\n## 2. TP, FP and FN by repeat class (PASS)\n",
           "| run | test | truth by class | TP by class | FP by class | FN by class | FP median AF (TP) |",
           "|---|---|---|---|---|---|---|"]
    rsn = ["\n## 3. Label reason of the false calls' tensors, and of all tensors >= the PASS threshold\n",
           "| run | test | false calls by reason of their tensor(s) | tensors >= threshold by label / reason |", "|---|---|---|---|"]
    grp = ["\n## 4. p_somatic of the test tensors by label group (R13 -1 = repeat_indel_without_truth)\n",
           "| run | test | somatic (1) | non (0) | germline (2) | R13 -1 | other -1 |", "|---|---|---|---|---|---|---|"]
    fmt = lambda c: ", ".join(f"{k} {c[k]:,}" for k in CLASSES if c.get(k)) or "-"  # noqa: E731
    for g in (4, 6):
        for p, (hp, cp) in PLAT.items():
            run = RUNS / f"multicancer{g}_{p}_INDEL"
            if not run.exists():
                continue
            for target, tset, chrom in (("HCC1395T", f"HCC1395T_{p}_d46", "chr1"), ("HCC1395T", f"HCC1395T_{p}_d46", "chr20"),
                                        ("HG008T", f"HG008T_{hp}_d46", "chr1"), ("COLO829T", f"COLO829T_{cp}_d46", "chr1")):
                v = run / f"vcf_{chrom}_{tset}_normfix" / "eval_bed"
                if not (v / "report.json").exists():
                    continue
                rep = json.load(open(v / "report.json"))
                thr, a, b = rep["threshold"], rep["pon"]["at_threshold"], rep["pon"]["best"]
                curve = roc(v / "vcfeval" / "pon_curve" / "weighted_roc.tsv.gz")
                n = rep["truth_records"]
                comp = []
                for c in ("deepsomatic", "clairs_to"):
                    x = (e12b.get((target, chrom, p)) or {}).get(c)
                    if not x:
                        comp.append("-")
                        continue
                    ctp = round(x["r"] * n)
                    cfp = round(ctp * (1 - x["p"]) / x["p"]) if x["p"] else 0
                    m = fp_at_tp(curve, ctp)
                    comp.append("-" if m is None else f"{m[0]} / {m[1]:,} ({cfp:,})")
                name, test = f"multicancer{g} {p}", f"{target} {chrom}"
                head.append(f"| {name} | {test} | {thr:.4f} | {a['tp_call']} / {a['fp']:,} ({a['fn']}) | {rep['pon']['ceiling']:.3f} | "
                            f"{b['p_somatic']:.3f}, {b['tp_call']} / {b['fp']:,} | {comp[0]} | {comp[1]} |")
                # repeat classes
                calls = v / "vcfeval" / "pon_calls"
                counts = {}
                afs = {}
                fp_ids = []
                for what, f in (("truth", v / "truth.INDEL.vcf.gz"), ("tp", calls / "tp.vcf.gz"), ("fp", calls / "fp.vcf.gz"),
                                ("fn", calls / "fn.vcf.gz")):
                    c = Counter()
                    af = []
                    for chrom_, pos0, r, x, ids, f_af in vcf_records(f):
                        c[repeat_class(genome, chrom_, pos0, r, x)] += 1
                        if f_af is not None:
                            af.append(f_af)
                        if what == "fp":
                            fp_ids.append(ids)
                    counts[what], afs[what] = c, sorted(af)
                med = lambda xs: f"{xs[len(xs) // 2]:.3f}" if xs else "-"  # noqa: E731
                cls.append(f"| {name} | {test} | {fmt(counts['truth'])} | {fmt(counts['tp'])} | {fmt(counts['fp'])} | "
                           f"{fmt(counts['fn'])} | {med(afs['fp'])} ({med(afs['tp'])}) |")
                # reasons
                pred = run / f"test_{chrom}_normfix" / f"{tset}.tensors.INDEL.predictions.ndjson.gz"
                reason, above, groups = {}, Counter(), {}
                if pred.exists():
                    for q in map(json.loads, gzip.open(pred, "rt")):
                        reason[q["candidate_id"]] = (q["label"], q["reason"])
                        if q["p_somatic"] >= thr:
                            above[f"{q['label']}/{q['reason']}"] += 1
                        key = {1: "somatic", 0: "non", 2: "germline"}.get(
                            q["label"], "R13" if q["reason"] == "repeat_indel_without_truth" else "other")
                        groups.setdefault(key, []).append(q["p_somatic"])
                cells = []
                for key in ("somatic", "non", "germline", "R13", "other"):
                    xs = sorted(groups.get(key, []))
                    ge = sum(x >= thr for x in xs)
                    cells.append("-" if not xs else f"{len(xs):,}, median {xs[len(xs) // 2]:.4f}, >= threshold {ge:,} "
                                                    f"({ge / len(xs):.1%})")
                grp.append(f"| {name} | {test} | " + " | ".join(cells) + " |")
                fpr = Counter()
                for ids in fp_ids:
                    hit = [reason.get(i) or reason.get(i.split("@")[0]) for i in ids]
                    got = sorted({f"{h[0]}/{h[1]}" for h in hit if h}) or ["no tensor found"]
                    fpr[got[0] if len(got) == 1 else "+".join(got)] += 1
                top = lambda c: ", ".join(f"{k} {v:,}" for k, v in c.most_common(6))  # noqa: E731
                rsn.append(f"| {name} | {test} | {top(fpr)} | {top(above)} |")
                print(name, test, "done", flush=True)
    (HERE / "tables.md").write_text("\n".join(L + head + cls + rsn + grp) + "\n")
    print("\n".join(L + head + cls + rsn + grp))


if __name__ == "__main__":
    main()

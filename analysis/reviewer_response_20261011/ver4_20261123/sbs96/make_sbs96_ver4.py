"""Response Figure 3.2 (revised Figure 4a) of letter ver4: SBS96 spectra on chr1 for the three test tumors, with
Pansoma = the six-tumor current SNV models on HPRC v2.1 (snv_newpipe_20261109/multicancer6_<P>_SNV_newpipe, seed 0).

Per tumor (HG008T, COLO829T, HCC1395T) one figure, 5 rows x 3 platforms (layout and style of
supplementary/R3.5/fig_R3.5a_<tumor>.png, panel code of supplementary/R3.5/make_sbs96.py):
  row 1  chr1 somatic truth in the benchmark regions (rtg tp-baseline + fn of the Pansoma call set),
  row 2  the matched normal's germline SNVs on chr1 in the same regions (E18 sbs96.tsv rows "germline <tumor>"),
  rows 3-5  PASS SNV calls after the population filter (rtg tp + fp, benchmark regions) of Pansoma (six-tumor model,
            HPRC v2.1), DeepSomatic-TO and ClairS-TO (E8 score, pass_pon: the tool's PASS after the same filter).
All counts are taken directly from the rtg vcfeval VCFs; the trinucleotide context is read from the GRCh38 FASTA
(no pysam: chr1 through the .fai offsets).

Checks (the script stops when one fails):
  - truth n = TP + FN of report.json (pon.at_threshold, pon.baseline) and the truth spectrum is identical in the nine
    Pansoma call sets of a tumor's three platforms and in the 2 x 3 comparator call sets (their tp-baseline + fn);
  - Pansoma n = tp_call + fp of report.json pon.at_threshold; comparator n = tp_call + fp of their report.json
    eval_bed.pass_pon.at_threshold; no record is skipped (every tp / fp / truth record is a single-base SNV on chr1);
  - HG008T, COLO829T: the comparator, truth and germline spectra (96 counts) equal the rows of
    supplementary/R3.5/sbs96.tsv, and the comparator rows of supplementary/R3.5/cosine.tsv are reproduced exactly
    (n, cosine to truth, cosine to germline, UV share);
  - with --check-germline: the germline spectrum of each tumor is recounted from the germline VCF and the BED
    (E18's rule) and compared with the E18 row.

Writes next to this script: fig_3.2_<tumor>.png (15 x 9 in, dpi 200 = 3000 x 1800 px), sbs96.tsv, cosine.tsv,
cosine_table.md, facts.json (provenance, checks and the numbers the caption / text need).
usage: /wanglab/jshen/anaconda3/envs/sambcfenv/bin/python make_sbs96_ver4.py [--check-germline] [--no-figures]
"""
import argparse
import csv
import gzip
import json
import math
import sys
from bisect import bisect_right
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
AR = HERE.parents[1]
OLD = AR / "supplementary" / "R3.5"                      # October outputs, read only (check)
E18 = AR / "experiments" / "E18_sbs96" / "sbs96.tsv"      # germline spectra
E18_SUMMARY = AR / "experiments" / "E18_sbs96" / "summary.json"
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs/snv_newpipe_20261109")
E8 = Path("/scratch/jshen/data/reviewer_experiments_20261013/E8_comparators/score")
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
CHROM = "chr1"
COMP = str.maketrans("ACGT", "TGCA")
CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
CHANNELS = [f"{a}[{s}]{b}" for s in CLASSES for a in "ACGT" for b in "ACGT"]
COLOR = {"C>A": "#1ebff0", "C>G": "#050708", "C>T": "#e62725", "T>A": "#cbcacb", "T>C": "#a1cf64", "T>G": "#edc8c5"}
# tumor -> [(platform label, dataset, model platform)]
TUMORS = {"HG008T": [("Illumina", "HG008T_Illumina", "Illumina"), ("PacBio HiFi", "HG008T_PacBio", "PacBio"),
                     ("ONT", "HG008T_ONT", "ONT")],
          "COLO829T": [("Illumina", "COLO829T_Illumina", "Illumina"), ("Fiber-seq", "COLO829T_fiberseq", "PacBio"),
                       ("ONT", "COLO829T_ONT", "ONT")],
          "HCC1395T": [("Illumina", "HCC1395T_Illumina", "Illumina"), ("PacBio HiFi", "HCC1395T_PacBio", "PacBio"),
                       ("ONT", "HCC1395T_ONT", "ONT")]}
CALLERS = ["Pansoma", "DeepSomatic-TO", "ClairS-TO"]
E8_NAME = {"DeepSomatic-TO": "deepsomatic", "ClairS-TO": "clairs_to"}
USE = "Response Figure 3.2"


def load_chrom():
    with open(FASTA + ".fai") as f:
        for line in f:
            name, ln, off, lb, lw = line.split("\t")[:5]
            if name == CHROM:
                ln, off, lb, lw = int(ln), int(off), int(lb), int(lw)
                with open(FASTA, "rb") as g:
                    g.seek(off)
                    raw = g.read(ln + ln // lb * (lw - lb) + lw)
                seq = raw.replace(b"\n", b"")[:ln].decode().upper()
                assert len(seq) == ln
                return seq
    raise SystemExit(f"{CHROM} not in {FASTA}.fai")


def sbs96_key(seq, pos, ref, alt):
    """COSMIC SBS96 channel of a single-base substitution at 1-based pos (pyrimidine strand); None if not usable."""
    ctx = seq[pos - 2:pos + 1]
    if len(ctx) != 3 or any(b not in "ACGT" for b in ctx) or ctx[1] != ref or alt not in "ACGT" or alt == ref:
        return None
    if ref in "GA":
        ctx, ref, alt = ctx.translate(COMP)[::-1], ref.translate(COMP), alt.translate(COMP)
    return f"{ctx[0]}[{ref}>{alt}]{ctx[2]}"


def vcf_counts(seq, vcfs):
    """SBS96 counts of the records of rtg output VCFs; returns (Counter, records, skipped)."""
    c, records, skipped = Counter(), 0, 0
    for v in vcfs:
        with gzip.open(v, "rt") as f:
            for line in f:
                if line[0] == "#":
                    continue
                x = line.split("\t", 6)
                assert x[0] == CHROM, (v, x[0])
                records += 1
                ref = x[3].upper()
                for alt in x[4].upper().split(","):
                    k = sbs96_key(seq, int(x[1]), ref, alt) if len(ref) == 1 and len(alt) == 1 else None
                    if k:
                        c[k] += 1
                    else:
                        skipped += 1
    return c, records, skipped


def cosine(a, b):
    dot = sum(a.get(k, 0) * b.get(k, 0) for k in CHANNELS)
    na, nb = math.sqrt(sum(a.get(k, 0) ** 2 for k in CHANNELS)), math.sqrt(sum(b.get(k, 0) ** 2 for k in CHANNELS))
    return dot / (na * nb) if na and nb else float("nan")


def uv(c):
    n = sum(c.values())
    return sum(v for k, v in c.items() if k[2:5] == "C>T" and k[0] in "CT") / n if n else float("nan")


def class_share(c, cls):
    n = sum(c.values())
    return sum(v for k, v in c.items() if k[2:5] == cls) / n if n else float("nan")


def germline_spectra():
    out = {}
    with open(E18) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["spectrum"].startswith("germline ") and row["spectrum"].split()[1] in TUMORS:
                out[row["spectrum"].split()[1]] = Counter({k: int(row[k]) for k in CHANNELS if int(row[k])})
    assert set(out) == set(TUMORS), sorted(out)
    return out


def germline_vcfs():
    out = {}
    for t in json.loads(E18_SUMMARY.read_text())["truth"]:
        if t.get("sample") in TUMORS and t.get("germline_vcf"):
            out[t["sample"]] = t["germline_vcf"]
    return out


def recount_germline(seq, vcf, bed):
    """E18's germline rule, independently: chr1 SNVs in the BED, FILTER PASS / . / GAP1 / GAP2, single-base ALT in GT."""
    starts, ends = [], []
    with open(bed) as f:
        for line in f:
            x = line.split("\t")
            if x[0] == CHROM:
                starts.append(int(x[1]))
                ends.append(int(x[2]))
    order = sorted(range(len(starts)), key=starts.__getitem__)
    starts, ends = [starts[i] for i in order], [ends[i] for i in order]
    c, seen = Counter(), False
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            if not line.startswith(CHROM + "\t"):
                if seen:
                    break
                continue
            seen = True
            x = line.rstrip("\n").split("\t")
            if not (x[6] in ("PASS", ".") or set(x[6].split(";")) <= {"GAP1", "GAP2"}):
                continue
            gt = x[9].split(":")[0].replace("|", "/").split("/") if len(x) > 9 else ["1"]
            pos, ref = int(x[1]), x[3].upper()
            for i, alt in enumerate(x[4].split(","), 1):
                if str(i) in gt and len(ref) == 1 and len(alt) == 1:
                    j = bisect_right(starts, pos - 1) - 1
                    if j >= 0 and pos - 1 < ends[j]:
                        k = sbs96_key(seq, pos, ref, alt.upper())
                        if k:
                            c[k] += 1
    return c


def tool_version(vcf):
    try:
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if not line.startswith("##"):
                    break
                low = line.lower()
                if low.startswith("##clairs_to_version=") or low.startswith("##deepvariant_version="):
                    return line.strip().split("=", 1)[1]
    except OSError:
        pass
    return None


def panel(ax, c, title, pad=4):
    n = sum(c.values())
    ax.bar(range(96), [c.get(k, 0) / n * 100 if n else 0 for k in CHANNELS], color=[COLOR[k[2:5]] for k in CHANNELS],
           width=0.8)
    ax.set_title(title, fontsize=8, loc="left", pad=pad)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for i in range(1, 6):
        ax.axvline(16 * i - 0.5, color="#bbbbbb", lw=0.5)
    ax.tick_params(axis="y", labelsize=6.5)
    ax.set_xticks([])


def class_header(ax):
    for i, c in enumerate(CLASSES):
        ax.text(16 * i + 7.5, 1.0, c, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7,
                fontweight="bold", color=COLOR[c] if c not in ("T>A", "T>G") else "#555555")


def read_old():
    """October supplementary/R3.5 outputs: (tumor, platform, spectrum) -> counts, and -> cosine row."""
    counts, cos = {}, {}
    with open(OLD / "sbs96.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["use"] == "Figure 4a":
                counts[(r["tumor"], r["platform"], r["spectrum"])] = Counter({k: int(r[k]) for k in CHANNELS if int(r[k])})
    with open(OLD / "cosine.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cos[(r["use"], r["tumor"], r["platform"], r["spectrum"])] = r
    return counts, cos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-germline", action="store_true", help="recount the germline spectra from the germline VCFs")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    seq = load_chrom()
    germ = germline_spectra()
    old_counts, old_cos = read_old()
    data, prov, checks = {}, {}, []      # data[(tumor, plat)] = {"truth": Counter, caller: Counter}
    n_truth_checks = n_n_checks = 0
    for tumor, plats in TUMORS.items():
        truth0 = None
        for plat, ds, mp in plats:
            vdir = RUNS / f"multicancer6_{mp}_SNV_newpipe" / f"vcf_chr1_{ds}_d46"
            rep = json.loads((vdir / "eval_bed" / "report.json").read_text())
            at = rep["pon"]["at_threshold"]
            w = vdir / "eval_bed" / "vcfeval" / "pon_calls"
            truth, trec, tskip = vcf_counts(seq, [w / "tp-baseline.vcf.gz", w / "fn.vcf.gz"])
            nt = sum(truth.values())
            assert tskip == 0 and trec == nt, (ds, "truth records skipped", tskip)
            assert nt == at["tp_baseline"] + at["fn"] == rep["pon"]["baseline"] == rep["raw"]["baseline"], (ds, nt, at)
            if truth0 is None:
                truth0 = truth
            assert truth == truth0, (ds, "truth spectrum differs between platforms")
            n_truth_checks += 1
            d = data[(tumor, plat)] = {"truth": truth}
            calls, crec, cskip = vcf_counts(seq, [w / "tp.vcf.gz", w / "fp.vcf.gz"])
            assert cskip == 0 and crec == sum(calls.values()) == at["tp_call"] + at["fp"], (ds, crec, cskip, at)
            n_n_checks += 1
            d["Pansoma"] = calls
            sg = rep["steps"]["graph_vcf"]
            assert Path(sg["checkpoint"]).parent == vdir.parent, sg["checkpoint"]
            prov[ds] = {"tumor": tumor, "platform": plat,
                        "pansoma": {"vcf_dir": str(vdir), "model": vdir.parent.name, "checkpoint": sg["checkpoint"],
                                    "call_rule": sg["call"], "threshold": rep["threshold"], "calls": rep["calls"],
                                    "tp_call": at["tp_call"], "fp": at["fp"], "tp_baseline": at["tp_baseline"],
                                    "fn": at["fn"], "precision": at["precision"], "recall": at["recall"],
                                    "f1": at["f1"], "pon_files": rep["pon_files"]},
                        "truth": rep["truth"], "bed": rep["bed"], "truth_n": nt}
            for lab in CALLERS[1:]:
                sd = E8 / E8_NAME[lab] / ds / "SNV"
                crep = json.loads((sd / "report.json").read_text())
                cat = crep["eval_bed"]["pass_pon"]["at_threshold"]
                cw = sd / "eval_bed" / "pass_pon"
                ctruth, _, cts = vcf_counts(seq, [cw / "tp-baseline.vcf.gz", cw / "fn.vcf.gz"])
                assert cts == 0 and ctruth == truth0, (lab, ds, "comparator truth differs from Pansoma's")
                assert sum(ctruth.values()) == cat["tp_baseline"] + cat["fn"], (lab, ds)
                assert crep["truth"] == rep["truth"] and crep["bed"] == rep["bed"], (lab, ds, "truth / BED file differs")
                assert crep["pon_files"] == rep["pon_files"], (lab, ds, "population resources differ")
                n_truth_checks += 1
                cc, ccrec, ccskip = vcf_counts(seq, [cw / "tp.vcf.gz", cw / "fp.vcf.gz"])
                assert ccskip == 0 and ccrec == sum(cc.values()) == cat["tp_call"] + cat["fp"], (lab, ds, ccrec, ccskip, cat)
                n_n_checks += 1
                d[lab] = cc
                prov[ds][lab] = {"score_dir": str(sd), "input_vcfs": crep["input_vcfs"],
                                 "version": tool_version(crep["input_vcfs"][0]), "tp_call": cat["tp_call"],
                                 "fp": cat["fp"], "precision": cat["precision"], "recall": cat["recall"], "f1": cat["f1"]}
    checks.append(f"truth n = TP + FN of report.json and truth spectrum identical across platforms and callers: "
                  f"{n_truth_checks} call sets (9 Pansoma + 18 comparator), all equal")
    checks.append(f"n of every call set = tp_call + fp of its report.json (Pansoma pon.at_threshold; comparators "
                  f"eval_bed.pass_pon.at_threshold), no record skipped: {n_n_checks} call sets, all equal")

    # reproduce the October comparator / truth / germline numbers for HG008T and COLO829T
    n_old = 0
    for tumor in ("HG008T", "COLO829T"):
        for plat, ds, _ in TUMORS[tumor]:
            d = data[(tumor, plat)]
            assert d["truth"] == old_counts[(tumor, plat, "chr1 somatic truth")], (ds, "truth vs October sbs96.tsv")
            assert germ[tumor] == old_counts[(tumor, plat, "matched-normal germline")], (ds, "germline vs October")
            for lab in CALLERS[1:]:
                assert d[lab] == old_counts[(tumor, plat, lab)], (ds, lab, "96 counts differ from October sbs96.tsv")
                o = old_cos[("Figure 4a", tumor, plat, lab)]
                new = (str(sum(d[lab].values())), f"{cosine(d[lab], d['truth']):.3f}",
                       f"{cosine(d[lab], germ[tumor]):.3f}", f"{uv(d[lab]):.3f}")
                assert new == (o["n"], o["cosine_to_truth"], o["cosine_to_normal_germline"], o["uv_type_ct_share"]), (ds, lab, new, o)
                n_old += 1
        o = old_cos[("reference", tumor, "-", "matched-normal germline")]
        t = data[(tumor, TUMORS[tumor][0][0])]["truth"]
        assert (str(sum(germ[tumor].values())), f"{cosine(germ[tumor], t):.3f}") == (o["n"], o["cosine_to_truth"]), (tumor, o)
    checks.append(f"HG008T / COLO829T: comparator spectra (96 counts), n, cosine to truth, cosine to germline and UV "
                  f"share equal supplementary/R3.5/sbs96.tsv and cosine.tsv in {n_old} of {n_old} comparator rows; "
                  f"truth and germline spectra and the truth-germline cosine also equal")

    if args.check_germline:
        gv = germline_vcfs()
        for tumor, plats in TUMORS.items():
            bed = prov[plats[0][1]]["bed"]
            c = recount_germline(seq, gv[tumor], bed)
            assert c == germ[tumor], (tumor, sum(c.values()), sum(germ[tumor].values()))
            checks.append(f"germline {tumor}: recounted from {gv[tumor]} in {bed} (chr1, FILTER PASS / . / GAP1 / GAP2, "
                          f"single-base ALT in the genotype): {sum(c.values()):,} SNVs, 96 counts equal the E18 row")
            for _, ds, _ in plats:
                prov[ds]["germline_vcf"] = gv[tumor]

    # tables
    tsv, cos = [], []
    for tumor, plats in TUMORS.items():
        for plat, ds, _ in plats:
            d, g = data[(tumor, plat)], germ[tumor]
            rows = [("chr1 somatic truth", d["truth"]), ("matched-normal germline", g)] + [(lab, d[lab]) for lab in CALLERS]
            for lab, c in rows:
                tsv.append([USE, tumor, plat, lab, sum(c.values())] + [c.get(k, 0) for k in CHANNELS])
            for lab in CALLERS:
                c = d[lab]
                cos.append([USE, tumor, plat, lab, sum(c.values()), f"{cosine(c, d['truth']):.3f}",
                            f"{cosine(c, g):.3f}", f"{uv(c):.3f}"])
    for tumor, plats in TUMORS.items():
        t, g = data[(tumor, plats[0][0])]["truth"], germ[tumor]
        cos.append(["reference", tumor, "-", "matched-normal germline", sum(g.values()), f"{cosine(g, t):.3f}", "1.000",
                    f"{uv(g):.3f}"])
        cos.append(["reference", tumor, "-", "chr1 somatic truth", sum(t.values()), "1.000", f"{cosine(t, g):.3f}",
                    f"{uv(t):.3f}"])
    with open(HERE / "sbs96.tsv", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["use", "tumor", "platform", "spectrum", "n"] + CHANNELS)
        w.writerows(tsv)
    with open(HERE / "cosine.tsv", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["use", "tumor", "platform", "spectrum", "n", "cosine_to_truth", "cosine_to_normal_germline",
                    "uv_type_ct_share"])
        w.writerows(cos)

    # Markdown table and the numbers a caption / sentence needs (all formatted here, from the counts)
    crow = {(r[1], r[2], r[3]): r for r in cos}
    md = ["Cosine similarity of SBS96 spectra to the chr1 somatic truth (chr1, benchmark regions; PASS SNV calls after "
          "the population filter, TP + FP; Pansoma = six-tumor model on HPRC v2.1, seed 0).", "",
          "| Tumor | Platform | chr1 truth SNVs | Pansoma, n | Pansoma, cosine to truth | DeepSomatic-TO, n | "
          "DeepSomatic-TO, cosine to truth | ClairS-TO, n | ClairS-TO, cosine to truth | Truth vs matched-normal "
          "germline, cosine |", "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    md2 = ["Cosine similarity of the same call sets to the matched normal's germline spectrum (chr1, same regions) and "
           "the share of T>C substitutions.", "",
           "| Tumor | Platform | Germline SNVs | Pansoma, cosine to germline | DeepSomatic-TO, cosine to germline | "
           "ClairS-TO, cosine to germline | T>C share: truth | germline | Pansoma | DeepSomatic-TO | ClairS-TO |",
           "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    facts = {"per_dataset": {}, "per_tumor": {}}
    closest, farthest = [], []
    for tumor, plats in TUMORS.items():
        g = germ[tumor]
        t0 = data[(tumor, plats[0][0])]["truth"]
        tg = crow[(tumor, "-", "matched-normal germline")][5]
        ft = facts["per_tumor"][tumor] = {"truth_n": sum(t0.values()), "germline_n": sum(g.values()),
                                          "cosine_truth_germline": tg, "truth_uv_share": f"{uv(t0):.3f}",
                                          "germline_uv_share": f"{uv(g):.3f}",
                                          "truth_TtoC_share": f"{class_share(t0, 'T>C'):.3f}",
                                          "germline_TtoC_share": f"{class_share(g, 'T>C'):.3f}"}
        for lab in CALLERS:
            vt = [crow[(tumor, p, lab)][5] for p, _, _ in plats]
            vg = [crow[(tumor, p, lab)][6] for p, _, _ in plats]
            ft[lab] = {"cosine_to_truth": " / ".join(vt), "cosine_to_truth_range": f"{min(vt)}-{max(vt)}",
                       "cosine_to_germline": " / ".join(vg), "cosine_to_germline_range": f"{min(vg)}-{max(vg)}",
                       "n": " / ".join(f"{crow[(tumor, p, lab)][4]:,}" for p, _, _ in plats)}
        for plat, ds, _ in plats:
            d = data[(tumor, plat)]
            r = {lab: crow[(tumor, plat, lab)] for lab in CALLERS}
            md.append(f"| {tumor} | {plat} | {sum(d['truth'].values()):,} | " +
                      " | ".join(f"{r[lab][4]:,} | {r[lab][5]}" for lab in CALLERS) + f" | {tg} |")
            md2.append(f"| {tumor} | {plat} | {sum(g.values()):,} | " + " | ".join(r[lab][6] for lab in CALLERS) +
                       f" | {class_share(d['truth'], 'T>C'):.3f} | {class_share(g, 'T>C'):.3f} | " +
                       " | ".join(f"{class_share(d[lab], 'T>C'):.3f}" for lab in CALLERS) + " |")
            exact_t = {lab: cosine(d[lab], d["truth"]) for lab in CALLERS}
            exact_g = {lab: cosine(d[lab], g) for lab in CALLERS}
            best_t = max(exact_t, key=exact_t.get)
            low_g = min(exact_g, key=exact_g.get)
            closest.append((ds, best_t))
            farthest.append((ds, low_g))
            facts["per_dataset"][ds] = {
                "truth_n": sum(d["truth"].values()),
                "n": {lab: sum(d[lab].values()) for lab in CALLERS},
                "cosine_to_truth": {lab: r[lab][5] for lab in CALLERS},
                "cosine_to_germline": {lab: r[lab][6] for lab in CALLERS},
                "uv_share": {lab: r[lab][7] for lab in CALLERS},
                "TtoC_share": {lab: f"{class_share(d[lab], 'T>C'):.3f}" for lab in CALLERS},
                "closest_to_truth_unrounded": best_t,
                "pansoma_minus_best_comparator_cosine_to_truth": f"{exact_t['Pansoma'] - max(exact_t[x] for x in CALLERS[1:]):+.4f}",
                "pansoma_closest_or_tied_at_3_decimals": r["Pansoma"][5] >= max(r[x][5] for x in CALLERS[1:]),
                "farthest_from_germline_unrounded": low_g,
                "pansoma_minus_lowest_comparator_cosine_to_germline": f"{exact_g['Pansoma'] - min(exact_g[x] for x in CALLERS[1:]):+.4f}"}
    facts["pansoma_closest_to_truth_unrounded"] = f"{sum(1 for _, b in closest if b == 'Pansoma')} of {len(closest)}"
    facts["pansoma_closest_or_tied_at_3_decimals"] = \
        f"{sum(1 for v in facts['per_dataset'].values() if v['pansoma_closest_or_tied_at_3_decimals'])} of {len(closest)}"
    facts["not_pansoma_closest_to_truth"] = [f"{ds}: {b}" for ds, b in closest if b != "Pansoma"]
    facts["pansoma_farthest_from_germline_unrounded"] = f"{sum(1 for _, b in farthest if b == 'Pansoma')} of {len(farthest)}"
    facts["not_pansoma_farthest_from_germline"] = [f"{ds}: {b}" for ds, b in farthest if b != "Pansoma"]
    facts["checks"] = checks
    facts["provenance"] = prov
    facts["germline_source"] = str(E18)
    (HERE / "cosine_table.md").write_text("\n".join(md + [""] + md2) + "\n")
    (HERE / "facts.json").write_text(json.dumps(facts, indent=1) + "\n")
    print("\n".join(md + [""] + md2))
    print()
    for c in checks:
        print("CHECK OK:", c)

    if args.no_figures:
        return 0
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for tumor, plats in TUMORS.items():
        g = germ[tumor]
        fig, axes = plt.subplots(5, 3, figsize=(15, 9.0), sharey="row")
        for j, (plat, ds, _) in enumerate(plats):
            d = data[(tumor, plat)]
            truth = d["truth"]
            rows = [("chr1 somatic truth", truth), ("matched-normal germline", g)] + [(lab, d[lab]) for lab in CALLERS]
            for i, (lab, c) in enumerate(rows):
                n = sum(c.values())
                if i == 0:
                    t = f"{tumor} {plat}: {lab} (n = {n:,})"
                elif i == 1:
                    t = f"{lab} (n = {n:,}; cosine to truth {cosine(c, truth):.3f})"
                else:
                    t = f"{lab} (n = {n:,}; cosine to truth {cosine(c, truth):.3f}, to germline {cosine(c, g):.3f})"
                panel(axes[i, j], c, t, pad=14 if i == 0 else 4)
            class_header(axes[0, j])
        for i in range(5):
            axes[i, 0].set_ylabel("% of SNVs", fontsize=7)
        for j in range(3):
            axes[4, j].set_xticks(range(96))
            # 3.2 pt (October figure: 3.8 pt, neighbouring labels touched at 96 labels per 4.4 in panel)
            axes[4, j].set_xticklabels([k[0] + k[2] + k[6] for k in CHANNELS], rotation=90, fontsize=3.2,
                                       family="monospace")
            axes[4, j].tick_params(axis="x", length=2.5, width=0.5, pad=1.5)
        fig.suptitle(f"{tumor}, chr1, PASS SNVs after population filtering, benchmark regions "
                     f"(Pansoma: six-tumor model, HPRC v2.1)", fontsize=9, x=0.01, ha="left")
        fig.tight_layout()
        fig.align_ylabels(axes[:, 0])
        out = HERE / f"fig_3.2_{tumor}.png"
        fig.savefig(out, dpi=200)
        plt.close(fig)
        w, h = (int(v) for v in fig.get_size_inches() * 200)
        print("wrote", out, f"{w} x {h} px")
    return 0


if __name__ == "__main__":
    sys.exit(main())

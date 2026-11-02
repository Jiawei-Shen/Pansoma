"""COLO829T chr1 somatic truth by expected VAF bin: truth counts, candidate stage (tensor / filtered / no candidate) per
tensor set, and final PASS recall (rtg pon_calls tp-baseline / fn) per local call set; plus the 200-row HapMap SNV
model on the P100 low-AF rebuild (E15b, tensor level, from its predictions file).

Expected VAF of a truth record: truth VCF VAF (VAF_Ill for Illumina, VAF_PB for fiberseq / ONT) x the measured purity
of the set (purity.md of the 7fc7396 purity benchmark: median somatic SNV tensor AF relative to P100; Illumina /
fiberseq 0.66 / 0.40 / 0.20, ONT 0.70 / 0.44 / 0.22 for P75 / P50 / P25). Truth = chr1, PASS, in the BED (the
somatic.recall.tsv 'passed' / 'in_bed' columns, as the pipeline's truth metrics).
Bins: [0, 0.05), [0.05, 0.10), [0.10, 0.20), [0.20, 1].
usage: python dilution_vaf.py > dilution_vaf.md
"""
import gzip
import json
from pathlib import Path

D = Path("/scratch/jshen/data/pansoma_v2_tensors")
R = Path("/scratch/jshen/data/pansoma_net_v2_runs")
TRUTH = D / "COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz"
FACTOR = {"Illumina": {"P100": 1.0, "P75": 0.66, "P50": 0.40, "P25": 0.20},
          "fiberseq": {"P100": 1.0, "P75": 0.66, "P50": 0.40, "P25": 0.20},
          "ONT": {"P100": 1.0, "P75": 0.70, "P50": 0.44, "P25": 0.22}}
BINS = ((0.0, 0.05), (0.05, 0.10), (0.10, 0.20), (0.20, 9.0))
BIN_NAMES = ("<0.05", "0.05-0.10", "0.10-0.20", ">=0.20")
SNV_KINDS = {"SNP", "SNV"}


def vaf_field(platform):
    return "VAF_Ill" if platform == "Illumina" else "VAF_PB"


def bin_of(v):
    for i, (lo, hi) in enumerate(BINS):
        if lo <= v < hi:
            return i
    return None


# truth records in VCF order (truth_id = record index, checked against somatic.recall.tsv below)
records = []
with gzip.open(TRUTH, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
        records.append(dict(chrom=c[0], pos=int(c[1]), ref=c[3], alt=c[4],
                            VAF_Ill=float(info["VAF_Ill"]) if info.get("VAF_Ill", ".") != "." else None,
                            VAF_PB=float(info["VAF_PB"]) if info.get("VAF_PB", ".") != "." else None))
by_key = {(r["chrom"], r["pos"], r["ref"], r["alt"]): r for r in records}


def recall_table(set_name):
    out = {}
    with open(D / set_name / "tensors/somatic.recall.tsv") as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            tid = int(r["truth_id"])
            assert records[tid]["chrom"] == r["chrom"] and records[tid]["pos"] == int(r["vcf_pos"]), (set_name, tid)
            if r["chrom"] == "chr1" and r["passed"] == "True" and r["in_bed"] == "True":
                out[tid] = (("SNV" if r["kind"] in SNV_KINDS else "INDEL"), r["status"])
    return out


def stage(status):
    if status.startswith("tensor"):
        return "tensor"
    return status  # filtered / no_candidate


def candidate_rows(platform, purity, set_name, label):
    tab = recall_table(set_name)
    fac = FACTOR[platform][purity]
    rows = []
    for kind in ("SNV", "INDEL"):
        cnt = [[0, 0, 0, 0] for _ in BINS]  # truth, tensor, filtered, no_candidate
        unknown = 0
        for tid, (k, st) in tab.items():
            if k != kind:
                continue
            v = records[tid][vaf_field(platform)]
            if v is None:
                unknown += 1
                continue
            b = bin_of(v * fac)
            cnt[b][0] += 1
            cnt[b][{"tensor": 1, "filtered": 2, "no_candidate": 3}[stage(st)]] += 1
        tot = [sum(c[i] for c in cnt) for i in range(4)]
        rows.append((label, platform, purity, kind, cnt, tot, unknown))
    return rows


def rtg_recall(eval_dir, platform, purity):
    """per bin TP / truth from rtg pon_calls (tp-baseline, fn) and the overall P / R / F1 from report.json."""
    fac = FACTOR[platform][purity]
    vc = Path(eval_dir) / "vcfeval/pon_calls"
    tp, n = [0] * len(BINS), [0] * len(BINS)
    for fname, hit in (("tp-baseline.vcf.gz", 1), ("fn.vcf.gz", 0)):
        with gzip.open(vc / fname, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split("\t", 5)
                r = by_key.get((c[0], int(c[1]), c[3], c[4]))
                if r is None or r[vaf_field(platform)] is None:
                    continue
                b = bin_of(r[vaf_field(platform)] * fac)
                n[b] += 1
                tp[b] += hit
    rep = json.load(open(Path(eval_dir) / "report.json"))
    a = rep["pon"]["at_threshold"]
    return tp, n, (a["precision"], a["recall"], a["f1"], rep.get("unplaced_pass"))


def fmt(tp, n):
    return f"{tp / n:.3f} ({tp}/{n})" if n else "-"


print("# COLO829T chr1 by expected truth VAF (computed 2026-10-30, dilution_vaf.py)\n")
print("Expected VAF = truth VCF VAF (VAF_Ill Illumina, VAF_PB fiberseq / ONT) x measured purity factor "
      "(purity.md: 0.66 / 0.40 / 0.20; ONT 0.70 / 0.44 / 0.22). chr1, PASS truth in the BED.\n")
print("## 1. Candidate stage per tensor set (cells: truth n; tensor / filtered / no candidate)\n")
print("| set | platform | purity | kind | " + " | ".join(BIN_NAMES) + " | all | VAF unknown |")
print("|---|---|---|---|" + "---|" * (len(BINS) + 2))
sets = []
for p, base in (("Illumina", "COLO829T_Illumina"), ("fiberseq", "COLO829T_fiberseq"), ("ONT", "COLO829T_ONT")):
    for q in ("P100", "P75", "P50", "P25"):
        sets.append((p, q, base + ("" if q == "P100" else f"_{q}"), "default build (d9, AF 0.08)"))
for q in ("P100", "P75", "P50", "P25"):
    sets.append(("Illumina", q, "COLO829T_Illumina" + ("" if q == "P100" else f"_{q}") + "_lowaf_chr1",
                 "low-AF rebuild (d9, AF 0.015, node-alt 0.01)"))
for p, q, s, label in sets:
    for lab, plat, pur, kind, cnt, tot, unk in candidate_rows(p, q, s, label):
        cells = [f"{c[0]}; {c[1]} / {c[2]} / {c[3]}" for c in cnt]
        print(f"| {s} ({lab}) | {plat} | {pur} | {kind} | " + " | ".join(cells) +
              f" | {tot[0]}; {tot[1]} / {tot[2]} / {tot[3]} | {unk} |")

print("\n## 2. Final PASS recall after the PoN per expected-VAF bin (rtg pon_calls; cells: R (TP/truth))\n")
print("| call set | model | platform | purity | " + " | ".join(BIN_NAMES) + " | P / R / F1 (report.json pon) | unplaced PASS |")
print("|---|---|---|---|" + "---|" * (len(BINS) + 2))
calls = []
for q in ("P100", "P75", "P50", "P25"):
    s = "COLO829T_fiberseq" + ("" if q == "P100" else f"_{q}")
    calls.append((R / f"reviewer_cpu_20261013/COLO829_fiberseq_SNV_joint_7fc7396_gpu4/vcf_chr1_{s}/eval_bed",
                  "joint 4-purity SNV, d9 (COLO829T chr2-22 of all purities)", "fiberseq", q))
for m, lab in (("COLO829T_ONT_SNV_joint_control_tq", "joint 4-purity SNV, d9, E14b control"),
               ("COLO829T_ONT_SNV_joint_readnorm_tq", "joint 4-purity SNV, d9, E14b per-read BQ")):
    for q in ("P100", "P75", "P50", "P25"):
        calls.append((R / f"{m}/vcf_chr1_{q}/eval_bed", lab, "ONT", q))
for p, ds, d46 in (("Illumina", "COLO829T_Illumina_d46", "COLO829_Illumina_SNV_d46_7fc7396"),
                   ("fiberseq", "COLO829T_fiberseq_d46", "COLO829_fiberseq_SNV_d46_7fc7396"),
                   ("ONT", "COLO829T_ONT_d46", "COLO829_ONT_SNV_d46_7fc7396")):
    pl = {"Illumina": "Illumina", "fiberseq": "PacBio", "ONT": "ONT"}[p]
    calls.append((R / f"reviewer_cpu_20261013/{d46}/vcf_chr1/eval_bed", "single-tumor COLO829T SNV, d46 (in-sample)", p, "P100"))
    calls.append((R / f"multicancer6_{pl}_SNV/vcf_chr1_{ds}_normfix/eval_bed", "multicancer6 SNV, d46 (COLO829T chr2-22 in training)", p, "P100"))
    calls.append((R / f"multicancer4_{pl}_SNV/vcf_chr1_{ds}_normfix/eval_bed", "multicancer4 SNV, d46 (COLO829T unseen)", p, "P100"))
for d, lab, p, q in calls:
    tp, n, (P, Rr, F, u) = rtg_recall(d, p, q)
    print(f"| `{d}` | {lab} | {p} | {q} | " + " | ".join(fmt(a, b) for a, b in zip(tp, n)) +
          f" | {P:.3f} / {Rr:.3f} / {F:.3f} | {u} |")

# E15b: HapMap 200-row SNV model on the P100 low-AF rebuild, tensor level (truth called if any of its tensors
# reaches the threshold); thresholds from the predictions' metrics.json
pred = R / "HapMap_Illumina_SNV_full/test_chr1_colo829t_illumina/COLO829T_Illumina_lowaf_chr1.tensors.SNV.predictions.ndjson.gz"
met = json.load(open(str(pred).replace("predictions.ndjson.gz", "metrics.json")))
thr = {"pipeline (val recall 0.9)": met["rules"]["pipeline"]["threshold"], "threshold_t": met["rules"]["threshold_t"]["threshold"]}
best = {}
with gzip.open(pred, "rt") as f:
    for line in f:
        if '"label": 1' not in line:
            continue
        r = json.loads(line)
        if r["label"] != 1 or r["excluded"]:
            continue
        for t in r["truth_ids"]:
            best[t] = max(best.get(t, 0.0), r["p_somatic"])
tab = recall_table("COLO829T_Illumina_lowaf_chr1")
print("\n## 3. E15b: HapMap 200-row SNV model (HapMap_Illumina_SNV_full) on COLO829T Illumina P100, low-AF rebuild, "
      "tensor level, no PoN (cells: R (called/truth))\n")
print("| threshold | value | " + " | ".join(BIN_NAMES) + " | all | truth-level P / R / F1 (metrics.json) |")
print("|---|---|" + "---|" * (len(BINS) + 2))
for name, t in thr.items():
    tp, n = [0] * len(BINS), [0] * len(BINS)
    for tid, (k, st) in tab.items():
        if k != "SNV":
            continue
        v = records[tid]["VAF_Ill"]
        if v is None:
            continue
        b = bin_of(v)
        n[b] += 1
        tp[b] += int(best.get(tid, 0.0) >= t)
    key = "pipeline" if name.startswith("pipeline") else "threshold_t"
    tr = met["rules"][key]["truth"]
    print(f"| {name} | {t:.4f} | " + " | ".join(fmt(a, b) for a, b in zip(tp, n)) +
          f" | {fmt(sum(tp), sum(n))} | {tr['precision']:.3f} / {tr['recall']:.3f} / {tr['f1']:.3f} |")

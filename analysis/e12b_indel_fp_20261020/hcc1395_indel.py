"""HCC1395T INDEL (E12b): Pansoma's multi-cancer models against DeepSomatic-TO and ClairS-TO, call by call. Standard
library only (RIS login node):

    python3 hcc1395_indel.py

Input (copied from tequila 2026-10-20 into $PANSOMA_RUNS/from_tequila_20261020/): runs/multicancer{4,6}_<P>_INDEL
(vcf_<chrom>_HCC1395T_<P>_d46_normfix, test_<chrom>_normfix), E8_score/{deepsomatic,clairs_to}/HCC1395T_<P>[_chr20]/INDEL
(score_comparator.py output), tensors/<set>/tensors/INDEL/labels.manifest.json and HCC1395T chr1 / chr20 labels.
Every caller is scored with the same truth (SEQC2 v1.2.1), BED, PoN and rtg vcfeval (pansoma_net_v2/vcfeval.py);
here: Pansoma PASS after the PoN (vcfeval pon_calls), comparators PASS after the PoN (pass_pon).
  1. training labels of the four lines (labels.manifest.json totals; R13 = repeat_indel_without_truth);
  2. per test (chr1, chr20) and platform, per caller: TP / FP / FN, and TP (truth side), FP and FN by GRCh38 repeat
     class (indel_fp.repeat_class = truth_labels RepeatContext);
  3. Pansoma's false calls in HP>=7 / STR: are they records of the comparators' raw output, and with which FILTER;
  4. check: indel_fp.repeat_class against the labels' own repeat_class for the chr1 tensors R13 judged.
Writes hcc1395.md next to this script.
"""
import gzip
import importlib.util
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("indel_fp", HERE / "indel_fp.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
BASE = Path("/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs/from_tequila_20261020")
RUNS, SCORE, TENSORS = BASE / "runs", BASE / "E8_score", BASE / "tensors"
PLATFORMS = ("Illumina", "PacBio", "ONT")
LINES = ("HCC1395T", "HCC1954T", "H1437T", "H2009T")
CLASSES = F.CLASSES
REPEAT = ("HP>=7", "STR")


def key(rec):
    return rec[0], rec[1], rec[2], rec[3]


def classes(genome, recs):
    return Counter(F.repeat_class(genome, r[0], r[1], r[2], r[3]) for r in recs)


def raw_filters(path):
    """(chrom, pos0, ref, alt) -> FILTER of every record of a comparator's calls.all.vcf.gz (anchor removed)."""
    out = {}
    for line in gzip.open(path, "rt"):
        if line[0] == "#":
            continue
        f = line.split("\t", 8)
        ref, alt = f[3].upper(), f[4].upper()
        k = 0
        while k < min(len(ref), len(alt)) and ref[k] == alt[k]:
            k += 1
        out[(f[0], int(f[1]) - 1 + k, ref[k:], alt[k:])] = f[6]
    return out


def main():
    genome = F.Fasta(F.FASTA)
    L = ["# HCC1395T INDEL: the multi-cancer models against DeepSomatic-TO and ClairS-TO, call by call\n",
         "(`hcc1395_indel.py`; SEQC2 v1.2.1 truth, SEQC2 high-confidence BED, the same PoN and rtg vcfeval for every "
         "caller; PASS after the PoN. Repeat class: GRCh38 RepeatContext of truth_labels. DeepSomatic-TO was trained on "
         "HCC1395 chr2-20 (chr1 held out, E8b D2 / D10): its chr20 rows are training data; ClairS-TO never saw HCC1395.)\n"]

    # 1. training labels
    L += ["\n## 1. Training labels of the multi-cancer lines (INDEL, all chromosomes; R13 = `repeat_indel_without_truth`)\n",
          "| set | somatic | non (0) | germline | R13 -1 | R13 / (non + R13) |", "|---|---:|---:|---:|---:|---:|"]
    for p in PLATFORMS:
        for line in LINES:
            m = json.load(open(TENSORS / f"{line}_{p}_d46" / "tensors" / "INDEL" / "labels.manifest.json"))
            t, r = m["totals"], m["reasons"]
            r13 = r.get("repeat_indel_without_truth", 0)
            L.append(f"| {line} {p} | {t.get('somatic', 0):,} | {t.get('non', 0):,} | {t.get('germline', 0):,} | {r13:,} | "
                     f"{r13 / max(1, t.get('non', 0) + r13):.0%} |")

    # 2. per caller
    L += ["\n## 2. TP / FP / FN by repeat class (PASS after the PoN)\n",
          "| test | platform | caller | TP / FP / FN | P / R / F1 | truth found by class (of the truth's) | FP by class | "
          "FP not HP>=7 / STR |", "|---|---|---|---|---|---|---|---:|"]
    overlap = ["\n## 3. Pansoma's false calls in HP>=7 / STR in the comparators' raw output (any FILTER)\n",
               "| test | platform | model | its HP>=7 / STR false calls | in DeepSomatic's records: FILTER | "
               "in ClairS-TO's records: FILTER |", "|---|---|---|---:|---|---|"]
    for chrom in ("chr1", "chr20"):
        for p in PLATFORMS:
            sfx = "" if chrom == "chr1" else "_chr20"
            callers = {f"multicancer{g}": RUNS / f"multicancer{g}_{p}_INDEL" / f"vcf_{chrom}_HCC1395T_{p}_d46_normfix" /
                       "eval_bed" / "vcfeval" / "pon_calls" for g in (4, 6)}
            callers.update({c: SCORE / c / f"HCC1395T_{p}{sfx}" / "INDEL" / "eval_bed" / "pass_pon"
                            for c in ("deepsomatic", "clairs_to")})
            truth_v = RUNS / f"multicancer4_{p}_INDEL" / f"vcf_{chrom}_HCC1395T_{p}_d46_normfix" / "eval_bed" / "truth.INDEL.vcf.gz"
            truth = classes(genome, F.vcf_records(truth_v))
            raws = {c: raw_filters(SCORE / c / f"HCC1395T_{p}{sfx}" / "INDEL" / "calls.all.vcf.gz")
                    for c in ("deepsomatic", "clairs_to")}
            for name, d in callers.items():
                if not (d / "fp.vcf.gz").exists():
                    L.append(f"| {chrom} | {p} | {name} | - | - | - | - | - |")
                    continue
                tpb, fp, fn = (F.vcf_records(d / f) for f in ("tp-baseline.vcf.gz", "fp.vcf.gz", "fn.vcf.gz"))
                tp = len(tpb)
                tpc = len(F.vcf_records(d / "tp.vcf.gz"))
                P = tpc / max(1, tpc + len(fp))
                R = tp / max(1, tp + len(fn))
                f1 = 2 * P * R / (P + R) if P + R else 0
                ct, cf = classes(genome, tpb), classes(genome, fp)
                found = ", ".join(f"{k} {ct.get(k, 0)}/{truth[k]}" for k in CLASSES if truth.get(k))
                fpc = ", ".join(f"{k} {cf[k]:,}" for k in CLASSES if cf.get(k)) or "-"
                other = sum(v for k, v in cf.items() if k not in REPEAT)
                label = {"deepsomatic": "DeepSomatic-TO", "clairs_to": "ClairS-TO"}.get(name, name)
                L.append(f"| {chrom} | {p} | {label} | {tp} / {len(fp):,} / {len(fn)} | {P:.3f} / {R:.3f} / {f1:.3f} | "
                         f"{found} | {fpc} | {other:,} |")
                if name.startswith("multicancer"):
                    rep = [r for r in fp if F.repeat_class(genome, r[0], r[1], r[2], r[3]) in REPEAT]
                    cells = []
                    for c in ("deepsomatic", "clairs_to"):
                        got = Counter(raws[c].get(key(r), "not a record") for r in rep)
                        cells.append(", ".join(f"{k} {v:,}" for k, v in got.most_common(5)))
                    overlap.append(f"| {chrom} | {p} | {name} | {len(rep):,} | {cells[0]} | {cells[1]} |")
            print(chrom, p, "done", flush=True)

    # 4. check of the repeat classification
    chk = ["\n## 4. Check: indel_fp.repeat_class against the labels' repeat_class (HCC1395T chr1 tensors R13 judged)\n",
           "| set | tensors with repeat_class | same class | differ |", "|---|---:|---:|---:|"]
    for p in PLATFORMS:
        same = diff = 0
        for i, line in enumerate(open(TENSORS / f"HCC1395T_{p}_d46" / "tensors" / "INDEL" / "chr1_labels.ndjson")):
            q = json.loads(line)
            rc, lin = q.get("repeat_class"), q.get("grch38")
            if rc is None or not lin:
                continue
            mine = F.repeat_class(genome, lin["chrom"], lin["pos0"], lin["ref"], lin["alt"])
            same += mine == rc
            diff += mine != rc
            if same + diff >= 20000:
                break
        chk.append(f"| HCC1395T {p} chr1 (first 20,000) | {same + diff:,} | {same:,} | {diff:,} |")
    out = L + overlap + chk
    (HERE / "hcc1395.md").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()

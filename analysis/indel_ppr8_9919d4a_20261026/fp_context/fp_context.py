"""PASS false positives of the d46 and GRCh38-only (linear) ppr8 runs by repeat context, seeds 0-2 (chr1, after the PoN,
in the benchmark BED).

    ris/submit.sh -g 0 -n 2 -m 8 pnv2_fp_context py.sh $PWD/analysis/indel_ppr8_9919d4a_20261026/fp_context/fp_context.py

FPs: <run>/vcf_chr1/eval_bed/vcfeval/pon_calls/fp.vcf.gz. Each FP of one arm is
  shared        the other arm has a PASS FP within 10 bp;
  other cand.   no FP of the other arm there, but the other arm has a candidate (a record of its post-PoN curve input,
                i.e. any prediction left by the PoN) within 10 bp;
  own only      the other arm has no candidate within 10 bp; split by whether the other arm had one before the PoN
                (its raw curve input, every prediction): "PoN removed" or "no candidate".
Classes as E17a (../../reviewer_response_20261011/experiments/E17a_indel_strata): type INS / DEL / complex
(truth_labels.split_alleles) and repeat context truth_labels.RepeatContext.classify (HP>=7, HP4-6, STR unit 2-6 >= 3
units, VNTR unit > 6 >= 2 units, none) on the GRCh38 no-alt FASTA; the label pipeline's code, imported unchanged.
For reference: the context of the chr1 truth INDELs in the BED (tp-baseline + fn of the d46 seed-0 run).
Writes tables.md and summary.json next to this script.
"""
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
from indexed_gam_pipeline_v4.tensor_postprocessing.truth_labels import RepeatContext, split_alleles  # noqa: E402

RUNS = Path(os.environ.get("PANSOMA_RUNS", "/storage2/fs1/epigenome/Active/jiawei.s/pansoma_net_v2_runs"))
FASTA = Path(os.environ.get("PANSOMA_DATA", "/storage2/fs1/epigenome/Active/jiawei.s/Dataset")) / \
    "REF/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
SETS = (("COLO829", "Illumina"), ("COLO829", "fiberseq"), ("HG008", "Illumina"), ("HG008", "PacBio"))
SEEDS = ("", "_s1", "_s2")
CONTEXTS = ("HP>=7", "HP4-6", "STR", "VNTR", "none", "complex")
GROUPS = ("shared", "other cand.", "own only")
OWN = ("own only: PoN removed", "own only: no candidate")
WINDOW = 10


def records(path):
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 6)
            out.append((c[0], int(c[1]), c[3].upper(), c[4].split(",")[0].upper()))
    return out


def near(pos, sorted_positions):
    import bisect
    i = bisect.bisect_left(sorted_positions, pos - WINDOW)
    return i < len(sorted_positions) and sorted_positions[i] <= pos + WINDOW


class Classes:
    def __init__(self):
        self.repeats = RepeatContext(FASTA)

    def __call__(self, chrom, pos, ref, alt):
        for pos0, r, a, kind, _ in split_alleles(pos, ref, [alt]):
            if kind in ("INS", "DEL"):
                return kind, self.repeats.classify(chrom, pos0, r, a)
            if kind == "COMPLEX":
                return "complex", "complex"
        return "other", "other"


def main():
    cls = Classes()
    summary, L = {}, ["# ppr8 PASS false positives by repeat context: d46 vs GRCh38-only\n",
                      "(`fp_context.py`; chr1, after the PoN, in the benchmark BED; seeds 0-2 summed. shared = the other "
                      f"arm has a PASS FP within {WINDOW} bp; other cand. = no FP there, but the other arm has a candidate "
                      f"left by the PoN within {WINDOW} bp; own only = the other arm has no candidate there. Context as "
                      "E17a: HP>=7, HP4-6, STR (unit 2-6, >= 3 units), VNTR, none. Cells: count (% of the row).)\n"]
    for s, p in SETS:
        counts = defaultdict(Counter)  # (arm, group) -> context
        types = defaultdict(Counter)
        for x in SEEDS:
            fp, cand, raw = {}, {}, {}
            for arm in ("d46", "linear"):
                v = RUNS / f"{s}_{p}_INDEL_{arm}_ppr8{x}" / "vcf_chr1" / "eval_bed" / "vcfeval"
                fp[arm] = records(v / "pon_calls" / "fp.vcf.gz")
                cand[arm] = sorted(r[1] for r in records(v / "pon_curve_input.vcf.gz"))
                raw[arm] = sorted(r[1] for r in records(v / "raw_curve_input.vcf.gz"))
            for arm, other in (("d46", "linear"), ("linear", "d46")):
                other_fp = sorted(r[1] for r in fp[other])
                for r in fp[arm]:
                    g = "shared" if near(r[1], other_fp) else "other cand." if near(r[1], cand[other]) else "own only"
                    t, c = cls(*r)
                    counts[(arm, g)][c] += 1
                    counts[(arm, "all")][c] += 1
                    types[(arm, g)][t] += 1
                    if g == "own only":
                        counts[(arm, OWN[0] if near(r[1], raw[other]) else OWN[1])][c] += 1
        v = RUNS / f"{s}_{p}_INDEL_d46_ppr8" / "vcf_chr1" / "eval_bed" / "vcfeval" / "pon_calls"
        truth = Counter(cls(*r)[1] for r in records(v / "tp-baseline.vcf.gz") + records(v / "fn.vcf.gz"))
        L += [f"\n## {s}T {p}\n", "| FPs | n | " + " | ".join(CONTEXTS) + " | INS / DEL |",
              "|---" * (len(CONTEXTS) + 3) + "|"]
        for arm in ("d46", "linear"):
            for g in ("all",) + GROUPS + OWN:
                c = counts[(arm, g)]
                n = sum(c.values())
                ty = sum((types[(arm, h)] for h in GROUPS), Counter()) if g == "all" else types.get((arm, g), Counter())
                label = f"{'d46' if arm == 'd46' else 'GRCh38-only'} {g}"
                L.append(f"| {label} | {n} | " + " | ".join(f"{c[k]} ({100 * c[k] / n:.0f}%)" if n else "-" for k in CONTEXTS)
                         + (f" | {ty['INS']} / {ty['DEL']} |" if ty else " | |"))
                summary[f"{s}_{p}/{arm}/{g}"] = dict(n=n, context=dict(c), type=dict(ty))
        n = sum(truth.values())
        L.append(f"| truth INDELs (chr1, BED) | {n} | " + " | ".join(f"{truth[k]} ({100 * truth[k] / n:.0f}%)" for k in CONTEXTS)
                 + " | |")
        summary[f"{s}_{p}/truth"] = dict(n=n, context=dict(truth))
    (HERE / "tables.md").write_text("\n".join(L) + "\n")
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

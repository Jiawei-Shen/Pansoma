"""chr1 INDEL calls of pansoma_net_v2 runs under each PoN rule of pon_features.py (HG008T Illumina, chr1 in the BED).

Calls and their TP/FP status: <run>/vcf_chr1/eval_bed/vcfeval/raw_curve/{tp,fp}.vcf.gz (every call with its score,
rtg vcfeval against the 587 chr1 truth records in the BED). A rule removes the calls it tags. Reported per rule:
- at the run's PASS threshold (validation recall 0.9): TP, FP, P, R, F1;
- the best F1 over all score thresholds after the rule.
Recall = kept TP calls / 587 (rtg's TP calls and TP baseline are equal in these runs).
Writes calls_chr1.md and calls_chr1.tsv (one row per call: run, status, PASS, score, features).
"""
import sys
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pon_features import FIELDS, PON_DIR, PON_FILES, RULES, as_dict, features, pon_record  # noqa: E402

pysam.set_verbosity(0)
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
RUN_NAMES = ["HG008_Illumina_INDEL_nopartial_b1024", "HG008_Illumina_INDEL_nopartial_b1024_nopc"]
TRUTH = 587
pons = [pysam.VariantFile(f"{PON_DIR}/{f}") for f in PON_FILES]


def feat(chrom, pos, ref, alt):
    recs = {}
    for p, pon in enumerate(pons):
        contig = chrom if chrom in pon.header.contigs else chrom[3:]
        recs[p] = [pon_record(r, p) for r in pon.fetch(contig, pos - 1, pos) if r.pos == pos]
    return features(recs, ref.upper(), alt.upper())


def f1(p, r):
    return 2 * p * r / (p + r) if p + r else 0.0


rows, md = [], ["# chr1 INDEL calls (HG008T Illumina, chr1 in the BED, 587 truth records) under each INDEL PoN rule\n"]
for run in RUN_NAMES:
    base = RUNS / run / "vcf_chr1/eval_bed/vcfeval/raw_curve"
    calls = []
    for status in ("tp", "fp"):
        with pysam.VariantFile(str(base / f"{status}.vcf.gz")) as vcf:
            for v in vcf:
                calls.append(dict(status=status, passed="PASS" in v.filter.keys(), score=float(v.qual), chrom=v.chrom,
                                  pos=v.pos, ref=v.ref, alt=v.alts[0], feat=feat(v.chrom, v.pos, v.ref, v.alts[0])))
    for c in calls:
        rows.append([run, c["status"], c["passed"], c["score"], c["chrom"], c["pos"], c["ref"], c["alt"], *c["feat"]])
    table = []
    for name, rule in [("no PoN", lambda f: False)] + list(RULES.items()):
        kept = [c for c in calls if not rule(as_dict(c["feat"]))]
        tp = sum(c["status"] == "tp" and c["passed"] for c in kept)
        fp = sum(c["status"] == "fp" and c["passed"] for c in kept)
        p = tp / (tp + fp) if tp + fp else 0.0
        r = tp / TRUTH
        best = (0.0, None, 0.0, 0.0)
        ntp = nfp = 0
        for c in sorted(kept, key=lambda c: -c["score"]):
            ntp += c["status"] == "tp"
            nfp += c["status"] == "fp"
            bp, br = ntp / (ntp + nfp), ntp / TRUTH
            if f1(bp, br) > best[0]:
                best = (f1(bp, br), c["score"], bp, br)
        tagged_tp = sum(c["status"] == "tp" for c in calls) - sum(c["status"] == "tp" for c in kept)
        tagged_fp = sum(c["status"] == "fp" for c in calls) - sum(c["status"] == "fp" for c in kept)
        table.append([name, f"{tagged_tp:,} / {tagged_fp:,}", f"{tp:,}", f"{fp:,}", f"{p:.3f}", f"{r:.3f}",
                      f"**{f1(p, r):.3f}**", f"{best[0]:.3f} (P {best[2]:.3f} R {best[3]:.3f})"])
    total_tp = sum(c["status"] == "tp" for c in calls)
    total_fp = sum(c["status"] == "fp" for c in calls)
    md.append(f"## {run}\n")
    md.append(f"All scored calls: TP {total_tp:,}, FP {total_fp:,}; PASS = the run's validation recall-0.9 threshold.\n")
    md.append("| rule | tagged TP / FP (all calls) | chr1 PASS TP | chr1 PASS FP | P | R | F1 | chr1 best F1 over thresholds |")
    md.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    md += ["| " + " | ".join(r) + " |" for r in table]
    md.append("")
(HERE / "calls_chr1.md").write_text("\n".join(md) + "\n")
with open(HERE / "calls_chr1.tsv", "w") as out:
    out.write("\t".join(["run", "status", "pass", "score", "chrom", "pos", "ref", "alt", *FIELDS]) + "\n")
    out.writelines("\t".join(map(str, r)) + "\n" for r in rows)
print((HERE / "calls_chr1.md").read_text())

"""Where the COLO829T chr1 truth alleles end up in a DeepSomatic tumor-only run (v1.9.0, PACBIO_TUMOR_ONLY model on the
COLO829T fiberseq BAM, default PoN filtering, chr1 only; run 2025-08 in COLO829T_fiberseq/deeepsomatic_TO_results).

Per truth allele (exact POS/REF/ALT): DeepSomatic FILTER (PASS / GERMLINE = the model's non-somatic class / PON / RefCall)
or no record; whether the repo PoN rule tags it (../truth_pon_COLO829T.tsv); median DeepSomatic VAF, truth VAF_PB and
population AF (max of gnomAD and CoLoRSdb exact-allele AF) per group. Writes deepsomatic_colo829t.txt.
"""
import csv, statistics, sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

pysam.set_verbosity(0)
csv.field_size_limit(sys.maxsize)
HERE = Path(__file__).resolve().parent
DS = "/scratch/jshen/data/COLO829T/COLO829T_fiberseq/deeepsomatic_TO_results/COLO829T_Pacbio_deepsomatic.vcf.gz"
TRUTH = "/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz"
info = {(r["chrom"], r["vcf_pos"], r["vcf_ref"], r["vcf_alt"]): r
        for r in csv.DictReader(open(HERE.parent / "truth_pon_COLO829T.tsv"), delimiter="\t")}
calls, bypos = {}, defaultdict(list)
with pysam.VariantFile(DS) as f:
    for v in f:
        flt, vaf = ";".join(v.filter.keys()), v.samples[0].get("VAF")
        for i, a in enumerate(v.alts or ()):
            calls[(v.chrom, v.pos, v.ref, a)] = (flt, vaf[i] if vaf else None)
        bypos[(v.chrom, v.pos)].append(flt)
groups = defaultdict(list)
with pysam.VariantFile(TRUTH) as f:
    for v in f.fetch("chr1"):
        kind = "SNV" if len(v.ref) == 1 and len(v.alts[0]) == 1 else "INDEL"
        k = (v.chrom, v.pos, v.ref, v.alts[0])
        if k in calls:
            st, vaf = calls[k]
        else:
            st, vaf = ("same POS, other allele" if (v.chrom, v.pos) in bypos else "no record"), None
        r = info[(v.chrom, str(v.pos), v.ref, v.alts[0])]
        pop = max([float(x) for x in (r["gnomAD_AF"], r["CoLoRSdb_AF"]) if x] or [0.0])
        groups[(kind, st)].append((vaf, float(v.info.get("VAF_PB") or 0), pop, r["rule"] != "none"))
lines = ["COLO829T chr1 truth in DeepSomatic tumor-only (PacBio model, fiberseq BAM, default PoN)", ""]
for kind in ("SNV", "INDEL"):
    n = sum(len(x) for (k, _), x in groups.items() if k == kind)
    lines.append(f"{kind}: {n:,} truth alleles")
    for (k, st), xs in sorted(groups.items(), key=lambda g: -len(g[1])):
        if k != kind:
            continue
        vafs = [x[0] for x in xs if x[0] is not None]
        med = lambda a: statistics.median(a) if a else float("nan")  # noqa: E731
        lines.append(f"  {st:24s} {len(xs):6,} ({100 * len(xs) / n:5.1f}%)  repo-PoN tagged {sum(x[3] for x in xs):5,}  "
                     f"median DeepSomatic VAF {med(vafs):.2f}  truth VAF_PB {med([x[1] for x in xs]):.2f}  "
                     f"population AF {med([x[2] for x in xs]):.1e} (>= 0.01: {sum(x[2] >= 0.01 for x in xs):,})")
    lines.append("")
(HERE / "deepsomatic_colo829t.txt").write_text("\n".join(lines))
print("\n".join(lines))

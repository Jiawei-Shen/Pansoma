"""INDEL_base chr1 calls at predict's best-F1 threshold (vcf_chr1/, rtg eval_bed): TP / FP by repeat context and AF,
FP by label reason; plus the training positives by label reason. Read-only."""
import collections
import gzip
import json

import pysam

R = "/scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_INDEL_base"
FASTA = pysam.FastaFile("/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
T = "/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors/INDEL"


def period(s):
    for p in range(1, len(s) + 1):
        if len(s) % p == 0 and s[:p] * (len(s) // p) == s:
            return s[:p]
    return s


def context(chrom, pos, ref, alt):
    ev = alt[1:] if len(alt) > len(ref) else ref[1:]
    unit = period(ev.upper())
    after = FASTA.fetch(chrom, pos, pos + 400).upper()
    n = 0
    while after[n * len(unit):(n + 1) * len(unit)] == unit:
        n += 1
    return "hp>=5" if len(unit) == 1 and n >= 5 else "STR>=3" if len(unit) > 1 and n >= 3 else "other"


def af_bin(af):
    return "<0.15" if af < 0.15 else "0.15-0.3" if af < 0.3 else "0.3-0.6" if af < 0.6 else ">=0.6"


pred = {}
for line in gzip.open(f"{R}/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.INDEL.predictions.ndjson.gz", "rt"):
    p = json.loads(line)
    pred[p["candidate_id"]] = p
C = collections.Counter
for name in ("tp", "fp"):
    ctx, afs, reasons, both = C(), C(), C(), C()
    for line in gzip.open(f"{R}/vcf_chr1/eval_bed/vcfeval/raw_calls/{name}.vcf.gz", "rt"):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        af = float(fmt["AF"])
        c = context(f[0], int(f[1]), f[3], f[4])
        ctx[c] += 1
        afs[af_bin(af)] += 1
        both[(c, af_bin(af))] += 1
        ids = f[2].split(";")
        reasons[pred[ids[0]]["reason"] if ids[0] in pred else "?"] += 1
    print(f"== rtg {name} (PASS at predict's threshold, in BED): {sum(ctx.values())}")
    print("  context", dict(ctx))
    print("  AF", dict(afs))
    print("  context x AF", sorted(both.items()))
    print("  label reason", reasons.most_common(8))

# unplaced PASS: label reasons
un = C()
for line in gzip.open(f"{R}/vcf_chr1/HG008T_Illumina.INDEL.linear.unplaced.vcf.gz", "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    if f[6] == "PASS":
        un[pred[f[2]]["reason"]] += 1
print("== unplaced PASS label reasons", un.most_common(6))

# training positives (chr2-22) by reason and AF
pos = C()
for chrom in [f"chr{i}" for i in range(2, 23)]:
    with open(f"{T}/{chrom}_labels.ndjson") as f:
        for line in f:
            r = json.loads(line)
            if r.get("label") == 1:
                pos[r.get("reason")] += 1
print("== training somatic INDEL tensors by reason (chr2-22, current labels)", pos.most_common())

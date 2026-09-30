"""INDEL_base chr1, in BED: which PoN subsets remove germline false calls but keep somatic truths. Uses the rtg
curve run (--all-records, PASS+LowQual) of vcf_chr1/eval_bed and the PoN tags of the calls and of the truth.
Approximate: removing calls does not change the rtg status of the others. Read-only."""
import gzip

R = "/scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_INDEL_base/vcf_chr1/eval_bed"


def tags(path):
    out = {}
    for line in gzip.open(path, "rt"):
        if line[0] == "#":
            continue
        f = line.split("\t")
        info = dict(x.partition("=")[::2] for x in f[7].split(";"))
        out[(f[0], f[1], f[3], f[4])] = set(info.get("PANSOMA_PON", "").split(",")) - {""}
    return out


calls = tags(f"{R}/HG008T_Illumina.INDEL.linear.pon.vcf.gz")
truth = tags(f"{R}/truth.INDEL.pon.vcf.gz")
recs = []
for name in ("tp", "fp"):
    for line in gzip.open(f"{R}/vcfeval/raw_curve/{name}.vcf.gz", "rt"):
        if line[0] == "#":
            continue
        f = line.split("\t")
        recs.append((float(f[5]), name == "tp", calls[(f[0], f[1], f[3], f[4])]))
fn_base = sum(1 for line in gzip.open(f"{R}/vcfeval/raw_curve/fn.vcf.gz", "rt") if line[0] != "#")
tp_base = sum(1 for line in gzip.open(f"{R}/vcfeval/raw_curve/tp-baseline.vcf.gz", "rt") if line[0] != "#")
N = tp_base + fn_base
thr = -10 * __import__("math").log10(1 - 0.28608837723731995)
G, D, K, C = "PoN1_gnomAD", "PoN2_dbSNP", "PoN3_1000G", "PoN4_CoLoRSdb"
print(f"truths in BED {N}; truth tags: " + ", ".join(f"{n} {sum(n in t for t in truth.values())}" for n in (G, D, K, C)))
print(f"{'PoN subset':<28} {'truth cap':>9} | at best-F1 t: {'TP':>4} {'FP':>5} {'P':>6} {'R':>6} {'F1':>6} | best F1 (P, R)")
for subset in ((), (G,), (G, D), (G, K), (G, D, K), (C,), (G, D, K, C)):
    keep = [(s, tp) for s, tp, t in recs if not (t & set(subset))]
    cap = 1 - sum(bool(t & set(subset)) for t in truth.values()) / len(truth)
    tp = sum(1 for s, t in keep if t and s >= thr)
    fp = sum(1 for s, t in keep if not t and s >= thr)
    p, r = tp / max(tp + fp, 1), tp / N
    best, cum_tp, cum_fp = (0, 0, 0), 0, 0
    for s, t in sorted(keep, reverse=True):
        cum_tp += t
        cum_fp += not t
        pp, rr = cum_tp / (cum_tp + cum_fp), cum_tp / N
        f1 = 2 * pp * rr / (pp + rr) if pp + rr else 0
        best = max(best, (f1, pp, rr))
    print(f"{'+'.join(x.split('_')[1] for x in subset) or 'none':<28} {cap:>9.3f} | {tp:>17} {fp:>5} {p:>6.3f} {r:>6.3f} "
          f"{2 * p * r / (p + r) if p + r else 0:>6.3f} | {best[0]:.3f} ({best[1]:.3f}, {best[2]:.3f})")

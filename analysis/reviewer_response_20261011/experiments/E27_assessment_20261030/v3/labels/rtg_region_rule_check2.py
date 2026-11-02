#!/usr/bin/env python3
"""INDEL only: test alternative rtg -e in-region rules (anchor-trimmed span [s, e); INS zero-length at s)."""
import bisect, gzip, json
from collections import Counter
from pathlib import Path
R = Path("/scratch/jshen/data/pansoma_net_v2_runs/HapMap_Illumina_INDEL/vcf_chr1"); X = Path("/scratch/jshen/data/reviewer_experiments_20261013/E23_nograph")
bed = [tuple(map(int, l.split("\t")[1:3])) for l in open(X / "eval_chr1.bed")]; bs = [b[0] for b in bed]
def inb(x):
    i = bisect.bisect_right(bs, x) - 1
    return i >= 0 and bed[i][0] <= x < bed[i][1]
def recs(d):
    out = {}
    for n in ("tp-baseline.vcf.gz", "fn.vcf.gz"):
        for line in gzip.open(d / n, "rt"):
            if line[0] == "#": continue
            c = line.rstrip("\n").split("\t"); v = dict(zip(c[8].split(":"), c[9].split(":")))
            af = v.get("AF", "."); out[(c[0], int(c[1]), c[3], c[4])] = None if af in (".", "") else float(af)
    return out
F = {k for k, a in recs(R / "eval_bed/vcfeval/raw_calls").items() if a is not None and 0.0099 <= a <= 0.1651}
E = {k for k, a in recs(R / "eval_bed_nograph_all/vcfeval/raw_calls").items() if a is not None and 0.0099 <= a <= 0.1651}
rules = Counter(); tot = Counter()
for k in F:
    _, p, r, a = k
    n = 0
    while n < min(len(r), len(a)) and r[n] == a[n]: n += 1
    s, e = p - 1 + n, p - 1 + len(r)  # trimmed REF span
    kind = "INS" if len(a) > len(r) else "DEL" if len(a) < len(r) else "other"
    cand = {
        "first": inb(s),
        "any_trimmed": any(inb(x) for x in range(s, max(e, s + 1))),
        "all_trimmed": all(inb(x) for x in range(s, max(e, s + 1))),
        "any_full": any(inb(x) for x in range(p - 1, p - 1 + len(r))),
        "all_full": all(inb(x) for x in range(p - 1, p - 1 + len(r))),
        "first_and_last": inb(s) and inb(max(e, s + 1) - 1),
    }
    kept = k in E
    for name, v in cand.items():
        rules[(kind, name)] += (v == kept)
    tot[kind] += 1
print(dict(tot)); print({f"{k[0]}:{k[1]}": v for k, v in sorted(rules.items())})
ins = Counter(); bad = []
for k in F:
    _, p, r, a = k
    if len(a) <= len(r): continue
    n = 0
    while n < min(len(r), len(a)) and r[n] == a[n]: n += 1
    s = p - 1 + n
    kept = k in E
    for name, v in (("s-1 or s", inb(s - 1) or inb(s)), ("s-1 and s", inb(s - 1) and inb(s)), ("s-1", inb(s - 1)), ("s", inb(s))):
        ins[name] += (v == kept)
    if (inb(s)) != kept: bad.append((k, inb(s - 1), inb(s), kept))
print(dict(ins)); print(bad[:8])

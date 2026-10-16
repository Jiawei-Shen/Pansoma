#!/usr/bin/env python3
"""Sanity check of the simulation truth: for reads with no HG008-N indel / SV and no inserted bases, the read must have
exactly (ref_end - ref_start) = read length, and its mismatches against GRCh38[ref_start:ref_end] (on the truth strand)
must equal the haplotype's SNVs in the read + the simulated sequencing errors."""
import gzip
import sys

import numpy as np
import pysam

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import e11b_common as C  # noqa: E402

simdir, ref_fa = sys.argv[1], sys.argv[2]
T = dict(np.load(f"{simdir}/truth.npz"))
n = len(T["read"])
paired = (T["mate"] > 0).any()
seqs = {}
files = [f"{simdir}/sim_R1.fq.gz", f"{simdir}/sim_R2.fq.gz"] if paired else [f"{simdir}/sim.fq.gz"]
for fi, f in enumerate(files):
    with gzip.open(f, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].strip()
            elif i % 4 == 1:
                seqs[name] = line.strip()
fa = pysam.FastaFile(ref_fa)
REF = {c: fa.fetch(c).upper() for c in C.CHROMS}
ok = bad = skipped = 0
examples = []
for r in range(min(n, 20000)):
    if T["n_indel"][r] or T["n_sv"][r] or T["inserted_bases"][r]:
        skipped += 1
        continue
    name = f"i{T['read'][r]}/{T['mate'][r]}" if paired else f"h{T['read'][r]}"
    s = seqs[name]
    c = C.CHROMS[T["chrom"][r]]
    ref = REF[c][int(T["ref_start"][r]):int(T["ref_end"][r])]
    if T["ref_strand"][r] < 0:
        s = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    if len(ref) != len(s):
        bad += 1
        examples.append((name, "length", len(ref), len(s)))
        continue
    mm = sum(1 for x, y in zip(s, ref) if x != y)
    if mm == T["n_snv"][r] + T["n_err"][r]:
        ok += 1
    else:
        bad += 1
        if len(examples) < 10:
            examples.append((name, c, int(T["ref_start"][r]), mm, int(T["n_snv"][r]), int(T["n_err"][r])))
print(dict(checked=ok + bad, exact=ok, mismatch=bad, skipped_with_indel=skipped))
for e in examples[:10]:
    print(e)

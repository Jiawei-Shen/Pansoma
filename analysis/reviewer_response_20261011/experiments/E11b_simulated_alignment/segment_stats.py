#!/usr/bin/env python3
"""E11b: how much of GRCh38 chr1 / chr20 the one-to-one haplotype segments cover, and their SMaHT region composition
(the simulation can only sample where the assembly-to-GRCh38 truth is one-to-one). Writes results/segment_stats.json."""
import json
import os
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e11b_common as C  # noqa: E402

sim, ref_fa, out = sys.argv[1], sys.argv[2], sys.argv[3]
segs = C.read_tsv(f"{sim}/segments.tsv")
cums = C.smaht_cumsums()
fa = pysam.FastaFile(ref_fa)
res = {}
for c in C.CHROMS:
    seq = np.frombuffer(fa.fetch(c).upper().encode(), np.uint8)
    nonN = seq != ord("N")
    cd, cx = cums[c]
    code = np.diff(cd).astype(np.int8) + 2 * np.diff(cx).astype(np.int8)
    allc = {k: int(((code == i) & nonN).sum()) for i, k in enumerate(["easy", "difficult", "extreme"])}
    r = {"non_N_bases": int(nonN.sum()), "smaht_non_N": allc}
    for h in (1, 2):
        cov = np.zeros(len(seq), bool)
        for s in segs:
            if s["chrom"] == c and int(s["hap"]) == h:
                cov[int(s["bed_start"]):int(s["bed_end"])] = True
        r[f"hap{h}_covered_bases"] = int((cov & nonN).sum())
        r[f"hap{h}_covered_frac_of_non_N"] = round(float((cov & nonN).sum() / nonN.sum()), 4)
        r[f"hap{h}_smaht_covered_frac"] = {k: round(float(((code == i) & nonN & cov).sum() / max(allc[k], 1)), 4)
                                         for i, k in enumerate(["easy", "difficult", "extreme"])}
    res[c] = r
os.makedirs(out, exist_ok=True)
json.dump(res, open(f"{out}/segment_stats.json", "w"), indent=1)
print(json.dumps(res, indent=1))

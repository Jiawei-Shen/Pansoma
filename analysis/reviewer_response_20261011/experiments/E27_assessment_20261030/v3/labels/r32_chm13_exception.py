#!/usr/bin/env python3
"""R3.2 / Q4: the 'represented in the graph' stratum of Table R3.2 is E7's carrier-count proxy (>= 46 of the 463
non-GRCh38 HPRC v2.1 haplotypes of the vcfbub/vcfwave VCF carry the allele; CHM13 = bit 0, counted as one of the 463).
In the real d46 graph the CHM13 (and GRCh38) reference paths are never filtered, so a CHM13 allele with < 46 carriers
is still spelled by the CHM13 path. Count those, per sample x role x kind (PASS, in BED, chr1-22; E7 v2.1 joined.npz)."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E7_graph_composition")
import os; os.environ["E7_GRAPH"] = "v2.1"
from common import FINE, GROUPS, SAMPLES, ROLES  # noqa: E402
J = np.load("/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition/v2.1/joined.npz")
W = J["masks"]
cnt = sum(np.bitwise_count(W[:, w]).astype(np.int64) for w in range(W.shape[1]))
found = J["found"]
cnt[~found] = 0
chm13 = ((W[:, 0] & np.uint64(1)) == 1) & found
inbed = J["passed"] & J["in_bed"]
fine = J["fine"]
out = {}
for si, s in enumerate(SAMPLES[:2]):
    for ri, r in enumerate(ROLES):
        for g in ("SNV", "INDEL"):
            codes = [FINE.index(f) for f in GROUPS[g]]
            sel = (J["sample"] == si) & (J["role"] == ri) & inbed & np.isin(fine, codes)
            n = int(sel.sum()); pred = int((sel & (cnt >= 46)).sum())
            ch = int((sel & (cnt < 46) & chm13).sum()); any_ = int((sel & (cnt >= 1)).sum())
            out[f"{s} {r} {g}"] = dict(n=n, carriers_ge46=pred, share_ge46=round(pred / n, 4) if n else None,
                                       chm13_carrier_below46=ch, share_chm13_below46=round(ch / n, 4) if n else None,
                                       any_carrier=any_)
            print(s, r, g, out[f"{s} {r} {g}"], flush=True)
Path(__file__).with_suffix(".json").write_text(json.dumps(out, indent=1) + "\n")

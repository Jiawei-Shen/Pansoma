"""Somatic share of training candidates by read AF and population status (SNV site features), per set (training
chromosomes: every autosome except chr1 and chr20). Labels: 1 somatic, 2 germline, 0 non, -1 ignored."""
import numpy as np, sys
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/tmp/snv_sitefeat_20261103/code/machine_learning")
from pansoma_net_v2.data import KindIndex, SCALARS
from pansoma_net_v2 import site_features as sf
D = "/scratch/jshen/data/pansoma_v2_tensors"; C = "/scratch/jshen/data/pansoma_net_v2_runs/index_cache"
SF = "/scratch/jshen/data/snv_sitefeat_20261103/site_features"
for s in sys.argv[1:]:
    ix = KindIndex(f"{D}/{s}/tensors/SNV", cache_dir=C)
    d = sf.load(SF, ix); a = ix.arrays
    lab, af, pon = a["label"], a["scalars"][:, SCALARS.index("af")], d["pon_tag"]
    chrom = np.array(ix.meta["chroms"])[a["chrom"].astype(int)]
    train = ~np.isin(chrom, ["chr1", "chr20"]) & (lab >= 0)
    print(f"== {s}")
    for lo, hi in ((0.08, 0.4), (0.4, 0.6), (0.6, 0.85), (0.85, 1.01)):
        m = train & (af >= lo) & (af < hi)
        cells = []
        for name, sel in (("all", m), ("not in population", m & (pon == 0)), ("population", m & (pon == 1))):
            n = int(sel.sum()); som = int((lab[sel] == 1).sum())
            cells.append(f"{name}: {som:,} somatic of {n:,} ({som / max(n, 1):.2%})")
        print(f"  AF [{lo}, {hi}): " + " | ".join(cells), flush=True)

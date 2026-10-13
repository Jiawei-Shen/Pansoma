"""How the tumor (COLO829T) and normal (COLO829BL) reads of the in-silico mixtures differ, per read group.
Run in the container (pysam): ris/submit.sh -g 0 -n 2 -m 16 NAME py.sh $PWD/read_source.py [BAM ...]
Default BAMs: the Illumina (WGS) mixtures. For each BAM: the header's @RG and @PG lines, then the first N primary,
mapped reads of chr2:20,000,000-: per read group the read count, distinct BQ values (the alphabet), BQ mean and the most
common values, MAPQ mean and the fraction MAPQ 60, read length mode, and the fraction of reads with soft clips.
Writes read_source.json next to this script."""
import collections
import json
import sys
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
MIX = Path("/storage2/fs1/epigenome/Active/jiawei.s/Dataset/COLO829T_Tumor_Normal_inSilico_Mix")
bams = sys.argv[1:] or [str(MIX / "COLO829_WGS" / f"MIX_{p}_SMHTCOLO829.mixed.sorted.bam") for p in ("p75", "p50", "p25")]
N = 200_000
out = {}
for b in bams:
    f = pysam.AlignmentFile(b)
    h = f.header.to_dict()
    rgs = h.get("RG", [])
    pgs = [{k: v for k, v in pg.items() if k in ("ID", "PN", "VN", "CL")} for pg in h.get("PG", [])]
    print(f"== {b}\n  @RG {len(rgs)}: " + "; ".join(f"{r.get('ID')} SM={r.get('SM')} PL={r.get('PL')} PU={r.get('PU', '')}" for r in rgs[:12]))
    for pg in pgs[:12]:
        print(f"  @PG {pg.get('ID')} {pg.get('PN', '')} {pg.get('VN', '')} {str(pg.get('CL', ''))[:160]}")
    st = collections.defaultdict(lambda: {"n": 0, "bq": collections.Counter(), "mapq": [], "len": collections.Counter(), "clip": 0})
    n = 0
    for r in f.fetch("chr2", 20_000_000):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_qualities is None:
            continue
        rg = r.get_tag("RG") if r.has_tag("RG") else "none"
        s = st[rg]; s["n"] += 1
        s["bq"].update(r.query_qualities.tolist() if hasattr(r.query_qualities, "tolist") else list(r.query_qualities))
        s["mapq"].append(r.mapping_quality); s["len"][r.query_length] += 1
        s["clip"] += any(op == 4 for op, _ in (r.cigartuples or []))
        n += 1
        if n >= N:
            break
    res = {}
    for rg, s in sorted(st.items(), key=lambda x: -x[1]["n"]):
        tot = sum(s["bq"].values()); mean = sum(q * c for q, c in s["bq"].items()) / tot
        res[rg] = {"reads": s["n"], "bq_distinct": len(s["bq"]), "bq_mean": round(mean, 2),
                   "bq_top": [(int(q), round(c / tot, 3)) for q, c in s["bq"].most_common(5)],
                   "mapq_mean": round(sum(s["mapq"]) / len(s["mapq"]), 2),
                   "mapq60": round(sum(m == 60 for m in s["mapq"]) / len(s["mapq"]), 3),
                   "len_mode": s["len"].most_common(1)[0][0], "softclip_frac": round(s["clip"] / s["n"], 3)}
        x = res[rg]
        print(f"  RG {rg}: reads {x['reads']:,} | BQ distinct {x['bq_distinct']}, mean {x['bq_mean']}, top {x['bq_top']} | "
              f"MAPQ mean {x['mapq_mean']}, =60 {x['mapq60']} | len {x['len_mode']} | soft-clipped {x['softclip_frac']}", flush=True)
    out[b] = {"RG": rgs, "PG": pgs, "per_read_group": res}
(HERE / "read_source.json").write_text(json.dumps(out, indent=1) + "\n")

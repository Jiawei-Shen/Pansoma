"""v2 false calls by AF on chr1 (same scorer as my_check.py): composition at matched TP, and what an AF floor does.
Read-only on inputs; run from machine_learning (imports my_check, which recomputes its curves in ~30 s)."""
import collections, io, contextlib, sys
import numpy as np
sys.path.insert(0, "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/synthesis")
with contextlib.redirect_stdout(io.StringIO()):
    import my_check as M
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning")
from pansoma_net_v2.data import KindIndex, SCALARS  # noqa: E402

ix = KindIndex(M.V2 + "/SNV", "/scratch/jshen/data/pansoma_net_v2_runs/index_cache"); a = ix.arrays
chr1 = np.flatnonzero(a["chrom"] == ix.meta["chroms"].index("chr1"))
AFT = a["scalars"][chr1, SCALARS.index("af")]
import gzip, json  # noqa: E402
SITE_AF = {}
with open(M.V2 + "/SNV/chr1_labels.ndjson") as fl:
    for k, line in enumerate(fl):
        g = json.loads(line)["grch38"]
        if g is not None:
            key = (g["pos0"], g["ref"].upper(), g["alt"].upper())
            SITE_AF[key] = max(SITE_AF.get(key, 0.0), float(AFT[k]))
BANDS = [0.0, 0.08, 0.10, 0.15, 0.20, 0.35, 0.60, 1.01]
def band(x):
    for lo, hi in zip(BANDS, BANDS[1:]):
        if lo <= x < hi:
            return f"{lo:.2f}-{hi if hi <= 1 else 1:.2f}"
def fp_sites(sites, thr, pon=None):
    return [k for k, s in sites.items() if s >= thr and k not in M.TRUTH and M.inconf(k[0])
            and not (pon is not None and pon.get(k, False))]
def table(name, keys):
    c = collections.Counter((band(SITE_AF.get(k, -1)) if k in SITE_AF else "no v2 AF", "germ" if k in M.GERM else "non-germ") for k in keys)
    bands = sorted({b for b, _ in c}, key=lambda s: (s == "no v2 AF", s))
    tot = collections.Counter(g for _, g in c.elements())
    print(f"  {name}: {len(keys)} false calls (germline {tot['germ']}, non-germline {tot['non-germ']})")
    for b in bands:
        print(f"      AF {b:>9}: germline {c[(b, 'germ')]:>5}   non-germline {c[(b, 'non-germ')]:>5}")

print("==== false calls by AF (v2 AF of the site's tensor; v1 calls looked up in v2's AF where the site has a v2 tensor) ====")
for run in ("HG008_Illumina_SNV_base", "HG008_Illumina_SNV_keephard_w100", "HG008_Illumina_SNV_allnon_w100"):
    s, o, nt, neu = M.V2R[run]
    r = run.replace("HG008_Illumina_SNV_", "v2_")
    t0 = M.at_tp(M.curve(s), 608); t1 = M.at_tp(M.curve(s, pon=M.PONF), 561)
    table(f"{r} no PoN @ TP 608 (p>={t0[0]:.3f})", fp_sites(s, t0[0]))
    table(f"{r} repo PoN @ TP {t1[1]} (p>={t1[0]:.3f})", fp_sites(s, t1[0], M.PONF))
table("v1 e053 no PoN (PROB>=0.5, TP 608)", fp_sites(M.V1["v1_e053"], 0.5))
table("v1 e053 repo PoN (TP 561)", fp_sites(M.V1["v1_e053"], 0.5, M.PONF))

print("\n==== an AF floor on v2 calls (sites below the floor dropped; truths lost counted) ====")
for run in ("HG008_Illumina_SNV_base", "HG008_Illumina_SNV_keephard_w100", "HG008_Illumina_SNV_allnon_w100"):
    s, o, nt, neu = M.V2R[run]
    r = run.replace("HG008_Illumina_SNV_", "v2_")
    for floor in (0.0, 0.08, 0.10, 0.15, 0.20, 0.25):
        sf = {k: v for k, v in s.items() if SITE_AF.get(k, 1.0) >= floor}
        c0, c1 = M.curve(sf), M.curve(sf, pon=M.PONF)
        m0, m1 = M.summ(c0), M.summ(c1)
        a0, a1 = M.at_tp(c0, 608), M.at_tp(c1, 561)
        p0 = f"{a0[1] / (a0[1] + a0[2]):.3f}" if a0 else "  n.a."; p1 = f"{a1[1] / (a1[1] + a1[2]):.3f}" if a1 else "  n.a."
        print(f"  {r:<18} floor {floor:.2f}: ceiling {m0['ceiling']:.3f} | no PoN: P@TP608 {p0}, max F1 {m0['maxF1'][0]:.3f} "
              f"| repo PoN: P@TP561 {p1}, max F1 {m1['maxF1'][0]:.3f}")
print("v1 e053 for reference: no PoN P@TP608 0.100 (max F1 0.183); own PoN P@TP561 0.522 (max F1 0.643)")

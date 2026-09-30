"""Where v2's non-germline false calls (after the repo PoN, at TP 561 on chr1) come from: v2 label/reason, AF, read
features, distance to the nearest somatic truth / germline INDEL, and whether v1 had a candidate there and its score.
Read-only on inputs. Run from machine_learning."""
import bisect, collections, contextlib, gzip, io, json, sys
import numpy as np
import pysam
sys.path.insert(0, "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/synthesis")
with contextlib.redirect_stdout(io.StringIO()):
    import my_check as M
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning")
from pansoma_net_v2.data import KindIndex, TensorDataset, SCALARS  # noqa: E402

ix = KindIndex(M.V2 + "/SNV", "/scratch/jshen/data/pansoma_net_v2_runs/index_cache"); a = ix.arrays
chr1 = np.flatnonzero(a["chrom"] == ix.meta["chroms"].index("chr1"))
lab = [json.loads(l) for l in open(M.V2 + "/SNV/chr1_labels.ndjson")]
assert len(lab) == len(chr1)
by_site = collections.defaultdict(list)
for k, r in enumerate(lab):
    g = r["grch38"]
    if g is not None:
        by_site[(g["pos0"], g["ref"].upper(), g["alt"].upper())].append(k)

# positions of somatic truth (all kinds) and germline INDELs / SNVs on chr1
som_all, germ_indel, germ_snv = [], [], []
for r in pysam.VariantFile(M.TRUTH_ALL).fetch("chr1"):
    som_all.append((r.pos - 1, "SNV" if len(r.ref) == 1 and all(len(x) == 1 for x in r.alts) else "INDEL"))
seen = False
for r in pysam.VariantFile(M.GVCF):
    if r.chrom != "chr1":
        if seen: break
        continue
    seen = True
    (germ_snv if len(r.ref) == 1 and all(len(x) == 1 for x in (r.alts or ())) else germ_indel).append(r.pos - 1)
som_pos = sorted(p for p, _ in som_all); som_indel = sorted(p for p, k in som_all if k == "INDEL")
germ_indel.sort(); germ_snv.sort()
def dist(sorted_pos, p, skip_self=False):
    i = bisect.bisect_left(sorted_pos, p); best = 10 ** 9
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(sorted_pos) and not (skip_self and sorted_pos[j] == p):
            best = min(best, abs(sorted_pos[j] - p))
    return best

# v1 test rows -> GRCh38 keys and e067 probability
V1TEST = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_testing_data_SNV/5ch_testing_data_SNV_chr1/variant_summary.ndjson"
V1VAL = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset/variant_summary_classified.ndjson"
val = {}
for line in open(V1VAL):
    r = json.loads(line); val[(r["node_id"], r["variant_key"])] = (r["genomic_position"] - 1, r["v_ref"].upper(), r["v_alt"].upper())
p1 = np.load("/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head/v1_full/v1_chr1_test_p_true.e067.npy")
v1p = {}
for i, line in enumerate(open(V1TEST)):
    r = json.loads(line); s = val.get((r["node_id"], r["variant_key"]))
    if s is not None:
        v1p[s] = max(v1p.get(s, 0.0), float(p1[i]))

def feats(ks):
    """per tensor: A1 reads, A1 strand fraction (forward), mean A1 BQ at the site, mean A1 MAPQ, mismatches per A1 read
    vs per REF read (read base != graph base over covered columns, site column excluded)."""
    ds = TensorDataset([(ix, chr1[np.asarray(ks)])])
    out = []
    for j in range(len(ds)):
        x, blocks, _, _ = ds[j]; x = x.numpy().astype(np.int16); a1, alt, ref, n = blocks.tolist()
        cov = x[4] > 0; diff = cov & (x[0] != x[5]); diff[:, 50] = False
        mm = diff.sum(1)
        out.append(dict(a1=a1, fwd=float((x[7, :a1, 50] == 1).mean()) if a1 else np.nan,
                        bq=float(x[1, :a1, 50].mean()) if a1 else np.nan, mq=float(x[3, :a1, 50].mean()) if a1 else np.nan,
                        mm_a1=float(mm[:a1].mean()) if a1 else np.nan, mm_ref=float(mm[alt:ref].mean()) if ref > alt else np.nan))
    return out

AF = SCALARS.index("af")
for run in ("HG008_Illumina_SNV_base", "HG008_Illumina_SNV_allnon_w100"):
    s, o, nt, neu = M.V2R[run]
    thr = M.at_tp(M.curve(s, pon=M.PONF), 561)[0]
    fp = [k for k, v in s.items() if v >= thr and k not in M.TRUTH and M.inconf(k[0]) and not M.PONF.get(k, False) and k not in M.GERM]
    tp = [k for k, v in s.items() if v >= thr and k in M.TRUTH and not M.PONF.get(k, False)]
    print(f"\n==== {run.replace('HG008_Illumina_SNV_', 'v2_')}: {len(fp)} non-germline false calls after the repo PoN (p >= {thr:.3f}); {len(tp)} true calls ====")
    ks = [max(by_site[k], key=lambda t: a["scalars"][chr1[t], AF]) for k in fp]
    lab_c = collections.Counter((lab[t]["label"], lab[t]["reason"]) for t in ks)
    print("  v2 label / reason of the called tensor:")
    for (l, r), n in lab_c.most_common():
        print(f"     label {l:>2}  {r:<40} {n}")
    afs = np.array([a["scalars"][chr1[t], AF] for t in ks])
    print(f"  AF: median {np.median(afs):.3f}; < 0.20: {int((afs < 0.2).sum())}")
    d_si = np.array([dist(som_indel, k[0]) for k in fp]); d_s = np.array([dist(som_pos, k[0]) for k in fp])
    d_gi = np.array([dist(germ_indel, k[0]) for k in fp]); d_gs = np.array([dist(germ_snv, k[0], True) for k in fp])
    for name, d in (("somatic INDEL truth", d_si), ("any somatic truth", d_s), ("germline INDEL (HG008-N)", d_gi), ("other germline SNV", d_gs)):
        print(f"  nearest {name:<26}: <=10 bp {int((d <= 10).sum()):>4}  <=50 bp {int((d <= 50).sum()):>4}  <=200 bp {int((d <= 200).sum()):>4}")
    lowaf = afs < 0.2
    near = (d_gi <= 50) | (d_si <= 50)
    print(f"  AF < 0.20 and within 50 bp of a germline or somatic INDEL: {int((lowaf & near).sum())} of {int(lowaf.sum())} low-AF")
    v1c = np.array([k in v1p for k in fp]); v1s = np.array([v1p.get(k, np.nan) for k in fp])
    print(f"  v1 has a candidate at the site: {int(v1c.sum())} of {len(fp)}; v1 p_true >= 0.5 there: {int((v1s >= 0.5).sum())}; "
          f"v1 median p where present: {np.nanmedian(v1s):.3f}")
    F = feats(ks); T = feats([max(by_site[k], key=lambda t: s.get(k, 0)) for k in tp[:300]])
    def med(L, key): return np.nanmedian([f[key] for f in L])
    print("  read features (median)          false calls    true calls")
    for key, name in (("a1", "A1 reads"), ("fwd", "A1 forward-strand fraction"), ("bq", "A1 base quality at site"),
                      ("mq", "A1 MAPQ"), ("mm_a1", "other mismatches per A1 read"), ("mm_ref", "other mismatches per REF read")):
        print(f"     {name:<30} {med(F, key):>10.2f}   {med(T, key):>10.2f}")
    fs = np.array([f["fwd"] for f in F]); print(f"  A1 reads all on one strand (fraction 0 or 1): false calls {int(((fs == 0) | (fs == 1)).sum())} of {len(F)}; "
          f"true calls {int(sum(1 for f in T if f['fwd'] in (0.0, 1.0)))} of {len(T)}")

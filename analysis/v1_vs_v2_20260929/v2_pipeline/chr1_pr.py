"""Truth-level PR of the v2 SNV_base chr1 predictions, with the false calls broken down by what they are."""
import collections, csv, gzip, json, sys
import numpy as np
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning")
from pansoma_net_v2.bed import Bed, in_region

D = "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors"
S = D + "/SNV"
RUNS = "/scratch/jshen/data/pansoma_net_v2_runs"
man = json.load(open(S + "/labels.manifest.json"))
conf = Bed.read(man["truth"]["somatic"]["bed"]).intersect(Bed.read(man["truth"]["germline"]["bed"]))

tt = set()
with open(D + "/somatic.recall.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["chrom"] == "chr1" and r["kind"] == "SNP" and r["passed"] == "True" and r["in_bed"] == "True":
            tt.add(int(r["truth_id"]))
assert len(tt) == 697

labs = []
with open(S + "/chr1_labels.ndjson") as fl:
    for line in fl:
        lab = json.loads(line)
        g = lab["grch38"]
        e = 0 if lab["reason"] == "off_reference_no_truth_match" else lab["label"]
        if not in_region(conf, "SNP", g, lab.get("anchor")):
            e = -1
        ids = set()
        if lab["label"] == 1:
            ids = {d["truth_id"] for d in lab.get("somatic") or [] if d.get("representative")}
            if lab.get("partial_truth"):
                ids.add(lab["partial_truth"]["truth_id"])
        key = (g["pos0"], g["alt"]) if g else ("off", lab["candidate_id"])
        labs.append((lab["candidate_id"], e, ids, key, g is None))


def curve(run):
    p = np.empty(len(labs), np.float32)
    with gzip.open(f"{RUNS}/{run}/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz", "rt") as f:
        for k, line in enumerate(f):
            r = json.loads(line)
            assert r["candidate_id"] == labs[k][0]
            p[k] = r["p_somatic"]
    best = {}
    neg = collections.defaultdict(dict)  # kind -> key -> max p
    for k, (cid, e, ids, key, off) in enumerate(labs):
        if e < 0:
            continue
        if e == 1:
            own = ids & tt
            if own:
                for t in own:
                    best[t] = max(best.get(t, -1), p[k])
            else:  # other-kind or non-697 truth
                neg["other_kind_label1"][key] = max(neg["other_kind_label1"].get(key, -1), p[k])
        else:
            kind = "germline" if e == 2 else ("non_offref" if off else "non_onref")
            neg[kind][key] = max(neg[kind].get(key, -1), p[k])
    pos = np.array(sorted(best.values()))
    negs = {k: np.array(list(v.values())) for k, v in neg.items()}
    print(f"\n== {run}: truths with a test tensor {len(pos)}/697 (ceiling {len(pos)/697:.3f}); negatives",
          {k: len(v) for k, v in negs.items()})
    thr = np.unique(np.concatenate([pos, *negs.values()]))[::-1]
    def at(t, fp_kinds):
        tp = int((pos >= t).sum())
        fp = sum(int((negs[k] >= t).sum()) for k in fp_kinds if k in negs)
        return tp, fp
    variants = {
        "strict (FP = non + germline + other-kind label-1 calls)": ("non_onref", "non_offref", "germline", "other_kind_label1"),
        "v2-default-like (other-kind calls not FP)": ("non_onref", "non_offref", "germline"),
        "no germline FP (hypothetical: patient germline absent)": ("non_onref", "non_offref", "other_kind_label1"),
        "germline FP only": ("germline",),
    }
    # vectorised over thresholds
    for name, kinds in variants.items():
        ts = thr
        tps = len(pos) - np.searchsorted(pos, ts, side="left")
        fps = np.zeros(len(ts), int)
        for k in kinds:
            if k in negs:
                s = np.sort(negs[k])
                fps += len(s) - np.searchsorted(s, ts, side="left")
        rec = tps / 697
        prec = np.where(tps + fps > 0, tps / np.maximum(1, tps + fps), 1.0)
        f1 = np.where(prec + rec > 0, 2 * prec * rec / np.maximum(1e-12, prec + rec), 0)
        i = int(np.argmax(f1))
        out = [f"maxF1 {f1[i]:.3f} (P {prec[i]:.3f} R {rec[i]:.3f} TP {tps[i]} FP {fps[i]} thr {ts[i]:.4f})"]
        for P in (0.05, 0.10, 0.15, 0.20, 0.30):
            ok = prec >= P
            out.append(f"R@P{P:.2f}={rec[ok].max():.3f}" if ok.any() else f"R@P{P:.2f}=na")
        for R in (0.5, 0.7, 0.8, 0.9):
            ok = rec >= R
            out.append(f"P@R{R:.1f}={prec[ok].max():.3f}" if ok.any() else f"P@R{R:.1f}=na (max R {rec.max():.3f})")
        print(f"  {name}:\n    " + "  ".join(out))
    # FP composition at a few recalls (strict)
    for R in (0.5, 0.7, 0.8, 0.9):
        t = np.sort(pos)[::-1][int(np.ceil(R * 697)) - 1] if int(np.ceil(R * 697)) <= len(pos) else None
        if t is None:
            continue
        comp = {k: int((v >= t).sum()) for k, v in negs.items()}
        print(f"  at truth recall {R}: threshold {t:.4f}, TP {int((pos >= t).sum())}, FP by kind {comp}")


for run in sys.argv[1:]:
    curve(run)

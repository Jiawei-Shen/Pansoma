"""Per-tensor (v1-style) metrics of the v2 chr1 predictions: positives = label-1 tensors whose representative
allele is a somatic SNP truth (exact), in the test region; variants of the negative set."""
import csv, gzip, json, sys
import numpy as np
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning")
from pansoma_net_v2.bed import Bed, in_region

D = "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors"
S = D + "/SNV"
RUNS = "/scratch/jshen/data/pansoma_net_v2_runs"
OUT = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/v2_pipeline"
RUN_NAMES = ["HG008_Illumina_SNV_base", "HG008_Illumina_SNV_scalars", "HG008_Illumina_SNV_keephard_w100",
             "HG008_Illumina_SNV_allnon_w100"]
man = json.load(open(S + "/labels.manifest.json"))
conf = Bed.read(man["truth"]["somatic"]["bed"]).intersect(Bed.read(man["truth"]["germline"]["bed"]))
cls = []  # 'pos', 'partial', 'germ', 'non', 'offnon', None
with open(S + "/chr1_labels.ndjson") as fl:
    for line in fl:
        lab = json.loads(line)
        g = lab["grch38"]
        if not in_region(conf, "SNP", g, lab.get("anchor")):
            cls.append(None); continue
        if lab["reason"] == "off_reference_no_truth_match":
            cls.append("offnon")
        elif lab["label"] == 1:
            cls.append("pos" if lab["reason"] == "representative_allele_is_somatic_truth" else "partial")
        elif lab["label"] == 2:
            cls.append("germ")
        elif lab["label"] == 0:
            cls.append("non")
        else:
            cls.append(None)
cls = np.array(cls, dtype=object)
pon = {}
with open(OUT + "/chr1_pon_annotations.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        pon[int(r["index"])] = any(r[k] == "1" for k in ("gnomad", "dbsnp", "pon1000g", "colors"))
ponhit = np.zeros(len(cls), bool)
for k, v in pon.items():
    ponhit[k] = v
print({c: int((cls == c).sum()) for c in ("pos", "partial", "germ", "non", "offnon")})


def metrics(p, posmask, negmask):
    s_pos, s_neg = np.sort(p[posmask]), np.sort(p[negmask])
    ts = np.unique(np.concatenate([s_pos, s_neg]))[::-1]
    tp = len(s_pos) - np.searchsorted(s_pos, ts, side="left")
    fp = len(s_neg) - np.searchsorted(s_neg, ts, side="left")
    rec, prec = tp / len(s_pos), tp / np.maximum(1, tp + fp)
    f1 = 2 * prec * rec / np.maximum(1e-12, prec + rec)
    i = int(np.argmax(f1))
    out = [f"pos {len(s_pos)} neg {len(s_neg)}", f"maxF1 {f1[i]:.3f} (P {prec[i]:.3f} R {rec[i]:.3f})"]
    out += [f"P@R{x}={prec[rec >= x].max():.3f}" for x in (0.5, 0.8, 0.9, 0.935)]
    out += [f"R@P{x}={rec[prec >= x].max():.3f}" if (prec >= x).any() else f"R@P{x}=na" for x in (0.0986, 0.2)]
    return "  ".join(out)


for run in RUN_NAMES:
    p = np.empty(len(cls), np.float32)
    with gzip.open(f"{RUNS}/{run}/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz", "rt") as f:
        for k, line in enumerate(f):
            p[k] = json.loads(line)["p_somatic"]
    pos = cls == "pos"
    print(f"\n== {run}")
    print("  neg = non+offnon+germ+partial      ", metrics(p, pos, np.isin(cls, ["non", "offnon", "germ", "partial"])))
    print("  neg = non+offnon+germ (partial out)", metrics(p, pos, np.isin(cls, ["non", "offnon", "germ"])))
    print("  neg = non+offnon (germline out)    ", metrics(p, pos, np.isin(cls, ["non", "offnon"])))
    print("  neg = non on-ref only              ", metrics(p, pos, cls == "non"))
    keep = ~ponhit
    print("  PoN-filtered, neg = non+offnon+germ+partial", metrics(p, pos & keep, np.isin(cls, ["non", "offnon", "germ", "partial"]) & keep))

"""chr1 v2 SNV test tensors against the ClairS-TO-style PoNs (the repo's scripts/filter_panel_of_normals.py rules),
then the truth-level PR with PoN-matching calls removed. Read-only."""
import collections, csv, gzip, json, sys
import numpy as np
import pysam
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning")
from pansoma_net_v2.bed import Bed, in_region

D = "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors"
S = D + "/SNV"
RUNS = "/scratch/jshen/data/pansoma_net_v2_runs"
RUN_NAMES = ["HG008_Illumina_SNV_base", "HG008_Illumina_SNV_scalars", "HG008_Illumina_SNV_keephard_w100",
             "HG008_Illumina_SNV_allnon_w100"]
PON = "/scratch/jshen/data/Pansoma/panel_of_normal_VCFs"
OUT = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/v2_pipeline"
man = json.load(open(S + "/labels.manifest.json"))
conf = Bed.read(man["truth"]["somatic"]["bed"]).intersect(Bed.read(man["truth"]["germline"]["bed"]))
tt = set()
with open(D + "/somatic.recall.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["chrom"] == "chr1" and r["kind"] == "SNP" and r["passed"] == "True" and r["in_bed"] == "True":
            tt.add(int(r["truth_id"]))

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
        labs.append((lab["candidate_id"], e, ids, (g["pos0"], g["ref"], g["alt"]) if g else None))
n = len(labs)
P = {}
for run in RUN_NAMES:
    p = np.empty(n, np.float32)
    with gzip.open(f"{RUNS}/{run}/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz", "rt") as f:
        for k, line in enumerate(f):
            p[k] = json.loads(line)["p_somatic"]
    P[run] = p
pmax = np.max(np.stack(list(P.values())), axis=0)

# which tensors to annotate: on-reference test tensors labelled 1 or 2, and non with any run's p_somatic >= 0.02
todo = [k for k, (cid, e, ids, g) in enumerate(labs) if e >= 0 and g is not None and (e in (1, 2) or pmax[k] >= 0.02)]
print("tensors to annotate", len(todo), collections.Counter(labs[k][1] for k in todo))
fa = pysam.FastaFile(man["fasta"])
vcfs = {name: pysam.VariantFile(f"{PON}/{f}") for name, f in (
    ("gnomad", "af-only-gnomad.hg38.vcf.gz"), ("dbsnp", "Homo_sapiens_assembly38.dbsnp138.vcf.gz"),
    ("pon1000g", "1000g_pon.hg38.vcf.gz"), ("colors", "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"))}
ann = {}
refbad = 0
for k in todo:
    pos0, ref, alt = labs[k][3]
    if fa.fetch("chr1", pos0, pos0 + 1).upper() != ref:
        refbad += 1
    a = dict(gnomad_af=0.0, gnomad=False, dbsnp=False, pon1000g=False, colors=False)
    for r in vcfs["gnomad"].fetch("chr1", pos0, pos0 + 1):
        if r.pos == pos0 + 1 and r.ref == ref and alt in (r.alts or ()):
            af = r.info["AF"][r.alts.index(alt)]
            a["gnomad_af"] = max(a["gnomad_af"], af)
    a["gnomad"] = a["gnomad_af"] >= 0.001
    for r in vcfs["dbsnp"].fetch("chr1", pos0, pos0 + 1):
        if r.pos == pos0 + 1 and r.ref == ref and alt in (r.alts or ()):
            a["dbsnp"] = True
    for r in vcfs["pon1000g"].fetch("chr1", pos0, pos0 + 1):
        if r.pos == pos0 + 1:
            a["pon1000g"] = True
    for r in vcfs["colors"].fetch("chr1", pos0, pos0 + 1):
        if r.pos == pos0 + 1 and max(r.info["AF"]) >= 0.001:
            a["colors"] = True
    ann[k] = a
print("grch38 REF != FASTA:", refbad)

# PoN membership by class
for lab_value, name in ((1, "somatic(eval 1)"), (2, "germline(eval 2)"), (0, "non, p_max>=0.02")):
    ks = [k for k in todo if labs[k][1] == lab_value]
    c = collections.Counter()
    for k in ks:
        a = ann[k]
        for f in ("gnomad", "dbsnp", "pon1000g", "colors"):
            c[f] += a[f]
        c["any_pon"] += a["gnomad"] or a["dbsnp"] or a["pon1000g"] or a["colors"]
        c["gnomad_any_af"] += a["gnomad_af"] > 0
        c["gnomad_or_1000g_or_colors (no dbSNP)"] += a["gnomad"] or a["pon1000g"] or a["colors"]
    print(f"{name}: n={len(ks)} " + " ".join(f"{f}={v} ({v/len(ks):.1%})" for f, v in c.items()))
# the 697 truths: how many would a PoN drop
lost = collections.Counter()
for k in todo:
    if labs[k][1] == 1 and labs[k][2] & tt:
        a = ann[k]
        lost["any"] += a["gnomad"] or a["dbsnp"] or a["pon1000g"] or a["colors"]
        lost["no_dbsnp"] += a["gnomad"] or a["pon1000g"] or a["colors"]
        lost["dbsnp"] += a["dbsnp"]
print("697-truth tensors hit by PoN:", dict(lost))


def pr(run, drop):
    p = P[run]
    best = {}
    neg = collections.defaultdict(dict)
    for k, (cid, e, ids, g) in enumerate(labs):
        if e < 0 or (k in ann and drop(ann[k])):
            continue
        if e == 1 and ids & tt:
            for t in ids & tt:
                best[t] = max(best.get(t, -1), p[k])
        else:
            kind = "other_kind_label1" if e == 1 else ("germline" if e == 2 else ("non_offref" if g is None else "non_onref"))
            key = g if g else cid
            neg[kind][key] = max(neg[kind].get(key, -1), p[k])
    pos = np.sort(np.array(list(best.values())))
    negs = np.sort(np.concatenate([np.array(list(v.values())) for v in neg.values()]))
    ts = np.unique(np.concatenate([pos, negs]))[::-1]
    tps = len(pos) - np.searchsorted(pos, ts, side="left")
    fps = len(negs) - np.searchsorted(negs, ts, side="left")
    rec, prec = tps / 697, tps / np.maximum(1, tps + fps)
    f1 = 2 * prec * rec / np.maximum(1e-12, prec + rec)
    i = int(np.argmax(f1))
    out = [f"ceiling {len(pos)/697:.3f}", f"maxF1 {f1[i]:.3f} (P {prec[i]:.3f} R {rec[i]:.3f} TP {tps[i]} FP {fps[i]})"]
    out += [f"R@P{x:.2f}={rec[prec >= x].max():.3f}" if (prec >= x).any() else f"R@P{x:.2f}=na" for x in (0.05, 0.1, 0.15, 0.2, 0.3, 0.5)]
    out += [f"P@R{x:.1f}={prec[rec >= x].max():.3f}" if (rec >= x).any() else f"P@R{x:.1f}=na" for x in (0.5, 0.7, 0.8, 0.9)]
    # FP composition at R 0.8
    j = np.flatnonzero(rec >= 0.8)
    if len(j):
        t = ts[j[0]]
        out.append("FP@R0.8 by kind " + str({kk: int((np.array(list(v.values())) >= t).sum()) for kk, v in neg.items()}))
    return "  ".join(out)


filters = {
    "no PoN": lambda a: False,
    "ClairS-TO PoN (gnomAD AF>=0.001 allele, dbSNP138 allele, 1000G PoN pos, CoLoRSdb AF>=0.001 pos)":
        lambda a: a["gnomad"] or a["dbsnp"] or a["pon1000g"] or a["colors"],
    "PoN without dbSNP": lambda a: a["gnomad"] or a["pon1000g"] or a["colors"],
}
for run in RUN_NAMES:
    print(f"\n== {run}")
    for name, fn in filters.items():
        print(f"  [{name}]\n    {pr(run, fn)}")
with open(OUT + "/chr1_pon_annotations.tsv", "w") as f:
    f.write("index\tcandidate_id\teval_label\tpos0\tref\talt\tgnomad_af\tgnomad\tdbsnp\tpon1000g\tcolors\t" + "\t".join(RUN_NAMES) + "\n")
    for k in todo:
        a = ann[k]
        f.write(f"{k}\t{labs[k][0]}\t{labs[k][1]}\t{labs[k][3][0]}\t{labs[k][3][1]}\t{labs[k][3][2]}\t{a['gnomad_af']:.3g}\t"
                f"{int(a['gnomad'])}\t{int(a['dbsnp'])}\t{int(a['pon1000g'])}\t{int(a['colors'])}\t"
                + "\t".join(f"{P[r][k]:.5f}" for r in RUN_NAMES) + "\n")
print("done")

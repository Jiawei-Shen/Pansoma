"""Read-only chr1 SNV statistics of the v2 (v3_tensors) HG008 Illumina set."""
import collections, gzip, json, re, sys, csv
import numpy as np
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning")
from pansoma_net_v2.bed import Bed, in_region

D = "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors"
S = D + "/SNV"
OUT = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/v2_pipeline"
man = json.load(open(S + "/labels.manifest.json"))
sbed = Bed.read(man["truth"]["somatic"]["bed"])
gbed = Bed.read(man["truth"]["germline"]["bed"])
conf = sbed.intersect(gbed)

NUM = {k: re.compile(rb'"' + k.encode() + rb'": ([-+0-9.eE]+)') for k in
       ("af", "coverage", "alt_count", "allele_count", "site_coverage", "second_allele_af")}
CAND = re.compile(rb'"candidate_id": "([^"]+)"')
DS = re.compile(rb'"downsampled_from"')

rows = []  # per tensor
with open(S + "/chr1_variant_summary.ndjson", "rb") as fs, open(S + "/chr1_labels.ndjson") as fl:
    for line in fs:
        lab = json.loads(fl.readline())
        cid = CAND.search(line).group(1).decode()
        assert cid == lab["candidate_id"]
        v = {k: float(p.search(line).group(1)) for k, p in NUM.items()}
        g = lab["grch38"]
        inreg = in_region(conf, "SNP", g, lab.get("anchor"))
        som = lab.get("somatic") or []
        rep_ids = sorted({d["truth_id"] for d in som if d.get("representative")})
        part = lab.get("partial_truth")
        rows.append(dict(cid=cid, label=lab["label"], reason=lab["reason"], af=v["af"], cov=v["coverage"],
                         alt=v["alt_count"], nall=v["allele_count"], a2=v["second_allele_af"],
                         inreg=inreg, off=g is None, rep_ids=rep_ids,
                         part_id=part["truth_id"] if part else None, partial=lab.get("partial"),
                         down=bool(DS.search(line)),
                         gpos=(g["chrom"], g["pos0"], g["ref"], g["alt"]) if g else None))
print("chr1 SNV tensors", len(rows))

# 1. labels x reasons
c = collections.Counter((r["label"], r["reason"]) for r in rows)
print("\n# label, reason, count, in_confident_region")
cin = collections.Counter((r["label"], r["reason"]) for r in rows if r["inreg"])
for (l, rs), n in sorted(c.items()):
    print(f"{l:>3} {rs:55s} {n:>8} {cin[(l, rs)]:>8}")
print("labels:", collections.Counter(r["label"] for r in rows))
print("in-region labels:", collections.Counter(r["label"] for r in rows if r["inreg"]))

# eval set as in pansoma_net_v2.data: label, off_reference_no_truth_match -> 0, not in region -> -1
ev = []
for r in rows:
    e = 0 if r["reason"] == "off_reference_no_truth_match" else r["label"]
    if not r["inreg"]:
        e = -1
    r["eval"] = e
    ev.append(e)
print("eval labels:", collections.Counter(ev))

# 2. AF distribution by class
def q(vals):
    a = np.array(vals)
    if len(a) == 0:
        return "n=0"
    qs = np.quantile(a, [0.05, 0.25, 0.5, 0.75, 0.95])
    return f"n={len(a)} mean={a.mean():.3f} q05/25/50/75/95=" + "/".join(f"{x:.3f}" for x in qs)
bins = [0, 0.06, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0001]
def hist(vals):
    h, _ = np.histogram(np.array(vals), bins=bins)
    return " ".join(f"{int(x)}" for x in h)
print("\n# AF by eval class (test set of v2), bins", bins)
for name, sel in (("somatic(eval=1)", lambda r: r["eval"] == 1), ("germline(eval=2)", lambda r: r["eval"] == 2),
                  ("non(eval=0)", lambda r: r["eval"] == 0),
                  ("non(eval=0, on-ref)", lambda r: r["eval"] == 0 and not r["off"]),
                  ("non off-ref (eval 0)", lambda r: r["eval"] == 0 and r["off"]),
                  ("somatic rep exact (eval=1, reason rep)", lambda r: r["eval"] == 1 and r["reason"] == "representative_allele_is_somatic_truth"),
                  ("somatic partial (eval=1)", lambda r: r["eval"] == 1 and r["reason"] != "representative_allele_is_somatic_truth"),
                  ("below_snv_min_af (-1)", lambda r: r["reason"] == "below_snv_min_af"),
                  ("all label 1", lambda r: r["label"] == 1), ("all label 2", lambda r: r["label"] == 2),
                  ("all label 0", lambda r: r["label"] == 0)):
    vals = [r["af"] for r in rows if sel(r)]
    print(f"{name:42s} {q(vals)}\n{'':42s} hist {hist(vals)}")
    cov = [r["cov"] for r in rows if sel(r)]
    print(f"{'':42s} coverage {q(cov)}")

# 3. the 697 truths
truth = {}
with open(D + "/somatic.recall.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["chrom"] == "chr1" and r["kind"] == "SNP":
            truth[int(r["truth_id"])] = r
tt = {t: r for t, r in truth.items() if r["passed"] == "True" and r["in_bed"] == "True"}
print("\nchr1 SNP truth alleles", len(truth), "PASS+in somatic BED", len(tt))
print("status:", collections.Counter((r["status"], r["detail"]) for r in tt.values()))
inconf = {t: conf.contains("chr1", int(r["vcf_pos"]) - 1, int(r["vcf_pos"])) for t, r in tt.items()}
ingb = {t: gbed.contains("chr1", int(r["vcf_pos"]) - 1, int(r["vcf_pos"])) for t, r in tt.items()}
print("of the 697 in germline BED:", sum(ingb.values()), " in somatic∩germline:", sum(inconf.values()))
print("status x in-confident:", collections.Counter((r["status"], inconf[t]) for t, r in tt.items()))
# tensors per truth
by_truth = collections.defaultdict(list)
for i, r in enumerate(rows):
    for t in r["rep_ids"]:
        by_truth[t].append((i, "rep"))
    if r["part_id"] is not None and r["label"] == 1:
        by_truth[r["part_id"]].append((i, "partial"))
has_rep = {t for t in tt if any(k == "rep" and rows[i]["label"] == 1 for i, k in by_truth[t])}
has_any = {t for t in tt if any(rows[i]["label"] == 1 for i, k in by_truth[t])}
has_eval = {t for t in tt if any(rows[i]["eval"] == 1 for i, k in by_truth[t])}
print("697: with a label-1 SNV tensor as representative:", len(has_rep), " any label-1 SNV tensor (rep or partial):",
      len(has_any), " a label-1 tensor in the test region (eval 1):", len(has_eval))
# labels of representative-matching tensors (any label)
rep_labels = collections.Counter()
for t in tt:
    for i, k in by_truth[t]:
        rep_labels[(k, rows[i]["label"], rows[i]["reason"], rows[i]["inreg"])] += 1
print("tensor-matches of 697 (kind, label, reason, inreg):", rep_labels)
# truth-matching tensors per truth
mult = collections.Counter(len(by_truth[t]) for t in tt)
print("number of label-1 tensors per truth (of 697):", sorted(mult.items()))
af_truth = [max(rows[i]["af"] for i, k in by_truth[t] if k == "rep") for t in has_rep]
print("AF of the representative tensor of each found truth:", q(af_truth), "hist", hist(af_truth))
cov_truth = [max(rows[i]["cov"] for i, k in by_truth[t] if k == "rep") for t in has_rep]
print("coverage of those tensors:", q(cov_truth))

# label-1 tensors whose truth is not SNP (other kind) or outside the 697
ids_all = set()
oth = collections.Counter()
for r in rows:
    if r["label"] == 1:
        for t in r["rep_ids"] + ([r["part_id"]] if r["part_id"] is not None else []):
            oth["snp697" if t in tt else ("snp_other" if t in truth else "non_snp_truth")] += 1
print("label-1 tensor->truth links by truth set:", oth)

# 4. v1-style binary true sets
def summarize(name, sel):
    s = [r for r in rows if sel(r)]
    ids = set()
    for r in s:
        ids.update(r["rep_ids"])
        if r["part_id"] is not None:
            ids.add(r["part_id"])
    print(f"{name:75s} tensors={len(s):>7} distinct_truth={len(ids):>5} of697={len(ids & set(tt)):>4}")
print("\n# v1-style 'true' sets on chr1")
summarize("label==1 (all)", lambda r: r["label"] == 1)
summarize("label==1 & in confident region", lambda r: r["label"] == 1 and r["inreg"])
summarize("label==1 & representative exact SNP match", lambda r: r["label"] == 1 and r["reason"] == "representative_allele_is_somatic_truth")
summarize("label==1 & rep exact & in region", lambda r: r["label"] == 1 and r["reason"] == "representative_allele_is_somatic_truth" and r["inreg"])
summarize("label==1 & rep exact & in region & AF>=0.07", lambda r: r["label"] == 1 and r["reason"] == "representative_allele_is_somatic_truth" and r["inreg"] and r["af"] >= 0.07)
summarize("label==1 & rep exact & AF>=0.1", lambda r: r["label"] == 1 and r["reason"] == "representative_allele_is_somatic_truth" and r["af"] >= 0.1)
summarize("label==1 & partial (not rep)", lambda r: r["label"] == 1 and r["reason"] != "representative_allele_is_somatic_truth")
summarize("somatic truth among site alleles (any allele, any label)", lambda r: bool(r["rep_ids"]) or r["reason"] == "truth_matches_non_representative_allele")
print("rows with downsampled_from:", sum(r["down"] for r in rows))
# multi-allele
print("allele_count distribution:", collections.Counter(int(r["nall"]) for r in rows))
print("off-reference tensors:", sum(r["off"] for r in rows), " off-ref in region:", sum(r["off"] and r["inreg"] for r in rows))

# save compact arrays for later
np.savez_compressed(OUT + "/chr1_snv_rows.npz", af=np.array([r["af"] for r in rows], np.float32),
                    cov=np.array([r["cov"] for r in rows], np.int32), label=np.array([r["label"] for r in rows], np.int8),
                    eval=np.array(ev, np.int8), inreg=np.array([r["inreg"] for r in rows]),
                    off=np.array([r["off"] for r in rows]))
with open(OUT + "/chr1_697_truth.tsv", "w") as f:
    f.write("truth_id\tvcf_pos\tstatus\tdetail\tin_germline_bed\tin_confident\tn_label1_tensors\trep_tensor_af\trep_tensor_cov\trep_tensor_label\trep_tensor_inreg\n")
    for t, r in sorted(tt.items(), key=lambda x: int(x[1]["vcf_pos"])):
        reps = [i for i, k in by_truth[t] if k == "rep"]
        best = max(reps, key=lambda i: rows[i]["af"]) if reps else None
        f.write(f"{t}\t{r['vcf_pos']}\t{r['status']}\t{r['detail']}\t{ingb[t]}\t{inconf[t]}\t{len(by_truth[t])}\t"
                + (f"{rows[best]['af']:.4f}\t{int(rows[best]['cov'])}\t{rows[best]['label']}\t{rows[best]['inreg']}" if best is not None else "\t\t\t")
                + "\n")
print("done")

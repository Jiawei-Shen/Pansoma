"""v1 chr1 SNV candidates (val set with genomic_position + label) vs the 697 truths and the germline truth (read-only)."""
import json
import sys
from collections import Counter, defaultdict

sys.path.insert(0, "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head")
from common import truth_697, confident_chr1, somatic_bed_chr1, COMP  # noqa

VAL = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset/variant_summary_classified.ndjson"
GERM = "/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/truth/germline.graph.tsv"
SCR = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"

truth = truth_697()
conf = confident_chr1()
sbed = somatic_bed_chr1()
tpos = defaultdict(list)
for tid, t in truth.items():
    tpos[t["pos0"]].append(tid)

# germline chr1 SNP truth (PASS) positions
germ = {}
with open(GERM) as f:
    hdr = next(f).rstrip("\n").split("\t")
    ic, ik, ip, ia, ipass = hdr.index("chrom"), hdr.index("kind"), hdr.index("pos0"), hdr.index("alt"), hdr.index("passed")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[ic] != "chr1":
            continue
        if p[ic] > "chr1" and not p[ic].startswith("chr1"):
            break
        if p[ik] == "SNP":
            germ[int(p[ip])] = p[ia]
print("germline chr1 SNP truth", len(germ))

found_same = defaultdict(list)   # tid -> list of (label, af, strand_match)
cnt = Counter()
recs = []
with open(VAL) as f:
    for line in f:
        r = json.loads(line)
        pos0 = r["genomic_position"] - 1
        alt = r["v_alt"]
        ref = r["v_ref"]
        lab = r["label_int"]
        inconf = conf.contains0(pos0)
        g = germ.get(pos0)
        isgerm = g is not None and (g == alt or g == alt.translate(COMP))
        cnt[("label", lab, "conf", inconf, "germ", isgerm)] += 1
        for tid in tpos.get(pos0, ()):
            t = truth[tid]
            m = "same" if t["alt"] == alt else ("comp" if t["alt"] == alt.translate(COMP) else "other_alt")
            found_same[tid].append((lab, r["alt_allele_frequency"], m, r["coverage_at_locus"]))
        recs.append((pos0, alt, lab, inconf, isgerm, r["shard_index"], r["index_within_shard"], r["alt_allele_frequency"]))
print(sorted(cnt.items(), key=lambda x: -x[1]))
any_c = sum(1 for tid in truth if found_same.get(tid))
exact = sum(1 for tid in truth if any(m in ("same", "comp") for _, _, m, _ in found_same.get(tid, [])))
true_lab = sum(1 for tid in truth if any(l == 1 for l, _, m, _ in found_same.get(tid, []) if m in ("same", "comp")))
mism = Counter(m for tid in truth for _, _, m, _ in found_same.get(tid, []))
print("v1 val chr1: truths with any candidate at pos", any_c, "with matching allele", exact, "labeled true", true_lab,
      "allele-match kinds", dict(mism))
# v1 true labels not matching the 697
v1_true = [(p, a) for p, a, l, *_ in recs if l == 1]
in697 = sum(1 for p, a in v1_true if any(truth[t]["alt"] in (a, a.translate(COMP)) for t in tpos.get(p, ())))
print("v1 true-labeled tensors", len(v1_true), "matching a 697 truth", in697,
      "in somatic BED", sum(sbed.contains0(p) for p, a in v1_true), "in confident", sum(conf.contains0(p) for p, a in v1_true))
# v2 status of truths missing in v1
miss = [tid for tid in truth if not any(m in ("same", "comp") for _, _, m, _ in found_same.get(tid, []))]
print("v2 status of truths v1 lacks:", Counter(truth[t]["status"] for t in miss))
has_v1 = set(truth) - set(miss)
print("v2 status of truths v1 has:", Counter(truth[t]["status"] for t in has_v1))
json.dump(dict(v1_any=any_c, v1_allele=exact, v1_true_label=true_lab, n=len(truth),
               missing=[dict(tid=t, **truth[t]) for t in miss]), open(f"{SCR}/v1_ceiling.json", "w"), indent=1)
import pickle
pickle.dump(recs, open(f"{SCR}/v1_val_recs.pkl", "wb"))

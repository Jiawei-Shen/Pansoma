"""Score full v1 chr1 predictions (one float per row of 5ch_testing_data_SNV_chr1/variant_summary.ndjson) against the
697 truths with the h2h_final.py definitions (region = somatic BED ∩ germline BED; FP unit = distinct (pos0, alt);
v2 label-1 non-697 sites neutral; off-reference rows = rows without a GRCh38 position, one FP unit per row).
usage: score_v1_full.py <p.npy> [<p.npy> ...]   |   score_v1_full.py --selftest  (rebuilds p from the e053 VCFs)"""
import gzip, json, sys
import numpy as np
HERE = "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/head_to_head"
sys.path.insert(0, HERE)
from common import V2, truth_697, confident_chr1  # noqa
from curve import curve, summary  # noqa
from pon_annot import annotate, rule_repo, rule_af001  # noqa
TEST = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_testing_data_SNV/5ch_testing_data_SNV_chr1/variant_summary.ndjson"
VAL = "/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset/variant_summary_classified.ndjson"
truth = truth_697(); NT = len(truth)
tkey = {(t["pos0"], t["alt"]): tid for tid, t in truth.items()}
conf = confident_chr1()
neutral = set()
for line in open(f"{V2}/SNV/chr1_labels.ndjson"):
    r = json.loads(line)
    if r["label"] == 1 and r.get("grch38") and (r["grch38"]["pos0"], r["grch38"]["alt"]) not in tkey:
        neutral.add((r["grch38"]["pos0"], r["grch38"]["alt"]))
val = {}
for line in open(VAL):
    r = json.loads(line)
    val[(r["node_id"], r["variant_key"])] = (r["genomic_position"] - 1, r["v_ref"].upper(), r["v_alt"].upper())
site = []
for line in open(TEST):
    r = json.loads(line)
    site.append(val.get((r["node_id"], r["variant_key"])))
    assert len(site) - 1 == r["shard_index"] * 32768 + r["index_within_shard"]
print("test rows", len(site), "GRCh38", sum(s is not None for s in site), "off-reference", sum(s is None for s in site))


def score(p, name):
    assert len(p) == len(site)
    idx = np.flatnonzero(p >= 0.005)
    ann = annotate({site[i] for i in idx if site[i] is not None})
    out = {}
    for rn, drop in (("noPoN", lambda k: False), ("PoN_repo", lambda k: rule_repo(ann[k])), ("PoN_af001", lambda k: rule_af001(ann[k]))):
        tp, fp, off = {}, {}, []
        for i, s in enumerate(site):
            if s is None:
                off.append(float(p[i])); continue
            k2 = (s[0], s[2])
            if k2 in tkey:
                if p[i] >= 0.005 and drop(s):
                    continue
                tp[tkey[k2]] = max(tp.get(tkey[k2], 0.0), float(p[i]))
            elif conf.contains0(s[0]) and k2 not in neutral:
                if p[i] >= 0.005 and drop(s):
                    continue
                fp[k2] = max(fp.get(k2, 0.0), float(p[i]))
        out[f"{rn} | {name} | GRCh38 FP"] = summary(curve(list(tp.values()), list(fp.values()), NT), NT, len(tp))
        out[f"{rn} | {name} | + off-ref FP"] = summary(curve(list(tp.values()), list(fp.values()) + off, NT), NT, len(tp))
        out[f"{rn} | {name} | GRCh38 FP"]["at_0.5"] = dict(TP=sum(v >= 0.5 for v in tp.values()), FP=sum(v >= 0.5 for v in fp.values()), offref=sum(v >= 0.5 for v in off))
    for k, s in out.items():
        m = s["maxF1"]
        f = lambda v: "  n/a" if v is None else f"{v[0]:.3f}"
        print(f"{k} | {s['ceiling']:.3f} | " + " ".join(f(s[f'R@P{x:.2f}']) for x in (0.05, 0.10, 0.15, 0.20)) + " | " +
              " ".join(f(s[f'P@R{x:.1f}']) for x in (0.5, 0.7, 0.8, 0.9)) + f" | {m[0]:.3f} ({m[1]:.3f}, {m[2]:.3f})" +
              (f" | at0.5 {s['at_0.5']}" if "at_0.5" in s else ""))
    return out


if sys.argv[1] == "--selftest":
    p = np.zeros(len(site), np.float32)
    R = "/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.{}.vcf.gz"
    for kind in ("linear", "graph"):
        for line in gzip.open(R.format(kind), "rt"):
            if line.startswith("#"):
                continue
            info = dict(kv.split("=", 1) for kv in line.split("\t")[7].split(";"))
            i = int(info["SHARD"]) * 32768 + int(info["IDX"])
            p[i] = max(p[i], float(info["PROB"]))
    res = score(p, "selftest e053 from VCF")
else:
    res = {}
    for path in sys.argv[1:]:
        res.update(score(np.load(path), path.split("/")[-1]))
json.dump(res, open(f"{HERE}/score_v1_full.{'selftest' if sys.argv[1] == '--selftest' else 'full'}.json", "w"), indent=1)

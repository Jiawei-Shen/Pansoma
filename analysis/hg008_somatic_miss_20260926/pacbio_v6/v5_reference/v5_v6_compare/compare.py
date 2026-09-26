"""Stream v5/v6 autosome summaries + labels; write per-(version, kind) stats and somatic truth tables."""
import glob, json, sys, collections, statistics
from pathlib import Path
P = Path("/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125")
OUT = Path("/scratch/jshen/Github/Pansoma/tmp/v5_v6_compare")
version, kind = sys.argv[1], sys.argv[2]
d = P / f"{version}_tensors" / kind
ev, ac, afs, cov, lens = collections.Counter(), collections.Counter(), [], [], collections.Counter()
lab, reason = collections.Counter(), collections.Counter()
lab_by_type = collections.Counter()
som_rep, som_any, germ_rep = {}, {}, set()
def lbin(n): return "1" if n == 1 else "2-5" if n <= 5 else "6-20" if n <= 20 else "21-50"
for summ in sorted(glob.glob(str(d / "chr*_variant_summary.ndjson"))):
    labf = summ.replace("_variant_summary.ndjson", "_labels.ndjson")
    with open(summ) as fs, open(labf) as fl:
        for ls, ll in zip(fs, fl):
            m, l = json.loads(ls), json.loads(ll)
            assert m["candidate_id"] == l["candidate_id"], (summ, m["candidate_id"], l["candidate_id"])
            t = m["event_type"]; ev[t] += 1
            ac[min(m.get("allele_count", 1), 4)] += 1
            afs.append(m["af"]); cov.append(m["coverage"])
            if kind == "INDEL": lens[(t, lbin(m["event_length"]))] += 1
            lab[l["label_name"]] += 1; reason[l["reason"]] += 1; lab_by_type[(t, l["label_name"])] += 1
            for s in l.get("somatic", []):
                key = s["truth_id"]
                rec = dict(site=l["site_id"], cand=l["candidate_id"], matched=s.get("matched_candidate_id"),
                           rep=s.get("representative"), label=l["label_name"], reason=l["reason"], chrom=l["chrom"],
                           vcf=(s.get("vcf_pos"), s.get("vcf_ref"), s.get("vcf_alt")), af=m["af"], cov=m["coverage"])
                som_any.setdefault(key, rec)
                if l["label_name"] == "somatic": som_rep[key] = rec
            for g in l.get("germline", []):
                if l["label_name"] == "germline": germ_rep.add(g["truth_id"])
def q(x): x = sorted(x); return {p: round(x[int(p/100*(len(x)-1))], 3) for p in (10, 25, 50, 75, 90)} if x else {}
res = dict(n=sum(ev.values()), event_types=dict(ev), allele_count=dict(ac), af_q=q(afs), coverage_q=q(cov),
           lengths={f"{a}:{b}": n for (a, b), n in sorted(lens.items())}, labels=dict(lab), reasons=dict(reason),
           labels_by_type={f"{a}:{b}": n for (a, b), n in sorted(lab_by_type.items())},
           somatic_truth_labelled=len(som_rep), somatic_truth_present_any=len(som_any), germline_truth_labelled=len(germ_rep))
(OUT / f"{version}_{kind}.json").write_text(json.dumps(res, indent=1))
(OUT / f"{version}_{kind}_somatic.json").write_text(json.dumps(dict(rep=som_rep, any=som_any)))
print(version, kind, "done", res["n"])

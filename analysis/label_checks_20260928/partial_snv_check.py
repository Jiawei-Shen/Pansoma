"""SNV tensors labelled 1 through a partial somatic match (truth-labels-v6): what are they? vs exact somatic SNVs."""
import csv, glob, json, collections, random, sys, numpy as np
csv.field_size_limit(sys.maxsize)
R = "/scratch/jshen/data/pansoma_v2_tensors"
SETS = {"HG008 Illumina": "Liss_lab_BCM_Illumina-WGS_20240313", "HG008 PacBio": "Liss_lab_PacBio_Revio_20240125",
        "HG008 ONT": "Liss_lab_Northeastern-ONT-UL-20241216", "COLO829T Illumina": "COLO829T_Illumina"}
rng = random.Random(11); examples = []
for name, d in SETS.items():
    T = f"{R}/{d}/v3_tensors"
    status = {r["truth_id"]: r["status"] for r in csv.DictReader(open(f"{T}/somatic.recall.tsv"), delimiter="\t")}
    c = collections.Counter(); shards = {}; rows = []
    for f in sorted(glob.glob(f"{T}/SNV/chr*_labels.ndjson")):
        shards.clear()
        for a, b in zip(open(f), open(f.replace("_labels.ndjson", "_variant_summary.ndjson"))):
            if '"label": 1,' not in a: continue
            l = json.loads(a); s = json.loads(b)
            g = {x["allele"]: (x["start_row"], x["end_row"]) for x in s["row_groups"]}.get("A1")
            if s["shard_file"] not in shards: shards[s["shard_file"]] = np.load(f"{T}/SNV/{s['shard_file']}", mmap_mode="r")
            strand = np.asarray(shards[s["shard_file"]][s["index_within_shard"], 7, g[0]:g[1], :]).max(axis=1) if g else np.array([])
            fwd, rev = int((strand == 1).sum()), int((strand == 2).sum())
            one = fwd == 0 or rev == 0
            if l["partial"]:
                t = l["partial_truth"]; dl = len(t["vcf_alt"]) - len(t["vcf_ref"])
                tk = "SNV" if dl == 0 and len(t["vcf_ref"]) == 1 else ("INS" if dl > 0 else "DEL" if dl < 0 else "MNV/complex")
                ov = l["overlap"]; ob = "0.45-0.6 (new in v6)" if ov <= 0.6 else ">0.6"
                pos = l["grch38"]["pos0"] + 1 if l["grch38"] else None
                dist = abs(pos - t["vcf_pos"]) if pos else None
                own = status.get(str(t["truth_id"])) == "tensor_representative"
                rows.append(dict(set=name, cid=l["candidate_id"], chrom=l["chrom"], pos=pos, anchor=l.get("anchor"), truth=f"{t['chrom']}:{t['vcf_pos']} {t['vcf_ref'][:12]}>{t['vcf_alt'][:12]}",
                                 tk=tk, ov=ov, ob=ob, af=s["af"], alt=s["alt_count"], cov=s["coverage"], fwd=fwd, rev=rev, one=one, dist=dist, own=own, shard=s["shard_file"], idx=s["index_within_shard"]))
                c[("partial", tk, ob, "one strand" if one else "both strands")] += 1
            else:
                c[("exact", "SNV", "", "one strand" if one else "both strands")] += 1
    n_exact = sum(v for k, v in c.items() if k[0] == "exact"); n_part = len(rows)
    print(f"\n=== {name}: exact somatic SNV {n_exact:,} (one strand {100*c[('exact','SNV','','one strand')]/max(1,n_exact):.1f}%), partial {n_part:,}")
    agg = collections.defaultdict(lambda: collections.Counter())
    for r in rows: agg[(r["tk"], r["ob"])]["n"] += 1; agg[(r["tk"], r["ob"])]["one"] += r["one"]; agg[(r["tk"], r["ob"])]["own"] += r["own"]
    for k in sorted(agg):
        v = agg[k]; afs = sorted(r["af"] for r in rows if (r["tk"], r["ob"]) == k); ds = sorted(r["dist"] for r in rows if (r["tk"], r["ob"]) == k and r["dist"] is not None)
        print(f"   truth {k[0]:4s} overlap {k[1]:22s} n={v['n']:5d}  one-strand {100*v['one']/v['n']:5.1f}%  truth has its own exact tensor {100*v['own']/v['n']:5.1f}%"
              f"  AF p50 {afs[len(afs)//2]:.3f}  dist to truth p50 {ds[len(ds)//2] if ds else '-'} bp")
    if name == "HG008 Illumina":
        for grp in (lambda r: r["one"] and r["ob"].startswith("0.45"), lambda r: not r["one"] and r["ob"].startswith("0.45"), lambda r: r["one"] and r["ob"] == ">0.6", lambda r: not r["one"] and r["ob"] == ">0.6"):
            pool = [r for r in rows if grp(r)]
            examples += rng.sample(pool, min(3, len(pool)))
json.dump(examples, open("partial_snv_examples.json", "w"), indent=1)
print("\nexamples:")
for e in examples: print(" ", {k: e[k] for k in ("cid", "chrom", "pos", "truth", "tk", "ov", "af", "alt", "cov", "fwd", "rev", "own", "dist")})

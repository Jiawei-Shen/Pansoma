"""Illumina SNV: would a strand filter (ALT reads on one strand only) or an AF cut hurt somatic tensors?
All label-1 tensors, and a 1-in-40 sample of label-0 tensors (AF, ALT reads on each strand from channel 7)."""
import json, glob, random, collections, numpy as np
D = "/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors/SNV"
rng = random.Random(7); rows = collections.defaultdict(list); shards = {}
for f in sorted(glob.glob(f"{D}/chr*_labels.ndjson")):
    shards.clear()
    for a, b in zip(open(f), open(f.replace("_labels.ndjson", "_variant_summary.ndjson"))):
        if '"label": 1,' in a: lab = 1
        elif '"label": 0,' in a and rng.random() < 1 / 40: lab = 0
        else: continue
        s = json.loads(b); g = {x["allele"]: (x["start_row"], x["end_row"]) for x in s["row_groups"]}.get("A1")
        if not g: continue
        name = s["shard_file"]
        if name not in shards: shards[name] = np.load(f"{D}/{name}", mmap_mode="r")
        strand = np.asarray(shards[name][s["index_within_shard"], 7, g[0]:g[1], :]).max(axis=1)
        fwd, rev = int((strand == 1).sum()), int((strand == 2).sum())
        rows[lab].append((s["af"], fwd, rev, json.loads(a)["reason"]))
def keep(r, rule):
    af, f, v, _ = r; both = f > 0 and v > 0
    return {"none": True, "AF>=0.10": af >= 0.10, "AF>=0.15": af >= 0.15, "both strands": both,
            "both strands or AF>=0.15": both or af >= 0.15, "min(fwd,rev)>=2 or AF>=0.15": min(f, v) >= 2 or af >= 0.15}[rule]
n1, n0 = len(rows[1]), len(rows[0])
print(f"label 1 tensors {n1:,} (all); label 0 sample {n0:,} (1/40 of ~3.99M)")
print(f"{'rule':30s} {'somatic kept':>14s} {'0 kept (est.)':>16s}")
for rule in ("none", "AF>=0.10", "AF>=0.15", "both strands", "both strands or AF>=0.15", "min(fwd,rev)>=2 or AF>=0.15"):
    k1 = sum(keep(r, rule) for r in rows[1]); k0 = sum(keep(r, rule) for r in rows[0])
    print(f"{rule:30s} {k1:6,} ({100*k1/n1:5.1f}%) {k0*40/1e6:8.2f}M ({100*k0/n0:5.1f}%)")
lost = collections.Counter(r[3] for r in rows[1] if not keep(r, "both strands"))
print("somatic tensors with ALT reads on one strand only, by reason:", dict(lost))

# usage: python check_tensors.py SET CHROM [CHROM ...] -- the SNV tensors of <SET>_snvaf05 against <SET> per chromosome:
# every old site is in the new set; a site with the same alleles has the same tensor (bytes; all such sites, or a
# random 30,000 per chromosome) and the same label; the
# sites only in the new set, and old sites that gained an allele, have the new allele(s) at AF 0.05-0.08.
import json, sys
import numpy as np
from collections import Counter
S, chroms = sys.argv[1], sys.argv[2:]
SAMPLE = 30000  # same-allele sites whose tensors are compared per chromosome (all when fewer)
OLD, NEW = (f"/scratch/jshen/data/pansoma_v2_tensors/{x}/tensors/SNV" for x in (S, f"{S}_snvaf05"))

def summary(d, c):
    out = {}
    for l in open(f"{d}/{c}_variant_summary.ndjson"):
        r = json.loads(l)
        out[r["site_id"]] = (r["shard_file"], r["index_within_shard"], tuple(sorted(r["allele_labels"].values())), r["af"])
    return out

def labels(d, c):
    return {r["site_id"]: (r["label"], r.get("reason")) for r in map(json.loads, open(f"{d}/{c}_labels.ndjson"))}

def allele_afs(d, c):
    """candidate_id -> AF for every allele of every site (alleles[] of the summary)."""
    import ast
    out = {}
    for l in open(f"{d}/{c}_variant_summary.ndjson"):
        r = json.loads(l)
        al = r["alleles"]; al = ast.literal_eval(al) if isinstance(al, str) else al
        for a in al: out[a["candidate_id"]] = a.get("af", r["af"])
    return out

ok = True
for c in chroms:
    o, n = summary(OLD, c), summary(NEW, c)
    lo, ln = labels(OLD, c), labels(NEW, c)
    missing = [s for s in o if s not in n]
    same = [s for s in o if s in n and o[s][2] == n[s][2]]
    changed = [s for s in o if s in n and o[s][2] != n[s][2]]
    added = [s for s in n if s not in o]
    heads = {}
    def row(d, f, i):
        """Tensor i of an .npy shard by a seek and a read (no mmap: mapped pages count as the job's memory)."""
        p = f"{d}/{f}"
        if p not in heads:
            with open(p, "rb") as h:
                version = np.lib.format.read_magic(h)
                read = np.lib.format.read_array_header_1_0 if version == (1, 0) else np.lib.format.read_array_header_2_0
                shape, _, dtype = read(h)
                heads[p] = (h.tell(), shape, dtype)
        off, shape, dtype = heads[p]
        n = int(np.prod(shape[1:])) * dtype.itemsize
        with open(p, "rb") as h:
            h.seek(off + i * n); return np.frombuffer(h.read(n), dtype).reshape(shape[1:])
    rng = np.random.default_rng(0)
    checked = same if len(same) <= SAMPLE else [same[i] for i in sorted(rng.choice(len(same), SAMPLE, replace=False))]
    checked = sorted(checked, key=lambda s: (o[s][0], o[s][1]))
    diff = sum(not np.array_equal(row(OLD, *o[s][:2]), row(NEW, *n[s][:2])) for s in checked)
    ldiff = sum(lo.get(s) != ln.get(s) for s in same)
    na = allele_afs(NEW, c); oa = set(a for s in o for a in o[s][2])
    new_alleles = [a for s in added + changed for a in n[s][2] if a not in oa]
    afs = [na[a] for a in new_alleles if a in na]
    out_of_band = sum(not (0.05 <= x < 0.08) for x in afs)
    lab_added = Counter(ln[s][0] for s in added if s in ln)
    print(f"{S} {c}: old sites {len(o):,}, new {len(n):,} (+{len(added):,} = +{len(added)/len(o):.1%}); missing {len(missing)}; "
          f"same alleles {len(same):,}: tensor differs {diff} of {len(checked):,} compared, label differs {ldiff}; gained an allele {len(changed):,}; "
          f"new alleles {len(new_alleles):,}, AF range {min(afs, default=0):.4f}-{max(afs, default=0):.4f}, outside [0.05, 0.08) {out_of_band}; "
          f"labels of added sites {dict(sorted(lab_added.items()))} (-1 ignored, 0 non, 1 somatic, 2 germline)")
    ok &= not missing and diff == 0 and ldiff == 0 and out_of_band == 0
print("OK" if ok else "CHECK FAILED"); sys.exit(0 if ok else 1)

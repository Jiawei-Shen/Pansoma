# usage: python check_config.py SET -- the new run's config.json equals the old one except paths of the new root, the
# SNV AF floor (0.08 -> 0.05), the frozen source and the times; exits 1 on any other difference.
import json, sys
S = sys.argv[1]
old = json.load(open(f"/scratch/jshen/data/pansoma_v2_tensors/{S}/run/config.json"))
new = json.load(open(f"/scratch/jshen/data/pansoma_v2_tensors/{S}_snvaf05/run/config.json"))
skip = {"created", "source_sha256", "python", "package", "tensors", "parts", "schedule", "record_index", "native_decoder"}
bad = []
def cmp(a, b, path):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            if path == "" and k in skip: continue
            cmp(a.get(k), b.get(k), f"{path}.{k}")
    else:
        if isinstance(a, str): a = a.replace(f"/{S}/", f"/{S}_snvaf05/")
        if a != b: bad.append((path, a, b))
cmp(old, new, "")
# must not differ: what decides the tensors (every builder parameter, partition, inputs, chromosomes, merge, INDEL floor)
strict = (".builder", ".tasks", ".processes", ".inputs.gam", ".inputs.index", ".inputs.graph_index",
          ".chromosome_selection.chromosomes", ".chromosome_selection.nodes_in", ".chromosome_selection.nodes_kept",
          ".postprocess.merge_shard_size", ".postprocess.reference_path", ".variant_outputs.INDEL")
fail = [p for p, _, _ in bad if p.startswith(strict) and p != ".variant_outputs.SNV"]
for p, a, b in bad: print(("DIFF " if p in fail else "note ") + f"{p}: {str(a)[:200]} -> {str(b)[:200]}")
print("builder", new["builder"]); print("variant_outputs", new["variant_outputs"])
print("native decoder", new.get("native_decoder"), "record index", (new.get("record_index") or {}).get("path") if isinstance(new.get("record_index"), dict) else new.get("record_index"))
assert new['variant_outputs'] == {'SNV': 0.05, 'INDEL': 0.08}
sys.exit(1 if fail else 0)

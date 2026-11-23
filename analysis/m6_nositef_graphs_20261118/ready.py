#!/usr/bin/env python3
"""Is a tensor set's copy on storage2 complete? Exit 0 when <SET>/tensors/<KIND> has a complete manifest, a label
manifest for the same number of tensors, and every shard of the manifest with its data file at the size the manifest
implies (npy header + tensors x shape) and its labels file; else exit 1 with the first reason on stdout. submit.sh
calls it per set, so a set that is still being copied is not trained on. Standard library only (login node python3).
usage: ready.py TENSORS_DIR KIND"""
import json
import os
import sys

HEADER = 128  # npy header of the shard files


def reason(d):
    try:
        m = json.load(open(os.path.join(d, "manifest.json")))
        lm = json.load(open(os.path.join(d, "labels.manifest.json")))
    except (OSError, ValueError) as e:
        return f"no readable manifest: {getattr(e, 'strerror', None) or e}"
    if m.get("status") != "complete":
        return f"manifest status {m.get('status')}"
    if lm.get("tensors") != m["tensors"]:
        return f"labels for {lm.get('tensors')} tensors, manifest {m['tensors']}"
    per = 1
    for x in m["shape"]:
        per *= x
    for chrom, c in m["chromosomes"].items():
        if c["tensors"] and not os.path.exists(os.path.join(d, f"{chrom}_labels.ndjson")):
            return f"{chrom}_labels.ndjson missing"
        for s in c["shards"]:
            data = os.path.join(d, s["file"])
            labels = data.replace("_data.npy", "_labels.npy")
            for f, size in ((data, HEADER + per * s["tensors"]), (labels, HEADER + s["tensors"])):
                if not os.path.exists(f):
                    return f"{os.path.basename(f)} missing"
                if os.path.getsize(f) != size:
                    return f"{os.path.basename(f)} is {os.path.getsize(f)} bytes, expected {size}"
    return None


if __name__ == "__main__":
    r = reason(os.path.join(sys.argv[1], sys.argv[2]))
    if r:
        print(r)
    sys.exit(1 if r else 0)

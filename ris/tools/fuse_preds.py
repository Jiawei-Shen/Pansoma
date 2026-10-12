"""Seed ensemble of predict outputs (ris experiment): the class probabilities of K models of one recipe, averaged as a
normalized geometric mean (mean log-probability per class, then softmax), no fitted parameters.

    python ris/tools/fuse_preds.py OUT_DIR MEMBER_DIR [MEMBER_DIR ...]

Each MEMBER_DIR is a predict --output of the same set and chromosomes (one *.predictions.ndjson.gz + .metrics.json).
The records must list the same candidates in the same order (checked). The fused file keeps member 0's fields
(labels, reasons, truth ids) with p_non / p_somatic / p_germline replaced, pred = the argmax class, and member 0's
.metrics.json beside it (graph_vcf reads the checkpoint path from it)."""
import glob
import gzip
import json
import math
import os
import shutil
import sys

CLASSES = ("p_non", "p_somatic", "p_germline")
NAMES = {"p_non": "non", "p_somatic": "somatic", "p_germline": "germline"}


def main():
    out, members = sys.argv[1], sys.argv[2:]
    files = []
    for m in members:
        f = glob.glob(f"{m}/*.predictions.ndjson.gz")
        if len(f) != 1:
            raise SystemExit(f"{m}: expected one predictions file, found {f}")
        files.append(f[0])
    names = {os.path.basename(f) for f in files}
    if len(names) != 1:
        raise SystemExit(f"members predict different sets: {names}")
    os.makedirs(out, exist_ok=True)
    name = names.pop()
    streams = [gzip.open(f, "rt") for f in files]
    n = 0
    with gzip.open(os.path.join(out, name), "wt") as o:
        for lines in zip(*streams):
            recs = [json.loads(x) for x in lines]
            cid = recs[0]["candidate_id"]
            if any(r["candidate_id"] != cid for r in recs):
                raise SystemExit(f"record {n}: candidate ids differ: {[r['candidate_id'] for r in recs]}")
            logp = {c: sum(math.log(max(r[c], 1e-6)) for r in recs) / len(recs) for c in CLASSES}
            mx = max(logp.values())
            z = sum(math.exp(v - mx) for v in logp.values())
            fused = {c: math.exp(logp[c] - mx) / z for c in CLASSES}
            rec = dict(recs[0])
            for c in CLASSES:
                rec[c] = round(fused[c], 5)
            rec["pred"] = NAMES[max(CLASSES, key=lambda c: fused[c])]
            o.write(json.dumps(rec) + "\n")
            n += 1
    for s in streams:
        if s.readline():
            raise SystemExit("members have different numbers of records")
        s.close()
    metrics = files[0][:-len(".predictions.ndjson.gz")] + ".metrics.json"
    shutil.copy(metrics, os.path.join(out, os.path.basename(metrics)))
    with open(os.path.join(out, "fused_from.json"), "w") as f:
        json.dump({"members": members, "rule": "normalized geometric mean of class probabilities"}, f, indent=1)
    print(f"{out}: {n:,} records fused from {len(members)} members")


if __name__ == "__main__":
    main()

"""E9j step 4: the COLO829T off-reference somatic SNV keys (s3_classify.py, offref_somatic_<graph>.tsv) of one tensor
set scored by that set's in-sample SNV model, three ways (as E1.1 section 5 / asm_scores.py, here with the current
predict code so the normalization is the checkpoint's own, frozen copy code_8580a9f):
  base      the tensors as they are; called at the checkpoint's validation recall-0.9 threshold (graph_vcf's PASS)
  rescue    predict's off-reference rescue (OffrefRescue, --offref-site90 10: ch6 cells below the every-haplotype code
            raised to it in the site columns +- 10), called by argmax (graph_vcf's offref_call)
  all_H     every valid ch6 cell set to the every-haplotype code ("every haplotype runs through this node"), at the
            threshold and by argmax
Keys are matched to tensors by exact candidate id over every autosome of the set. Normalization as predict.py: a
norm "dataset" checkpoint uses its training statistics for its own training set (all sets here are in-sample; the
checkpoint records the set by its RIS path, so it is matched by the set name).
GPU (tequila). usage: python s4_scores.py --run RUN_DIR --set SET --keys offref_somatic_<g>.tsv --out OUT.json
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

CODE = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9j_offref_colo829t/code_8580a9f")
sys.path.insert(0, str(CODE / "machine_learning"))
from pansoma_net_v2 import metrics  # noqa: E402
from pansoma_net_v2.data import (INDEX_VERSION, PATH_COUNT_ALL, EpochSampler, KindIndex, TensorDataset,  # noqa: E402
                                 _labels_meta)
from pansoma_net_v2.encode import DATASET_NORM  # noqa: E402
from pansoma_net_v2.model import PansomaNetV2  # noqa: E402
from pansoma_net_v2.predict import PATH_COUNT, OffrefRescue  # noqa: E402
from pansoma_net_v2.train import dataset_stats, make_loader, predict_probs  # noqa: E402

TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
SHARED = Path("/scratch/jshen/data/pansoma_net_v2_runs/index_cache")


class AllH(TensorDataset):
    """Every valid ch6 cell set to `floor` (the every-haplotype code)."""

    def __init__(self, parts, floor, labels="label"):
        super().__init__(parts, labels)
        self.floor = floor

    def __getitem__(self, i):
        x, blocks, scalars, label, part = super().__getitem__(i)
        pc = x[PATH_COUNT]
        pc[pc > 0] = self.floor
        return x, blocks, scalars, label, part


def open_index(d, own):
    """KindIndex from the shared cache when it is current (read only), else built into own."""
    import hashlib
    d = Path(d).resolve()
    name = f"{d.parent.parent.name}.{d.parent.name}.{d.name}.{hashlib.sha1(str(d).encode()).hexdigest()[:12]}"
    npz = SHARED / f"{name}.npz"
    if npz.exists() and (SHARED / f"{name}.candidates.txt").exists():
        with np.load(npz) as z:
            meta = json.loads(str(z["meta"]))
        if meta.get("index_version") == INDEX_VERSION and meta["labels"] == _labels_meta(d) and meta["directory"] == str(d):
            return KindIndex(d, SHARED)
    return KindIndex(d, own)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--keys", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    ck = torch.load(Path(args.run) / "best.pth", map_location="cpu", weights_only=False)
    model = PansomaNetV2.from_checkpoint(ck).to(device).eval()
    amp = device.type == "cuda"
    norm = model.config["norm"]
    if norm == "dataset" and not model.encoder.dataset_norm:
        model.encoder.dataset_norm = DATASET_NORM
    trained_stats = {d["directory"]: d["stats"] for d in ck.get("dataset_stats", [])}
    out = Path(args.out)
    index = open_index(TENSORS / args.set / "tensors" / "SNV", out.parent / "index_cache")
    part_stats, source = None, None
    if norm == "dataset":
        tail = f"/{args.set}/tensors/SNV"   # the checkpoint records the RIS path of its training set
        own = [k for k in trained_stats if k == str(index.dir) or k.endswith(tail)]
        if own:
            part_stats, source = trained_stats[own[0]], f"checkpoint ({own[0]})"
        else:
            part_stats = dataset_stats(index, index.select(None, labelled=False), 20000, 12345, 4, print, str(index.dir))
            source = "fitted"
    model.encoder.set_part_stats([part_stats] if part_stats else [])
    thr = float(ck["val"]["somatic_at_recall"]["0.9"]["threshold"])
    rows = [r for r in csv.DictReader(open(args.keys), delimiter="\t") if r["kind"] == "SNP"]
    by_key = {}
    for r in rows:
        by_key.setdefault(r["key"], r)
    cands = index.candidates()
    pos = np.array(sorted(p for p, c in enumerate(cands) if c in by_key), np.int64)
    print(f"{args.set}: {len(by_key)} SNP keys, {len(pos)} with a tensor; norm {norm} ({source}); threshold {thr:.5f}",
          flush=True)
    floor = PATH_COUNT_ALL[index.path_count]
    a, reasons, chroms = index.arrays, index.meta["reasons"], index.meta["chroms"]
    res = {}
    for mode, ds in (("base", TensorDataset([(index, pos)], labels="eval_label")),
                     ("rescue", OffrefRescue([(index, pos)], 10, floor, labels="eval_label")),
                     ("all_H", AllH([(index, pos)], floor, labels="eval_label"))):
        t0 = time.time()
        _, probs = predict_probs(model, make_loader(ds, EpochSampler(len(ds), shuffle=False), 4, 256, persistent=False),
                                 device, amp)
        res[mode] = probs
        print(f"  {mode}: {time.time() - t0:.0f}s", flush=True)
    keys_out = []
    for k, p in enumerate(pos.tolist()):
        r = by_key[cands[p]]
        keys_out.append(dict(key=cands[p], event=r["event"], grch38=r["grch38_records"], chrom=chroms[a["chrom"][p]],
                             label=int(a["label"][p]), reason=reasons[a["reason"][p]], off_reference=bool(a["off_reference"][p]),
                             **{f"p_{m}": round(float(res[m][k][metrics.SOMATIC]), 5) for m in res},
                             **{f"argmax_{m}": int(res[m][k].argmax()) for m in res}))
    summ = {}
    for m in res:
        pr = res[m]
        chr1 = np.array([r["chrom"] == "chr1" for r in keys_out], bool)
        summ[m] = dict(ge_thr=int((pr[:, metrics.SOMATIC] >= thr).sum()), argmax=int((pr.argmax(1) == metrics.SOMATIC).sum()),
                       chr1_ge_thr=int((pr[chr1, metrics.SOMATIC] >= thr).sum()) if chr1.any() else 0,
                       chr1_argmax=int((pr[chr1].argmax(1) == metrics.SOMATIC).sum()) if chr1.any() else 0)
    out.write_text(json.dumps(dict(run=args.run, set=args.set, keys_table=args.keys, threshold=thr, norm=norm,
                                   norm_source=source, n_keys=len(by_key), n_with_tensor=int(len(pos)),
                                   chr1_with_tensor=int(sum(r["chrom"] == "chr1" for r in keys_out)),
                                   summary=summ, keys=keys_out), indent=1) + "\n")
    print(json.dumps(summ), flush=True)


if __name__ == "__main__":
    main()

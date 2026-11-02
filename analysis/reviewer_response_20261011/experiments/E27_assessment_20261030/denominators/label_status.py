#!/usr/bin/env python3
"""Per-truth-allele label status on the HPRC v2.1 (d46) tensor sets of HG008T / COLO829T (R1.6 / R1.8 fact sheet,
2026-10-30). Read-only: the E21 index caches (pansoma_net_v2 KindIndex npz, built from the sets' current labels;
checked against labels.manifest.json before use, never rebuilt), the sets' somatic.recall.tsv, the E5b per-truth PoN
tags, the somatic BED named in labels.manifest.json.

Per set x kind, for the PASS somatic truth alleles of the kind inside the truth BED:
  exact      a tensor labelled 1 by representative_allele_is_somatic_truth stands for it (truth_ids)
  partial    no exact tensor, but a tensor labelled 1 by allele_partial_somatic_truth / residual_partial_somatic_truth
             (left out of training and validation by every final model's --ignore-reasons; INDEL tests count them)
  no label-1 tensor: by somatic.recall.tsv status (tensor_non_representative_allele / filtered / no_candidate)
Scopes: chr1 (the test chromosome; joined with the PoN tag) and chr2-22 (all blocks, and the training blocks of the
train.py split for seeds 0 / 1 / 2: block_split, --val-fraction 0.05, --val-block-nodes 20000).
Also: tensors of chr2-22 (training blocks, seed 0) by label and reason, by off_reference; the -1 tensors with reason
outside_confident_region inside the somatic BED (trained as 0 by the final INDEL models, --non-reasons-in-bed).
Writes label_status.json next to this script (tables are made by the fact sheet).
"""
import csv
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning")
from pansoma_net_v2.data import _labels_meta, block_split, validation_truth  # noqa: E402

HERE = Path(__file__).resolve().parent
T = Path("/scratch/jshen/data/pansoma_v2_tensors")
CACHE = Path("/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/index_cache")
PT = "/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46/per_truth.tsv.gz"
SETS = [("HG008T", p) for p in ("Illumina", "PacBio", "ONT")] + [("COLO829T", p) for p in ("Illumina", "fiberseq", "ONT")]
KINDS = {"SNV": ("SNP",), "INDEL": ("INS", "DEL")}
EXACT = "representative_allele_is_somatic_truth"
PARTIAL = ("allele_partial_somatic_truth", "residual_partial_somatic_truth")
G_EXACT = ("representative_allele_is_germline_truth", "representative_allele_is_germline_truth_gap_filtered")
G_PART = ("allele_partial_germline_truth", "residual_partial_germline_truth")
TRAIN = [f"chr{i}" for i in range(2, 23)]


def load(set_name, kind):
    d = (T / set_name / "tensors" / kind).resolve()
    hits = sorted(CACHE.glob(f"{set_name}.tensors.{kind}.*.npz"))
    assert len(hits) == 1, (set_name, kind, hits)
    z = np.load(hits[0])
    meta = json.loads(str(z["meta"]))
    assert meta["directory"] == str(d), (meta["directory"], d)
    assert meta["labels"] == _labels_meta(d), f"{set_name} {kind}: cache labels differ from the current labels"
    a = {k: z[k] for k in z.files if k != "meta"}
    return SimpleNamespace(arrays=a, meta=meta, dir=d)


def somatic_bed(d):
    m = json.loads((d / "labels.manifest.json").read_text())
    path = m["truth"]["somatic"]["bed"]
    iv = defaultdict(list)
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        for line in f:
            x = line.split("\t")
            if len(x) >= 3 and not line.startswith(("#", "track", "browser")):
                iv[x[0]].append((int(x[1]), int(x[2])))
    out = {}
    for c, items in iv.items():
        items.sort()
        merged = []
        for s, e in items:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = (np.array([m[0] for m in merged]), np.array([m[1] for m in merged]))
    return path, out


def in_bed(index, bed):
    a, names = index.arrays, index.meta["chroms"]
    inside = np.zeros(len(a["label"]), bool)
    for c, name in enumerate(names):
        if name not in bed:
            continue
        sel = np.flatnonzero((a["chrom"] == c) & (a["pos0"] >= 0))
        starts, ends = bed[name]
        k = np.searchsorted(starts, a["pos0"][sel], side="right") - 1
        ok = k >= 0
        ok[ok] = ends[k[ok]] > a["pos0"][sel][ok]
        inside[sel[ok]] = True
    return inside


def main():
    pon = {}
    for r in csv.DictReader(gzip.open(PT, "rt"), delimiter="\t"):
        pon[(r["sample"], int(r["truth_id"]))] = bool(r["pon"])
    out = {}
    only = sys.argv[1:]  # optional "SET|KIND" filters (test runs); output then goes to label_status.test.json
    for sample, plat in SETS:
        arm = os.environ.get("ARM", "d46")  # d46 (default) / d9 (no suffix) / linear
        set_name = f"{sample}_{plat}" + ("" if arm == "d9" else f"_{arm}")
        rec = list(csv.DictReader(open(T / set_name / "tensors" / "somatic.recall.tsv"), delimiter="\t"))
        for kind in ("SNV", "INDEL"):
            if only and f"{set_name}|{kind}" not in only:
                continue
            idx = load(set_name, kind)
            a, reasons, names = idx.arrays, idx.meta["reasons"], idx.meta["chroms"]
            rname = np.array(reasons)[a["reason"]]
            chrom = np.array(names)[a["chrom"]]
            lab = a["label"]
            off = a["off_reference"]
            # per truth id: tensors labelled 1 with their reason / off_reference / position
            som = np.flatnonzero(lab == 1)
            by_truth = defaultdict(list)
            for p in som.tolist():
                for t in a["truth_ids"][a["truth_ptr"][p]:a["truth_ptr"][p + 1]].tolist():
                    by_truth[int(t)].append(p)
            # sanity: exact matches off the GRCh38 path
            exact_off = int(((lab == 1) & (rname == EXACT) & off).sum())
            partial_off = int(((lab == 1) & np.isin(rname, PARTIAL) & off).sum())
            partial_on = int(((lab == 1) & np.isin(rname, PARTIAL) & ~off).sum())
            # training-block masks per seed (chr2-22, all tensors)
            pos_train = np.flatnonzero(np.isin(a["chrom"], [i for i, c in enumerate(names) if c in TRAIN]))
            in_train = {}
            val_truth = {}
            truth_train = {int(r["truth_id"]): (r["chrom"], int(r["vcf_pos"]) - 1) for r in rec
                           if r["passed"] == "True" and r["in_bed"] == "True" and r["kind"] in KINDS[kind]
                           and r["chrom"] in TRAIN}
            for seed in (0, 1, 2):
                tr, va = block_split(idx, pos_train, 0.05, 20000, seed)
                m = np.zeros(len(lab), bool)
                m[tr] = True
                in_train[seed] = m
                val_truth[seed] = validation_truth(idx, pos_train, truth_train, 0.05, 20000, seed)
            res = {"checks": dict(exact_somatic_off_reference=exact_off, partial_somatic_off_reference=partial_off,
                                  partial_somatic_on_grch38=partial_on)}
            for scope in ("chr1", "chr2-22"):
                chroms = ["chr1"] if scope == "chr1" else TRAIN
                truth = [r for r in rec if r["passed"] == "True" and r["in_bed"] == "True" and r["kind"] in KINDS[kind]
                         and r["chrom"] in chroms]
                c = Counter()
                c_untag = Counter()
                c_seed = {s: Counter() for s in (0, 1, 2)}
                for r in truth:
                    tid = int(r["truth_id"])
                    ps = by_truth.get(tid, [])
                    ex = [p for p in ps if rname[p] == EXACT]
                    pa = [p for p in ps if rname[p] in PARTIAL]
                    if ex:
                        cls = "exact"
                    elif pa:
                        cls = "partial_only_offref" if all(off[p] for p in pa) else "partial_only_grch38"
                    else:
                        st = r["status"]
                        cls = {"tensor_representative": "rep_but_no_label1(?)", "tensor_non_representative_allele":
                               "other_allele_of_site", "filtered": "filtered", "no_candidate": "no_candidate"}.get(st, st)
                    c[cls] += 1
                    tagged = pon.get((sample, tid))
                    c["_tagged"] += bool(tagged)
                    c["_pon_missing"] += tagged is None
                    if not tagged:
                        c_untag[cls] += 1
                    if scope == "chr2-22":
                        for s in (0, 1, 2):
                            if tid in val_truth[s]:
                                c_seed[s]["in_validation_region"] += 1
                                continue
                            c_seed[s]["training_region"] += 1
                            if any(in_train[s][p] for p in ex):
                                c_seed[s]["exact_training_example"] += 1
                            elif any(in_train[s][p] for p in pa):
                                c_seed[s]["partial_only_training_block"] += 1
                c["n"] = len(truth)
                res[scope] = dict(all=dict(c), untagged=dict(c_untag))
                if scope == "chr2-22":
                    res[scope]["seeds"] = {s: dict(v) for s, v in c_seed.items()}
                    # exact positives outside the BED or not PASS (training uses them; not in the denominators)
                    outside = [r for r in rec if r["kind"] in KINDS[kind] and r["chrom"] in TRAIN
                               and not (r["passed"] == "True" and r["in_bed"] == "True")]
                    n_out = sum(any(rname[p] == EXACT for p in by_truth.get(int(r["truth_id"]), [])) for r in outside)
                    res[scope]["truth_outside_bed_or_not_pass"] = dict(n=len(outside), with_exact_label1=n_out)
            # tensors of chr2-22, seed-0 training blocks, by label x reason x off_reference
            m0 = in_train[0]
            tab = Counter()
            for l, rn, o in zip(lab[m0].tolist(), rname[m0].tolist(), off[m0].tolist()):
                tab[f"{l}|{rn}|{'off' if o else 'grch38'}"] += 1
            res["train_blocks_seed0_tensors"] = dict(tab)
            # outside_confident_region inside the somatic BED (final INDEL models: trained as 0), chr2-22 seed-0 train
            path, bed = somatic_bed(idx.dir)
            inb = in_bed(idx, bed)
            ocr = (rname == "outside_confident_region")
            res["outside_confident_region_seed0_train"] = dict(somatic_bed=path, total=int((ocr & m0).sum()),
                                                               inside_somatic_bed=int((ocr & m0 & inb).sum()))
            out[f"{set_name}|{kind}"] = res
            print(set_name, kind, json.dumps(res["checks"]), json.dumps(res["chr1"]["all"]), flush=True)
    suffix = "" if os.environ.get("ARM", "d46") == "d46" else "." + os.environ["ARM"]
    (HERE / ("label_status.test.json" if only else f"label_status{suffix}.json")).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()

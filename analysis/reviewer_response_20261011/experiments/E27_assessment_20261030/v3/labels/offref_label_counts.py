#!/usr/bin/env python3
"""R1.8 / Q3: off-reference (off-GRCh38-path, summary grch38 null) tensors of the six HPRC v2.1 d46 HG008T / COLO829T
sets, chr2-22, by label reason, in the seed-0 training and validation blocks of train.py (block_split, 0.05, 20000).
Read-only: the E21 index caches (KindIndex npz, checked against the sets' current labels.manifest.json).

Per set x kind:
  offref_by_reason[train|val]            label reason of every off-reference tensor
  how each final model uses them         training label (SNV: --ignore-reasons partial somatic; INDEL release: also
                                         --non-reasons-in-bed outside_confident_region = 0 where the anchor middle is
                                         inside the somatic BED) and validation label (eval_label: off_reference_no_truth_match -> 0)
  step4                                  the off-reference tensors labelled 0 by truth_matches_non_representative_allele
                                         (the letter's "83-717"), split by which truth set holds the site's other allele
                                         and by whether A1 is a multi-node deletion ('@' key), from <chrom>_labels.ndjson
"""
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, "/scratch/jshen/Github/Pansoma/machine_learning")
from pansoma_net_v2.data import _labels_meta, block_split  # noqa: E402
from pansoma_net_v2.bed import Bed  # noqa: E402

T = Path("/scratch/jshen/data/pansoma_v2_tensors")
CACHE = Path("/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/index_cache")
SETS = [f"HG008T_{p}_d46" for p in ("Illumina", "PacBio", "ONT")] + [f"COLO829T_{p}_d46" for p in ("Illumina", "fiberseq", "ONT")]
TRAIN = [f"chr{i}" for i in range(2, 23)]
IGNORE = {"residual_partial_somatic_truth", "allele_partial_somatic_truth"}
STEP4 = "truth_matches_non_representative_allele"


class Idx:
    def __init__(self, set_name, kind):
        self.dir = (T / set_name / "tensors" / kind).resolve()
        hits = sorted(CACHE.glob(f"{set_name}.tensors.{kind}.*.npz"))
        assert len(hits) == 1, (set_name, kind, hits)
        z = np.load(hits[0])
        self.meta = json.loads(str(z["meta"]))
        assert self.meta["directory"] == str(self.dir), (self.meta["directory"], self.dir)
        assert self.meta["labels"] == _labels_meta(self.dir), f"{set_name} {kind}: cache labels differ"
        self.arrays = {k: z[k] for k in z.files if k != "meta"}


def in_somatic_bed(idx):
    """train.in_somatic_bed: pos0 (off-reference: the anchor's middle base) inside the somatic truth BED."""
    m = json.loads((idx.dir / "labels.manifest.json").read_text())
    bed = Bed.read(m["truth"]["somatic"]["bed"])
    a, names = idx.arrays, idx.meta["chroms"]
    inside = np.zeros(len(a["label"]), bool)
    for c, name in enumerate(names):
        if name not in bed.data:
            continue
        sel = np.flatnonzero((a["chrom"] == c) & (a["pos0"] >= 0))
        starts, ends = bed.data[name]
        k = np.searchsorted(starts, a["pos0"][sel], side="right") - 1
        ok = k >= 0
        ok[ok] = ends[k[ok]] > a["pos0"][sel][ok]
        inside[sel[ok]] = True
    return inside


def step4_detail(idx):
    """Off-reference step-4 zeros on chr2-22 (all blocks), from labels.ndjson: which truth set the site's other allele
    is in, and whether A1 is a multi-node deletion."""
    c = Counter()
    for chrom in TRAIN:
        p = idx.dir / f"{chrom}_labels.ndjson"
        with open(p) as f:
            for line in f:
                if STEP4 not in line or '"grch38": null' not in line:
                    continue
                r = json.loads(line)
                if r["reason"] != STEP4 or r.get("grch38") is not None:
                    continue
                which = ("somatic" if r.get("somatic") else "") + ("+" if r.get("somatic") and r.get("germline") else "") \
                    + ("germline" if r.get("germline") else "")
                c[f"site mate in {which or 'none'}"] += 1
                c["A1 multi-node DEL (@)" if "@" in r["candidate_id"] else "A1 single-node"] += 1
                c["total"] += 1
    return dict(c)


def main():
    out = {}
    for s in SETS:
        for kind in ("SNV", "INDEL"):
            idx = Idx(s, kind)
            a, reasons = idx.arrays, idx.meta["reasons"]
            names = idx.meta["chroms"]
            on_train = np.isin(a["chrom"], [i for i, c in enumerate(names) if c in TRAIN])
            pos = np.flatnonzero(on_train)
            tr, va = block_split(idx, pos, 0.05, 20000, 0)
            off = a["off_reference"]
            rname = np.array(reasons)[a["reason"]]
            sbed = in_somatic_bed(idx)
            rec = {}
            for name, blk in (("train", tr), ("val", va)):
                o = blk[off[blk]]
                rec[f"offref_by_reason_{name}"] = dict(Counter(rname[o].tolist()))
                rec[f"offref_total_{name}"] = int(len(o))
                rec[f"all_tensors_{name}"] = int(len(blk))
                # training labels of the final models
                lab = a["label"][o]
                ign = np.isin(rname[o], list(IGNORE))
                snv_train = (lab >= 0) & ~ign
                rec[f"offref_used_snvrecipe_{name}"] = dict(Counter(lab[snv_train].tolist()))
                occ = (rname[o] == "outside_confident_region") & sbed[o]
                rec[f"offref_outside_conf_in_somatic_bed_{name}"] = int(occ.sum())
                ev = a["eval_label"][o]
                rec[f"offref_eval_label_{name}"] = dict(Counter(ev.tolist()))
                rec[f"offref_step4_zero_{name}"] = int(((rname[o] == STEP4) & (lab == 0)).sum())
                rec[f"offref_no_truth_match_{name}"] = int((rname[o] == "off_reference_no_truth_match").sum())
                rec[f"offref_partial_somatic_{name}"] = int(ign.sum())
                rec[f"non_eval_{name}"] = int((a["eval_label"][blk] == 0).sum())
            if kind == "INDEL":
                rec["step4_detail_chr2_22_all_blocks"] = step4_detail(idx)
            out[f"{s}.{kind}"] = rec
            print(s, kind, json.dumps(rec), flush=True)
    Path(__file__).with_suffix(".json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()

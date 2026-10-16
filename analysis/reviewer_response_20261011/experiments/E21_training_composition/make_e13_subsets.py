"""E13 input (prepared in E21): nested random subsets of the training node blocks of one set, as
--exclude-candidates files for machine_learning/pansoma_net_v2/train.py (no code change needed).

train.py splits chr2-22 into 20,000-node blocks and sends 5 % of them to validation (data.block_split, seed 0). This
script keeps the validation blocks and chr1 untouched and keeps a fraction f of the *training* blocks: block (chrom, b)
is kept when u_r(chrom, b) < f, with u_r a SHA-1 uniform of ("e13", replicate r, chrom, b), so for one replicate the
subsets are nested (10 % within 25 % within 50 % within 75 %). The candidate ids of every labelled tensor (label >= 0)
in a dropped training block go to the exclude file; train.py sets them to -1 (as if never built), so training sees only
the kept blocks, validation and its PASS threshold are those of the full run, and predict (whose default exclude file
is the checkpoint's) excludes nothing on chr1.

    python make_e13_subsets.py SET KIND OUTDIR [fractions ...] [--replicates N]
    python make_e13_subsets.py HG008T_PacBio_d46 SNV /scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/e13 0.10 0.25 0.50 0.75 --replicates 2

Writes OUTDIR/<SET>.<KIND>.keep<pct>.rep<r>.exclude.txt and prints / writes e13_subsets_<SET>_<KIND>.tsv (here) with
the training tensors per class and the distinct somatic truth alleles left in each subset.
"""
import argparse
import csv
import hashlib
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from summarize import CACHE, IGNORE_REASONS, TENSORS, TRAIN_CHROMS, BLOCK_NODES, VAL_FRACTION, SEED  # noqa: E402
from summarize import KindIndex, block_split  # noqa: E402


def u(rep, chrom, block):
    return int.from_bytes(hashlib.sha1(f"e13:{rep}:{chrom}:{block}".encode()).digest()[:8], "big") / 2 ** 64


def main():
    p = argparse.ArgumentParser()
    p.add_argument("set")
    p.add_argument("kind")
    p.add_argument("outdir")
    p.add_argument("fractions", nargs="*", type=float, default=[0.10, 0.25, 0.50, 0.75])
    p.add_argument("--replicates", type=int, default=1)
    args = p.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    index = KindIndex(TENSORS / args.set / "tensors" / args.kind, cache_dir=CACHE)
    a, names = index.arrays, index.meta["chroms"]
    pos = index.select(TRAIN_CHROMS, labelled=False)
    tr_all, va_all = block_split(index, pos, VAL_FRACTION, BLOCK_NODES, SEED)
    codes = [k for k, r in enumerate(index.meta["reasons"]) if r in IGNORE_REASONS]
    lab = tr_all[a["label"][tr_all] >= 0]
    train = lab[~np.isin(a["reason"][lab], codes)]
    cands = index.candidates()
    # candidate ids must not repeat outside the training blocks (else an exclusion would reach validation / chr1)
    outside = set(cands[i] for i in np.setdiff1d(np.arange(len(index)), tr_all))
    blocks = a["node"][lab] // BLOCK_NODES
    chrom = a["chrom"][lab]
    keys = np.unique(np.stack([chrom.astype(np.int64), blocks], 1), axis=0)
    rows = []
    full = dict(fraction=1.0, replicate=0, kept_blocks=len(keys), blocks=len(keys), excluded=0,
                **{c: int((a["label"][train] == k).sum()) for k, c in enumerate(("non", "somatic", "germline"))})
    ids = set()
    for p_ in train[a["label"][train] == 1]:
        ids.update(index.truth_of(p_).tolist())
    full["somatic_truth_alleles"] = len(ids)
    rows.append(full)
    for rep in range(1, args.replicates + 1):
        uk = {(int(c), int(b)): u(rep, names[c], int(b)) for c, b in keys.tolist()}
        ub = np.array([uk[(int(c), int(b))] for c, b in zip(chrom.tolist(), blocks.tolist())])
        for f in sorted(args.fractions):
            drop = lab[ub >= f]
            excl = [cands[i] for i in drop]
            clash = len(set(excl) & outside)
            if clash:
                raise SystemExit(f"{clash} excluded candidate ids also occur outside the training blocks")
            path = out / f"{args.set}.{args.kind}.keep{int(round(f * 100)):02d}.rep{rep}.exclude.txt"
            with open(path, "w") as fh:
                fh.write("candidate_id\n")
                fh.write("\n".join(excl) + "\n")
            kept = train[np.isin(train, drop, invert=True)]
            ids = set()
            for p_ in kept[a["label"][kept] == 1]:
                ids.update(index.truth_of(p_).tolist())
            rows.append(dict(fraction=f, replicate=rep, kept_blocks=sum(v < f for v in uk.values()), blocks=len(keys),
                             excluded=len(excl), **{c: int((a["label"][kept] == k).sum())
                                                    for k, c in enumerate(("non", "somatic", "germline"))},
                             somatic_truth_alleles=len(ids), file=str(path)))
            print(rows[-1], flush=True)
    with open(HERE / f"e13_subsets_{args.set}_{args.kind}.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[-1]), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in rows[-1]})


if __name__ == "__main__":
    main()

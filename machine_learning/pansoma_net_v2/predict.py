"""Class probabilities of a trained PansomaNetV2 for every tensor of merged tensor sets, and the test metrics.

    cd machine_learning
    python -m pansoma_net_v2.predict --checkpoint <run>/best.pth --tensors <set> [...] --output <dir> --chroms chr1

Writes, per <sample>.<set>.<KIND>:
- .predictions.ndjson.gz: one JSON record per tensor of the chromosomes, in index order: chrom, candidate_id,
  label (the truth label), test_label (the label the test scores, null when not in the test), in_test, reason,
  off_reference (the node is off the GRCh38 path), truth_ids (the somatic truth alleles a tensor labelled 1
  stands for: representative or partial match), p_non, p_somatic, p_germline, pred (with the somatic
  threshold: the checkpoint's, or --threshold).
- .metrics.json: metrics.report over the test tensors. These are the tensors labelled 0/1/2 inside the
  confident region, and the off-reference tensors without a truth match (as non: test-time calling meets them
  inside the BED). Tensors outside the BED (labels 1 and 2 included), below the AF floors, or without a GRCh38
  position are left out, since calling drops them without truth.
  Also `truth`: metrics.truth_report against the somatic truth VCF (the set's somatic.recall.tsv: PASS, in the
  BED, the chromosomes, this kind), at the same threshold. Every truth allele counts once: it is found when its
  best tensor (a labelled-1 tensor matching it, exactly or partially) is called somatic, so duplicate tensors of
  one truth do not count twice; truth alleles without a tensor are misses; a false positive is a somatic call on
  a tensor labelled 0 or 2. Reasons in --ignore-reasons (default: those the model was trained without) are left out
  like -1; their records say ignored: true and keep label, reason and truth_ids. A found truth allele of the other kind (an SNV tensor matching an INDEL truth
  partially) is a true call too, once. `ceiling` is the recall of a perfect model. `combine` scores the SNV and
  INDEL predictions of a sample together against the whole truth VCF.

The encoder uses the checkpoint's statistics; they are not refitted on these tensors.
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import torch

from . import metrics
from .data import CLASSES, KINDS, EpochSampler, TensorDataset, load_parts, somatic_truth
from .model import PansomaNetV2
from .train import describe, make_loader, predict_probs, retry_workers


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--tensors", nargs="+", required=True)
    p.add_argument("--kinds", nargs="+", default=list(KINDS), choices=KINDS)
    p.add_argument("--chroms", nargs="+", help="default: all chromosomes of the sets")
    p.add_argument("--output", required=True)
    p.add_argument("--cache-dir", help="index cache (default: <output>/index_cache)")
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--num-workers", type=int, default=8)
    p.add_argument("--amp", choices=["bf16", "off"], default="bf16")
    p.add_argument("--threshold", type=float, help="somatic threshold (default: the checkpoint's)")
    p.add_argument("--ignore-reasons", nargs="*",
                   help="label reasons left out of the test, as if -1 (default: the checkpoint's train --ignore-reasons); "
                        "their records keep label, reason and truth_ids and get ignored: true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model = PansomaNetV2.from_checkpoint(checkpoint).to(device).eval()
    amp = args.amp == "bf16" and device.type == "cuda"
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    print(f"checkpoint {args.checkpoint} (epoch {checkpoint['epoch']}), statistics {json.dumps(model.encoder.stats())}",
          flush=True)
    for index, positions in load_parts(args.tensors, args.kinds, args.cache_dir or str(out / "index_cache"),
                                       chroms=args.chroms, labelled=False):
        dataset = TensorDataset([(index, positions)], labels="eval_label")
        ignore = args.ignore_reasons if args.ignore_reasons is not None else checkpoint.get("args", {}).get("ignore_reasons", [])
        ignored = np.isin(index.arrays["reason"][positions],
                          [k for k, r in enumerate(index.meta["reasons"]) if r in set(ignore)])
        dataset.label[ignored] = -1  # out of the test like the other -1; the record keeps its label and reason
        labels, probs = retry_workers(
            lambda: predict_probs(model, make_loader(dataset, EpochSampler(len(dataset), shuffle=False), args.num_workers,
                                                     args.batch_size, persistent=False), device, amp),
            lambda m: print(m, flush=True), f"{index.dir}")
        threshold = args.threshold if args.threshold is not None else checkpoint.get("somatic_threshold")
        pred = metrics.call(probs, threshold)
        name = f"{index.dir.parent.parent.name}.{index.dir.parent.name}.{index.kind}"
        a, chroms, candidates = index.arrays, index.meta["chroms"], index.candidates()
        with gzip.open(out / f"{name}.predictions.ndjson.gz", "wt") as f:
            for k, pos in enumerate(positions):
                f.write(json.dumps(dict(
                    chrom=chroms[a["chrom"][pos]], candidate_id=candidates[pos], label=int(a["label"][pos]),
                    test_label=int(labels[k]) if labels[k] >= 0 else None, in_test=bool(labels[k] >= 0),
                    reason=index.reason(pos), ignored=bool(ignored[k]), off_reference=bool(a["off_reference"][pos]),
                    truth_ids=index.truth_of(pos).tolist(),
                    **{f"p_{c}": round(float(probs[k][j]), 5) for j, c in enumerate(CLASSES)},
                    pred=CLASSES[pred[k]])) + "\n")
        report = metrics.report(labels, probs, threshold)
        on = {index.meta["chroms"][c] for c in set(index.arrays["chrom"][positions].tolist())}
        truth = somatic_truth(index, chroms=on)
        truth_args = None
        if truth is not None:  # every truth allele of the chromosomes counts once; see metrics.truth_report
            truth_args = (set(truth), [set(index.truth_of(p).tolist()) for p in positions], somatic_truth(index, on, True))
            report["truth"] = metrics.truth_report(truth_args[0], truth_args[1], labels, probs, threshold,
                                                   other=truth_args[2])
        # three decision rules: the threshold t above, argmax, and the checkpoint's validation recall-0.9 threshold
        r09 = ((checkpoint.get("val") or {}).get("somatic_at_recall") or {}).get("0.9")
        report["rules"] = metrics.decision_rules(labels, probs, {"threshold_t": threshold, "argmax": None,
                                                                 "val_recall_0.9": r09["threshold"] if r09 else None},
                                                 truth_args)
        in_test = labels >= 0
        off = a["off_reference"][positions]
        report.update(all_tensors=int(len(labels)), left_out=int((~in_test).sum()),
                      off_reference_in_test=int((in_test & off).sum()),
                      off_reference_somatic_in_test=int((in_test & off & (labels == 1)).sum()),
                      ignored_reasons=list(ignore), ignored=int(ignored.sum()),
                      directory=str(index.dir), labels=index.meta["labels"],
                      checkpoint=str(Path(args.checkpoint).resolve()), chroms=args.chroms)
        (out / f"{name}.metrics.json").write_text(json.dumps(report, indent=2) + "\n")
        print(f"{name}: {len(labels):,} tensors, {report['tensors']:,} in the test "
              f"({report['off_reference_in_test']:,} off-reference) | {describe(report)}", flush=True)
        for rule, r in report["rules"].items():
            t, tr = r["per_tensor"], r.get("truth")
            at = "" if r["threshold"] is None else f"p>={r['threshold']:.3f}"
            line = f"  {rule:<15} {at:<9} per tensor P {t['precision']:.3f} R {t['recall']:.3f} F1 {t['f1']:.3f} ({t['calls']:,} calls)"
            if tr:
                line += (f" | vs truth P {tr['precision']:.3f} R {tr['recall']:.3f} F1 {tr['f1']:.3f} "
                         f"({tr['tp']} of {tr['truth_alleles']}, {tr['fp']:,} false calls)")
            print(line, flush=True)


if __name__ == "__main__":
    main()

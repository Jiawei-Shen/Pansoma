"""Class probabilities of a trained PansomaNetV2 for every tensor of merged tensor sets, and the test metrics.

    cd machine_learning
    python -m pansoma_net_v2.predict --checkpoint <run>/best.pth --tensors <set> [...] --output <dir> --chroms chr1

Writes, per <sample>.<set>.<KIND>:
- .predictions.ndjson.gz: one JSON record per tensor of the chromosomes, in index order: chrom, candidate_id,
  label (the truth label), test_label (the label the test scores, null when not in the test), in_test, reason,
  off_reference (the node is off the GRCh38 path), truth_ids (the somatic truth alleles a tensor labelled 1
  stands for: representative or partial match), p_non, p_somatic, p_germline, pred (with the somatic
  threshold: the checkpoint's, or --threshold). SNV off-reference tensors also get p_offref_non,
  p_offref_somatic, p_offref_germline and offref_call (see the off-reference rescue below).
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

Off-reference rescue (SNV, --offref-site90 W, default 10; 0 turns it off). The model learned that a low path count
(ch6, how many HPRC haplotypes run through the node) is not somatic, so it does not call somatic on an
off-reference node, where every read has a low path count. Each SNV off-reference tensor is scored again with its
ch6 cells below the every-haplotype code raised to it in the site columns (the middle column +- W; uncovered cells
stay 0): 90 on int8-count-linear100-log2 sets (HPRC v1.1, 90 haplotype paths), 100 on int8-count-haplotypes100
sets (data.PATH_COUNT_ALL). It is called somatic when somatic is the most probable class of those probabilities
(offref_call). graph_vcf uses offref_call for these tensors. `rules.pipeline` in .metrics.json scores that call on
the off-reference tensors and the checkpoint's validation recall-0.9 threshold on the others; `offref_rescue`
counts the tensors and calls and gives the floor.

The encoder uses the checkpoint's statistics; they are not refitted on these tensors. If any set codes channel 6
otherwise than the checkpoint's training sets (checkpoint path_count_storage; int8-count-linear100-log2 when it
records none), predict exits before predicting anything (a model trained with --drop-planes path_count never reads
channel 6 and takes sets of any storage); a set built with another --haplotypes is predicted with a
note, and .metrics.json `path_count` gives the set's storage and H and the checkpoint's H.
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import torch

from . import metrics
from .data import (CHECKPOINT_PATH_COUNT, CLASSES, KINDS, PATH_COUNT_ALL, TENSOR_SHAPE, EpochSampler, TensorDataset,
                   load_parts, somatic_truth)
from .encode import CONTINUOUS
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
    p.add_argument("--offref-site90", type=int, default=10, metavar="W",
                   help="SNV off-reference rescue: score the off-reference tensors again with ch6 below the "
                        "every-haplotype code (90 on int8-count-linear100-log2 sets, 100 on int8-count-haplotypes100) "
                        "raised to it in the site columns +- W and call them by argmax (0: off)")
    p.add_argument("--ignore-reasons", nargs="*",
                   help="label reasons left out of the test, as if -1 (default: the checkpoint's train --ignore-reasons); "
                        "their records keep label, reason and truth_ids and get ignored: true")
    return p.parse_args(argv)


PATH_COUNT = dict((name, channel) for name, channel, _ in CONTINUOUS)["path_count"]
SITE = TENSOR_SHAPE[2] // 2


class OffrefRescue(TensorDataset):
    """Tensors with the ch6 (path count) cells below `floor` set to `floor` in the columns SITE +- halfwidth."""

    def __init__(self, parts, halfwidth, floor, labels="label"):
        super().__init__(parts, labels)
        self.halfwidth, self.floor = halfwidth, floor

    def __getitem__(self, i):
        x, blocks, scalars, label = super().__getitem__(i)
        window = x[PATH_COUNT, :, SITE - self.halfwidth:SITE + self.halfwidth + 1]
        window[(window > 0) & (window < self.floor)] = self.floor
        return x, blocks, scalars, label


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
    trained = checkpoint.get("path_count_storage", CHECKPOINT_PATH_COUNT)
    trained_haplotypes = checkpoint.get("path_count_haplotypes") or []
    reads_path_count = "path_count" not in model.config["drop_planes"]  # the only plane made from channel 6
    parts = load_parts(args.tensors, args.kinds, args.cache_dir or str(out / "index_cache"), chroms=args.chroms,
                       labelled=False)
    for index, _ in parts:  # before any prediction: every set must code channel 6 like the training sets
        if reads_path_count and index.path_count != trained:
            raise SystemExit(f"{index.dir}: channel 6 is {index.path_count}, the checkpoint was trained on {trained}")
    for index, positions in parts:
        if reads_path_count and index.haplotypes is not None and index.haplotypes not in trained_haplotypes:
            print(f"note: {index.dir} was built with --haplotypes {index.haplotypes}, the checkpoint's training sets "
                  f"with {trained_haplotypes or 'none recorded'}", flush=True)
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
        rescue = index.kind == "SNV" and args.offref_site90 > 0
        off_k = np.flatnonzero(a["off_reference"][positions]) if rescue else np.array([], np.int64)
        off_probs = {}
        if len(off_k):
            edited = OffrefRescue([(index, positions[off_k])], args.offref_site90, PATH_COUNT_ALL[index.path_count])
            _, probs2 = retry_workers(
                lambda: predict_probs(model, make_loader(edited, EpochSampler(len(edited), shuffle=False), args.num_workers,
                                                         args.batch_size, persistent=False), device, amp),
                lambda m: print(m, flush=True), f"{index.dir} off-reference rescue")
            off_probs = dict(zip(off_k.tolist(), probs2))
        with gzip.open(out / f"{name}.predictions.ndjson.gz", "wt") as f:
            for k, pos in enumerate(positions):
                f.write(json.dumps(dict(
                    chrom=chroms[a["chrom"][pos]], candidate_id=candidates[pos], label=int(a["label"][pos]),
                    test_label=int(labels[k]) if labels[k] >= 0 else None, in_test=bool(labels[k] >= 0),
                    reason=index.reason(pos), ignored=bool(ignored[k]), off_reference=bool(a["off_reference"][pos]),
                    truth_ids=index.truth_of(pos).tolist(),
                    **{f"p_{c}": round(float(probs[k][j]), 5) for j, c in enumerate(CLASSES)},
                    pred=CLASSES[pred[k]],
                    **({f"p_offref_{c}": round(float(off_probs[k][j]), 5) for j, c in enumerate(CLASSES)}
                       | {"offref_call": bool(off_probs[k].argmax() == metrics.SOMATIC)} if k in off_probs else {}))) + "\n")
        report = metrics.report(labels, probs, threshold)
        on = {index.meta["chroms"][c] for c in set(index.arrays["chrom"][positions].tolist())}
        truth = somatic_truth(index, chroms=on)
        truth_args = None
        if truth is not None:  # every truth allele of the chromosomes counts once; see metrics.truth_report
            truth_args = (set(truth), [set(index.truth_of(p).tolist()) for p in positions], somatic_truth(index, on, True))
            report["truth"] = metrics.truth_report(truth_args[0], truth_args[1], labels, probs, threshold,
                                                   other=truth_args[2])
        # three decision rules: the threshold t above, argmax, and the checkpoint's validation recall-0.9 threshold
        # and the pipeline's calls: the validation recall-0.9 threshold, off-reference SNV by the rescue's argmax
        r09 = ((checkpoint.get("val") or {}).get("somatic_at_recall") or {}).get("0.9")
        rules, preds = {"threshold_t": threshold, "argmax": None, "val_recall_0.9": r09["threshold"] if r09 else None}, {}
        if r09:
            preds["pipeline"] = metrics.call(probs, r09["threshold"])
            if off_probs:
                preds["pipeline"][off_k] = np.stack([off_probs[k] for k in off_k.tolist()]).argmax(1)
            rules["pipeline"] = r09["threshold"]
        report["rules"] = metrics.decision_rules(labels, probs, rules, truth_args, preds)
        in_test = labels >= 0
        off = a["off_reference"][positions]
        if rescue:
            calls = np.array([off_probs[k].argmax() == metrics.SOMATIC for k in off_k.tolist()], bool)
            report["offref_rescue"] = dict(site_halfwidth=args.offref_site90, floor=PATH_COUNT_ALL[index.path_count],
                                           tensors=int(len(off_k)),
                                           in_test=int(in_test[off_k].sum()), calls=int(calls.sum()),
                                           calls_in_test=int((calls & in_test[off_k]).sum()),
                                           calls_before=int((probs[off_k].argmax(1) == metrics.SOMATIC).sum()))
        report.update(all_tensors=int(len(labels)), left_out=int((~in_test).sum()),
                      off_reference_in_test=int((in_test & off).sum()),
                      off_reference_somatic_in_test=int((in_test & off & (labels == 1)).sum()),
                      ignored_reasons=list(ignore), ignored=int(ignored.sum()),
                      directory=str(index.dir), labels=index.meta["labels"],
                      path_count=dict(storage=index.path_count, haplotypes=index.haplotypes,
                                      checkpoint_haplotypes=trained_haplotypes),
                      checkpoint=str(Path(args.checkpoint).resolve()), chroms=args.chroms)
        (out / f"{name}.metrics.json").write_text(json.dumps(report, indent=2) + "\n")
        print(f"{name}: {len(labels):,} tensors, {report['tensors']:,} in the test "
              f"({report['off_reference_in_test']:,} off-reference) | {describe(report)}", flush=True)
        if "offref_rescue" in report:
            o = report["offref_rescue"]
            print(f"  off-reference rescue (site +- {o['site_halfwidth']}, ch6 floor {o['floor']}): {o['tensors']:,} "
                  f"tensors, argmax somatic {o['calls_before']:,} before, {o['calls']:,} after "
                  f"({o['calls_in_test']:,} in the test)", flush=True)
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

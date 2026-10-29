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
  a tensor labelled 0 or 2. Reasons in --ignore-reasons are left out like -1; their records say ignored: true and
  keep label, reason and truth_ids. The default is the reasons the model was trained without, except that INDEL
  sets keep the partial somatic matches (residual_partial_somatic_truth, allele_partial_somatic_truth) in the test
  as 1: an INDEL tensor that partly matches a somatic truth (another length in the same homopolymer, or a residual
  edit next to a graph allele) counts that truth as found when called. SNV sets leave them out. The calls (graph_vcf)
  and so rtg vcfeval do not depend on this. Candidates in --exclude-candidates (default: the checkpoint's train
  --exclude-candidates) are left out of the test too, as if never built (their truth alleles are misses), and say
  excluded: true; graph_vcf writes no record for them. A found truth allele of the other kind (an SNV tensor matching
  an INDEL truth partially) is a true call too, once. `ceiling` is the recall of a perfect model.

Off-reference rescue (SNV, --offref-site90 W; off by default, W = 10 turns it on). The model learned that a low path count
(ch6, how many HPRC haplotypes run through the node) is not somatic, so it does not call somatic on an
off-reference node, where every read has a low path count. Each SNV off-reference tensor is scored again with its
ch6 cells below the every-haplotype code raised to it in the site columns (the middle column +- W; uncovered cells
stay 0): 90 on int8-count-linear100-log2 sets (HPRC v1.1, 90 haplotype paths), 100 on int8-count-haplotypes100
sets (data.PATH_COUNT_ALL). It is called somatic when somatic is the most probable class of those probabilities
(offref_call). graph_vcf uses offref_call for these tensors. `rules.pipeline` in .metrics.json scores that call on
the off-reference tensors and the checkpoint's validation recall-0.9 threshold on the others; `offref_rescue`
counts the tensors and calls and gives the floor.

Z-score statistics: with a checkpoint of norm "dataset" (train --norm dataset, the default from 2026-10-13), its
dataset_norm channels (BQ; BQ and MAPQ for the checkpoints of 2026-10-13 to 10-19) are z-scored with the statistics of
each predicted dataset: the checkpoint's own for one of its training datasets (the same directory), otherwise fitted
as train does, on --stats-samples tensors of all the dataset's chromosomes whatever --chroms is (labels not used;
recorded in .metrics.json "normalization"). MAPQ, path count, and every channel of a norm "train" checkpoint use the
checkpoint's training statistics; a dataset whose mean of such a channel lies more than 0.5 SD (of the nearest
training dataset) from every training dataset's gets a warning (.metrics.json "normalization" "shift": the shift to
the nearest training dataset in its SDs). --norm overrides the checkpoint's choice (--norm dataset on a norm "train"
checkpoint: encode.DATASET_NORM). If any set codes channel 6
otherwise than the checkpoint's training sets (checkpoint path_count_storage; int8-count-linear100-log2 when it
records none), predict exits before predicting anything (a model trained with --drop-planes path_count never reads
channel 6 and takes sets of any storage); a set built with another --haplotypes is predicted with a
note, and .metrics.json `path_count` gives the set's storage and H and the checkpoint's H.

Site features: a checkpoint trained with --site-features reads every predicted set's <DIR>/<SET>.<KIND>.npz
(site_features.py; DIR: the checkpoint's, or --site-features), with the transform the checkpoint records; a set without
its file, or with a stale one (another tensor order, another PoN rule), stops predict before anything is predicted.
Each record then gets pon_tag, and .metrics.json `rules` adds the rules after the PoN (tagged tensors never called):
pipeline_pon (the validation recall-0.9 threshold + PoN), threshold_t_pon (the threshold t + PoN),
val_postpon_recall (the post-PoN validation threshold at the checkpoint's call_recall, 0.8, of the callable positives:
the PASS rule of a checkpoint with call_rule postpon-recall, the INDEL recipe) and, as reports, val_postpon (train's
somatic_threshold_postpon, the post-PoN best-F1 threshold) and val_postpon_truth (truth_threshold_postpon), each when
the checkpoint has it. The checkpoint records the sha256 of
each training set's feature file (feature_files); when the file predict reads for one of those sets differs, predict
prints a warning to stderr and .metrics.json has it (feature_files: the file, its sha256, the checkpoint's and changed;
feature_file_warnings, which also name a file built with another version of the PoN script).
--non-reasons / the checkpoint's train --non-reasons and --non-reasons-in-bed score those tensors as non inside the
benchmark, as train did; the calls (graph_vcf) do not depend on it.

An unlabelled set (no labels.manifest.json, e.g. a new sample; data.py) is predicted the same way: every record has
test_label null and reason "unlabelled", and the test metrics are empty. graph_vcf and linear_vcf need no labels.
"""
import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import torch

from . import metrics
from . import site_features as sitef
from .data import (CHECKPOINT_PATH_COUNT, CLASSES, KINDS, PATH_COUNT_ALL, TENSOR_SHAPE, EpochSampler, TensorDataset,
                   load_parts, read_candidate_ids, somatic_truth)
from .encode import CONTINUOUS, DATASET_NORM
from .model import PansomaNetV2
from .train import dataset_stats, describe, in_somatic_bed, make_loader, predict_probs, retry_workers


def stat_shift(stats, references):
    """{channel: (mean - reference mean) / reference std} of every CONTINUOUS channel against the nearest of the
    reference statistics (the checkpoint's training datasets); None without references."""
    if not references:
        return None
    out = {}
    for name, _, _ in CONTINUOUS:
        d = [(stats[name]["mean"] - r[name]["mean"]) / max(r[name]["std"], 1e-6) for r in references]
        out[name] = round(min(d, key=abs), 4)
    return out


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
    p.add_argument("--offref-site90", type=int, default=0, metavar="W",
                   help="SNV off-reference rescue: score the off-reference tensors again with ch6 below the "
                        "every-haplotype code (90 on int8-count-linear100-log2 sets, 100 on int8-count-haplotypes100) "
                        "raised to it in the site columns +- W and call them by argmax (default 0: off; 10 turns it on)")
    p.add_argument("--exclude-candidates", nargs="?", const="",
                   help="file of candidate ids taken out before the model (train --exclude-candidates): left out of "
                        "the test like -1, records get excluded: true and graph_vcf leaves them out of the calls "
                        "(default: the checkpoint's train --exclude-candidates; the option without a file: none)")
    p.add_argument("--norm", choices=["dataset", "train"],
                   help="BQ z-score: dataset (each predicted dataset's statistics) or train (the checkpoint's training "
                        "statistics); default: the checkpoint's config norm")
    p.add_argument("--stats-samples", type=int, default=20000,
                   help="tensors of each predicted dataset (all its chromosomes) used to fit its statistics (norm "
                        "dataset, for a dataset the checkpoint was not trained on)")
    p.add_argument("--ignore-reasons", nargs="*",
                   help="label reasons left out of the test, as if -1 (default: the checkpoint's train --ignore-reasons, "
                        "minus the partial somatic reasons for INDEL sets); their records keep label, reason and "
                        "truth_ids and get ignored: true")
    p.add_argument("--site-features", metavar="DIR",
                   help="site-feature directory of a checkpoint trained with --site-features (default: the one the "
                        "checkpoint records); <DIR>/<SET>.<KIND>.npz must exist for every predicted set")
    p.add_argument("--non-reasons", nargs="*",
                   help="label reasons scored as non (0) inside the benchmark, the confident region or the somatic "
                        "truth BED (default: the checkpoint's train --non-reasons); the calls (graph_vcf) do not "
                        "depend on it")
    return p.parse_args(argv)


PATH_COUNT = dict((name, channel) for name, channel, _ in CONTINUOUS)["path_count"]
PARTIAL_SOMATIC = ("residual_partial_somatic_truth", "allele_partial_somatic_truth")
SITE = TENSOR_SHAPE[2] // 2


class OffrefRescue(TensorDataset):
    """Tensors with the ch6 (path count) cells below `floor` set to `floor` in the columns SITE +- halfwidth."""

    def __init__(self, parts, halfwidth, floor, labels="label"):
        super().__init__(parts, labels)
        self.halfwidth, self.floor = halfwidth, floor

    def __getitem__(self, i):
        x, blocks, scalars, label, part = super().__getitem__(i)
        window = x[PATH_COUNT, :, SITE - self.halfwidth:SITE + self.halfwidth + 1]
        window[(window > 0) & (window < self.floor)] = self.floor
        return x, blocks, scalars, label, part


def main(argv=None):
    args = parse_args(argv)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model = PansomaNetV2.from_checkpoint(checkpoint).to(device).eval()
    amp = args.amp == "bf16" and device.type == "cuda"
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    norm = args.norm or model.config["norm"]
    if norm == "dataset" and not model.encoder.dataset_norm:  # a norm train checkpoint with --norm dataset
        model.encoder.dataset_norm = DATASET_NORM
    channels = list(model.encoder.dataset_norm) if norm == "dataset" else []
    trained_stats = {d["directory"]: d["stats"] for d in checkpoint.get("dataset_stats", [])}
    print(f"checkpoint {args.checkpoint} (epoch {checkpoint['epoch']}), statistics {json.dumps(model.encoder.stats())}"
          f"; z-score: {norm}" + (f" (per dataset: {', '.join(channels)})" if channels else ""), flush=True)

    def probs_of(dataset, what):
        """(labels, probabilities) of every sample of dataset, in order."""
        return retry_workers(
            lambda: predict_probs(model, make_loader(dataset, EpochSampler(len(dataset), shuffle=False), args.num_workers,
                                                     args.batch_size, persistent=False), device, amp),
            lambda m: print(m, flush=True), what)

    trained = checkpoint.get("path_count_storage", CHECKPOINT_PATH_COUNT)
    trained_haplotypes = checkpoint.get("path_count_haplotypes") or []
    reads_path_count = "path_count" not in model.config["drop_planes"]  # the only plane made from channel 6
    parts = load_parts(args.tensors, args.kinds, args.cache_dir or str(out / "index_cache"), chroms=args.chroms,
                       labelled=False)
    for index, _ in parts:  # before any prediction: every set must code channel 6 like the training sets
        if reads_path_count and index.path_count != trained:
            raise SystemExit(f"{index.dir}: channel 6 is {index.path_count}, the checkpoint was trained on {trained}")
    n_site, recorded = model.config.get("site_features", 0), checkpoint.get("site_features")
    if n_site:
        if not recorded or recorded.get("version") != sitef.VERSION or recorded.get("names") != list(sitef.SITE_FEATURES):
            raise SystemExit(f"{args.checkpoint}: site features {recorded and recorded.get('version')} "
                             f"{recorded and recorded.get('names')}; this code reads {sitef.VERSION} {list(sitef.SITE_FEATURES)}")
        site_dir = args.site_features or recorded["directory"]
    elif args.site_features:
        raise SystemExit(f"{args.checkpoint} was trained without site features")
    tagged = bool(n_site)
    recorded_files = (checkpoint.get("feature_files") or {}).get("site_features")  # train's sha256 of each file
    feature_files = {}  # index directory -> (metrics.json feature_files, warnings)
    for index, _ in parts:  # every set's file before any prediction: a missing or stale one stops predict here
        entry, warnings = None, []
        if n_site:
            try:
                d = sitef.attach(index, site_dir, True)
            except ValueError as e:
                raise SystemExit(f"site features: {e}") from None
            print(sitef.describe(index, d) + " (model input)", flush=True)
            compared, warnings = sitef.compare_files({sitef.key(index): (d["path"], d["sha256"])}, recorded_files,
                                                     "site-features")
            entry = compared[sitef.key(index)]
            warnings += [f"warning: {w}" for w in sitef.pon_check(d)]
        for w in warnings:
            print(w, file=sys.stderr, flush=True)
        feature_files[str(index.dir)] = (entry, warnings)
    exclude_path = args.exclude_candidates if args.exclude_candidates is not None else \
        checkpoint.get("args", {}).get("exclude_candidates")
    exclude = read_candidate_ids(exclude_path) if exclude_path else set()
    if exclude:
        print(f"--exclude-candidates {exclude_path}: {len(exclude):,} candidate ids", flush=True)
    for index, positions in parts:
        if reads_path_count and index.haplotypes is not None and index.haplotypes not in trained_haplotypes:
            print(f"note: {index.dir} was built with --haplotypes {index.haplotypes}, the checkpoint's training sets "
                  f"with {trained_haplotypes or 'none recorded'}", flush=True)
        dataset = TensorDataset([(index, positions)], labels="eval_label")
        trained_shape = tuple(checkpoint.get("tensor_shape", TENSOR_SHAPE))  # checkpoints before it: (8, 200, 101)
        if dataset.shape != trained_shape:
            raise SystemExit(f"{index.dir}: tensors of shape {dataset.shape}, the checkpoint was trained on {trained_shape}")
        part_stats, source, shift = None, None, None
        if norm == "dataset":  # this dataset's statistics: row 0 of the table, the part of every sample here
            if str(index.dir) in trained_stats:  # a training dataset: the statistics the model was trained with
                part_stats, source = trained_stats[str(index.dir)], "checkpoint"
            else:  # fitted as train does: tensors of all the dataset's chromosomes, not only --chroms
                part_stats = dataset_stats(index, index.select(None, labelled=False), args.stats_samples, 12345,
                                           args.num_workers, lambda m: print(m, flush=True), f"statistics of {index.dir}")
                source = "fitted"
            shift = stat_shift(part_stats, list(trained_stats.values()))
            print(f"{index.dir}: statistics ({source}, {part_stats['tensors']:,} tensors) "
                  + json.dumps({k: part_stats[k] for k in channels})
                  + (f"; shift to the nearest training dataset (SD): {json.dumps(shift)}" if shift else ""), flush=True)
            for name, d in (shift or {}).items():
                if name not in channels and abs(d) > 0.5:
                    print(f"warning: {index.dir}: {name} mean {d:+.2f} SD from every training dataset's, and {name} "
                          f"uses the training statistics", flush=True)
        model.encoder.set_part_stats([part_stats] if part_stats else [])
        if args.ignore_reasons is not None:
            ignore = args.ignore_reasons
        else:  # the training's; INDEL tests count the partial somatic matches as found truth (see above)
            ignore = checkpoint.get("args", {}).get("ignore_reasons", [])
            if index.kind == "INDEL":
                ignore = [r for r in ignore if r not in PARTIAL_SOMATIC]
        ignored = np.isin(index.arrays["reason"][positions],
                          [k for k, r in enumerate(index.meta["reasons"]) if r in set(ignore)])
        dataset.label[ignored] = -1  # out of the test like the other -1; the record keeps its label and reason
        excluded = np.array([index.candidates()[p] in exclude for p in positions], bool) if exclude else \
            np.zeros(len(positions), bool)
        dataset.label[excluded] = -1  # never built, as far as the test and the calls are concerned
        non = args.non_reasons if args.non_reasons is not None else checkpoint.get("args", {}).get("non_reasons") or []
        as_non = (np.isin(index.arrays["reason"][positions],
                          [k for k, r in enumerate(index.meta["reasons"]) if r in set(non)])
                  & (index.arrays["in_region"] | in_somatic_bed(index))[positions] & ~ignored & ~excluded)
        dataset.label[as_non] = 0  # e.g. rule (a) INDELs: test-time calling meets them inside the BED
        bed_non = checkpoint.get("args", {}).get("non_reasons_in_bed") or []
        if bed_non:  # e.g. outside_confident_region inside the somatic BED (train --non-reasons-in-bed)
            in_bed = (np.isin(index.arrays["reason"][positions],
                              [k for k, r in enumerate(index.meta["reasons"]) if r in set(bed_non)])
                      & in_somatic_bed(index)[positions] & ~ignored & ~excluded)
            dataset.label[in_bed] = 0
            as_non = as_non | in_bed
        labels, probs = probs_of(dataset, f"{index.dir}")
        threshold = args.threshold if args.threshold is not None else checkpoint.get("somatic_threshold")
        pred = metrics.call(probs, threshold)
        name = f"{index.dir.parent.parent.name}.{index.dir.parent.name}.{index.kind}"
        a, chroms, candidates = index.arrays, index.meta["chroms"], index.candidates()
        rescue = index.kind == "SNV" and args.offref_site90 > 0
        off_k = np.flatnonzero(a["off_reference"][positions]) if rescue else np.array([], np.int64)
        off_probs = {}
        if len(off_k):
            edited = OffrefRescue([(index, positions[off_k])], args.offref_site90, PATH_COUNT_ALL[index.path_count])
            _, probs2 = probs_of(edited, f"{index.dir} off-reference rescue")
            off_probs = dict(zip(off_k.tolist(), probs2))
        with gzip.open(out / f"{name}.predictions.ndjson.gz", "wt") as f:
            for k, pos in enumerate(positions):
                f.write(json.dumps(dict(
                    chrom=chroms[a["chrom"][pos]], candidate_id=candidates[pos], label=int(a["label"][pos]),
                    test_label=int(labels[k]) if labels[k] >= 0 else None, in_test=bool(labels[k] >= 0),
                    reason=index.reason(pos), ignored=bool(ignored[k]), excluded=bool(excluded[k]),
                    off_reference=bool(a["off_reference"][pos]),
                    **({"pon_tag": int(a["pon_tag"][pos])} if tagged else {}),
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
        if tagged:  # after the PoN: tagged tensors are never called (their calls go to non / germline)
            pon = a["pon_tag"][positions] == 1
            val = checkpoint.get("val") or {}
            call_recall = checkpoint.get("call_recall") or 0.8  # graph_vcf's postpon-recall call rule
            ppr = ((val.get("somatic_at_recall_postpon_callable") or {}).get(f"{call_recall:g}") or {}).get("threshold")
            for rule, t in (("pipeline_pon", r09["threshold"] if r09 else None),
                            ("threshold_t_pon", threshold),
                            ("val_postpon_recall", ppr),
                            ("val_postpon", val.get("somatic_threshold_postpon")),
                            ("val_postpon_truth", val.get("truth_threshold_postpon"))):
                if t is None:
                    continue
                pred_pon = preds["pipeline"].copy() if rule == "pipeline_pon" else metrics.call(probs, t)
                pred_pon[pon] = np.where(probs[pon, 2] > probs[pon, 0], 2, 0)
                rules[rule], preds[rule] = t, pred_pon
            report["pon_tagged"] = int(pon.sum())
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
                      non_reasons=list(non), as_non=int(as_non.sum()),
                      exclude_candidates=exclude_path or None, excluded=int(excluded.sum()),
                      directory=str(index.dir), labels=index.meta["labels"],
                      path_count=dict(storage=index.path_count, haplotypes=index.haplotypes,
                                      checkpoint_haplotypes=trained_haplotypes),
                      checkpoint=str(Path(args.checkpoint).resolve()), chroms=args.chroms,
                      site_features=dict(directory=str(Path(site_dir).resolve()), version=sitef.VERSION) if n_site else None,
                      feature_files=feature_files[str(index.dir)][0],
                      feature_file_warnings=feature_files[str(index.dir)][1],
                      normalization=dict(norm=norm, channels=channels, source=source, dataset_stats=part_stats,
                                         shift=shift, checkpoint_stats=model.encoder.stats()))
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

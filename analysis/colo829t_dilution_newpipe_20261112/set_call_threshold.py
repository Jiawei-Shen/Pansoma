#!/usr/bin/env python3
"""Move the PASS threshold a checkpoint records by OFFSET (the user, 2026-11-12: the Illumina and fiberseq dilution
models call at their validation threshold - 0.10). graph_vcf's default PASS threshold of a --best-by postpon checkpoint
is val.somatic_at_recall_postpon_callable[call_recall].threshold (call_rule postpon-recall); that record becomes
{threshold: validation threshold + OFFSET (6 decimals), manual: true, offset, note, validation: the record train wrote},
and the checkpoint gets call_threshold_manual with the same facts. Weights, the other validation records, call_rule and
call_recall are unchanged, so graph_vcf still names the rule "validation post-PoN recall 0.95 of the callable
positives" in the VCF header: the threshold it prints there is the manual one (the validation recall at it was not
measured). The file is replaced only when --original (an unchanged copy that stays) is byte-identical to it; a
checkpoint already changed is left as it is.
Usage (in the container, through set_call_threshold.sh): set_call_threshold.py CHECKPOINT OFFSET --original COPY"""
import argparse
import filecmp
import time
from pathlib import Path

import torch

p = argparse.ArgumentParser()
p.add_argument("checkpoint", type=Path)
p.add_argument("offset", type=float)
p.add_argument("--original", type=Path, required=True)
p.add_argument("--note", default="")
args = p.parse_args()

ck = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
rule, recall = ck.get("call_rule"), ck.get("call_recall")
if rule != "postpon-recall" or recall is None:
    raise SystemExit(f"{args.checkpoint}: call rule {rule} at {recall}, not postpon-recall")
at = ck["val"]["somatic_at_recall_postpon_callable"]
record = at[f"{recall:g}"]
if record.get("manual"):
    raise SystemExit(f"{args.checkpoint}: already changed: {record}")
if not filecmp.cmp(args.checkpoint, args.original, shallow=False):
    raise SystemExit(f"{args.original} is not a copy of {args.checkpoint}: the unchanged checkpoint would be lost")
threshold = round(record["threshold"] + args.offset, 6)
if not 0 < threshold < 1:
    raise SystemExit(f"threshold {threshold} out of range")
note = args.note or f"validation threshold {args.offset:+g}"
at[f"{recall:g}"] = dict(threshold=threshold, manual=True, offset=args.offset, note=note, validation=record)
ck["call_threshold_manual"] = dict(threshold=threshold, offset=args.offset, validation_threshold=record["threshold"],
                                   call_rule=rule, call_recall=recall, note=note, original=str(args.original),
                                   date=time.strftime("%Y-%m-%d %H:%M:%S"))
tmp = args.checkpoint.with_name(args.checkpoint.name + ".tmp")
torch.save(ck, tmp)
tmp.replace(args.checkpoint)
print(f"{args.checkpoint}: epoch {ck['epoch']}, {rule} at {recall:g}: PASS threshold {record['threshold']} -> "
      f"{threshold} ({note}); unchanged copy {args.original}")

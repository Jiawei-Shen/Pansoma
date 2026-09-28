"""The SNV and INDEL predictions of a sample together against its somatic truth VCF, as a pipeline's output VCF
(the calls of both models) is compared with it: every truth allele counts once.

    cd machine_learning
    python -m pansoma_net_v2.combine HG008_Illumina_SNV_base HG008_Illumina_INDEL_base [--output combined.json]

Arguments: prediction files (<sample>.<set>.<KIND>.predictions.ndjson.gz, from predict), directories holding
them, or run names / directories (their test_chr1/). Files are grouped by tensor set, one per kind; each group is
scored against the set's somatic.recall.tsv (from the directory in the file's .metrics.json): the truth alleles
that pass and lie in the BED, of the kinds predicted (SNV: SNP; INDEL: DEL, INS), on the predicted chromosomes.

- A truth allele is found when any test tensor standing for it (labelled 1: an exact or partial match, in either
  set) is called somatic by its model (pred, with that model's threshold). It counts once, however many of its
  tensors are called and in whichever set: an SNV tensor that partially matches an INDEL truth finds that truth.
- A found truth allele of a kind not predicted (only one set given) is a true call, not in the recall.
- A false call is a somatic call on a test tensor scored 0 or 2.
- precision = found / (found + false calls); recall = found / truth alleles. No AP: the two models' scores are
  not on one scale.
"""
import argparse
import csv
import glob
import gzip
import json
from pathlib import Path

from .data import TRUTH_KINDS
from .summary import run_dir

SUFFIX = ".predictions.ndjson.gz"


def prediction_files(paths):
    out = []
    for p in paths:
        p = Path(p)
        if not p.exists():
            p = run_dir(str(p))
        if p.is_dir():
            files = sorted(glob.glob(str(p / f"*{SUFFIX}"))) or sorted(glob.glob(str(p / "test_chr1" / f"*{SUFFIX}")))
            if not files:
                raise SystemExit(f"no *{SUFFIX} in {p} or {p / 'test_chr1'}")
            out += [Path(f) for f in files]
        else:
            out.append(p)
    return out


def groups(files):
    """{tensor set directory: {kind: prediction file}}, the set and kind from each file's .metrics.json."""
    out = {}
    for f in files:
        meta = json.loads(Path(str(f)[:-len(SUFFIX)] + ".metrics.json").read_text())
        kind_dir = Path(meta["directory"])
        by_kind = out.setdefault(kind_dir.parent, {})
        if kind_dir.name in by_kind:
            raise SystemExit(f"two {kind_dir.name} predictions of {kind_dir.parent}: {by_kind[kind_dir.name]} and {f}")
        by_kind[kind_dir.name] = f
    return out


def score(set_dir, by_kind):
    records = {}
    for kind, f in by_kind.items():
        with gzip.open(f, "rt") as fh:
            records[kind] = [r for r in map(json.loads, fh) if r["in_test"]]
    chroms = {r["chrom"] for rs in records.values() for r in rs}
    kinds = {k for kind in by_kind for k in TRUTH_KINDS[kind]}
    with open(set_dir / "somatic.recall.tsv") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["passed"] == "True" and r["in_bed"] == "True" and r["chrom"] in chroms]
    truth = {int(r["truth_id"]): r["kind"] for r in rows if r["kind"] in kinds}
    other = {int(r["truth_id"]): r["kind"] for r in rows if r["kind"] not in kinds}
    reached, found = {}, {}  # truth id -> the sets with a test tensor of it / with a somatic call on it
    calls = {}
    for kind, rs in records.items():
        c = dict(calls=0, on_truth=0, repeated=0, false=0, on_uncounted_truth=0)
        for r in rs:
            ids = [t for t in r["truth_ids"] if t in truth or t in other] if r["test_label"] == 1 else []
            for t in ids:
                reached.setdefault(t, set()).add(kind)
            if r["pred"] != "somatic":
                continue
            c["calls"] += 1
            if r["test_label"] == 1:
                c["on_truth" if ids else "on_uncounted_truth"] += 1
                c["repeated"] += bool(ids) and all(t in found for t in ids)  # its truth found by an earlier call
                for t in ids:
                    found.setdefault(t, set()).add(kind)
            else:
                c["false"] += 1
        calls[kind] = c
    tp = sum(t in truth for t in found)
    other_tp = sum(t in other for t in found)
    fp = sum(c["false"] for c in calls.values())
    prec = (tp + other_tp) / (tp + other_tp + fp) if tp + other_tp + fp else 0.0
    rec = tp / len(truth) if truth else 0.0
    per_kind = {}
    for k in sorted(set(truth.values())):
        ids = [t for t, v in truth.items() if v == k]
        per_kind[k] = dict(truth_alleles=len(ids), with_tensor=sum(t in reached for t in ids),
                           found=sum(t in found for t in ids))
        per_kind[k]["recall"] = per_kind[k]["found"] / len(ids)
    found_by = {}
    for t, sets in found.items():
        key = "+".join(sorted(sets))
        found_by[key] = found_by.get(key, 0) + 1
    return dict(tensor_set=str(set_dir), predictions={k: str(f) for k, f in by_kind.items()}, chroms=sorted(chroms),
                truth_alleles=len(truth), with_tensor=sum(t in reached for t in truth),
                ceiling=sum(t in reached for t in truth) / len(truth) if truth else 0.0,
                tp=tp, other_kind_tp=other_tp, fp=fp, fn=len(truth) - tp, precision=prec, recall=rec,
                f1=2 * prec * rec / (prec + rec) if prec + rec else 0.0,
                per_truth_kind=per_kind, found_by=found_by, calls=calls,
                repeated_calls=sum(c["repeated"] for c in calls.values()))


def show(r):
    print(f"== {r['tensor_set']}  ({', '.join(f'{k}: {Path(f).parent}' for k, f in r['predictions'].items())}; "
          f"{','.join(r['chroms'])})")
    print(f"  F1 {r['f1']:.3f}  P {r['precision']:.3f}  R {r['recall']:.3f}  | {r['tp']:,} of {r['truth_alleles']:,} "
          f"truth found, {r['fp']:,} false calls, ceiling {r['ceiling']:.3f}"
          + (f", {r['other_kind_tp']:,} found of kinds not predicted" if r["other_kind_tp"] else ""))
    for k, v in r["per_truth_kind"].items():
        print(f"  {k:>5}: {v['found']:,} of {v['truth_alleles']:,} found (R {v['recall']:.3f}), "
              f"{v['with_tensor']:,} with a test tensor")
    print(f"  found by: {', '.join(f'{k} {n:,}' for k, n in sorted(r['found_by'].items()))}; "
          f"{r['repeated_calls']:,} more calls on truth already found (not counted again)")
    for kind, c in r["calls"].items():
        print(f"  {kind} calls: {c['calls']:,} ({c['on_truth']:,} on truth, {c['false']:,} false"
              + (f", {c['on_uncounted_truth']:,} on truth outside the counted set" if c["on_uncounted_truth"] else "")
              + ")")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("predictions", nargs="+", help="prediction files, their directories, or run names / directories")
    p.add_argument("--output", help="write the reports as JSON")
    args = p.parse_args(argv)
    reports = [score(d, by_kind) for d, by_kind in groups(prediction_files(args.predictions)).items()]
    for r in reports:
        show(r)
    if args.output:
        Path(args.output).write_text(json.dumps(reports, indent=2) + "\n")
    return reports


if __name__ == "__main__":
    main()

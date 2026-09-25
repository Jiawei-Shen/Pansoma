#!/usr/bin/env python3
"""Command line for the indexed GAM tensor pipeline.

    discover  scan the GAM once and select nodes with an imperfect-mapping fraction above a threshold
    build     build SNV and INDEL site tensors for a node list (see build.py for the outputs)
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path

from .common import new_output, write_json
from .gam_reader import scan_gam


def discover(args):
    """Select nodes where more than `node_alt` of MAPQ-passing mappings carry an edit."""
    out = new_output(args.output)
    stats = defaultdict(lambda: [0, 0, 0])  # perfect, imperfect, max read length
    count = 0
    for alignment in scan_gam(args.gam, args.max_alignments):
        count += 1
        if alignment.mapping_quality <= args.min_mapq:
            continue
        for mapping in alignment.path.mapping:
            nid = mapping.position.node_id
            if not nid:
                continue
            imperfect = any(e.from_length != e.to_length or e.sequence for e in mapping.edit)
            stats[nid][int(imperfect)] += 1
            stats[nid][2] = max(stats[nid][2], len(alignment.sequence))
        if count % 100000 == 0:
            print(f"Scanned {count:,} alignments; {len(stats):,} nodes", flush=True)
    selected = sorted(nid for nid, (p, n, _) in stats.items() if n >= 1 and n / (p + n) > args.node_alt)
    (out / "target_nodes.txt").write_text("".join(f"{n}\n" for n in selected))
    write_json(out / "node_stats.json", {str(n): dict(perfect=p, not_perfect=q, max_read_length=r)
                                         for n, (p, q, r) in stats.items()})
    report = dict(gam=str(Path(args.gam).resolve()), alignments_scanned=count, nodes_observed=len(stats),
                  nodes_selected=len(selected), min_mapq=args.min_mapq, node_alt=args.node_alt,
                  exploratory=bool(args.max_alignments), max_alignments=args.max_alignments)
    write_json(out / "discovery_report.json", report)
    print(json.dumps(report, indent=2))


def build(args):
    from .build import build as build_tensors
    return build_tensors(args)


def positive(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def nonnegative(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return number


def fraction(value):
    number = float(value)
    if not 0 <= number <= 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return number


def add_build_arguments(sub, outputs=True):
    """Builder options shared with the orchestrator; `outputs=False` omits per-run output options."""
    sub.add_argument("--graph-index", required=True, help="unified GBZ sequence/path-count SQLite (graph_index.py build)")
    sub.add_argument("--snv-min-af", type=fraction, required=True, help="AF threshold for the SNV output")
    sub.add_argument("--indel-min-af", type=fraction, required=True, help="AF threshold for the INDEL output")
    if outputs:
        sub.add_argument("--snv-output", required=True, help="separate SNV output directory; one shared decoding pass")
        sub.add_argument("--indel-output", required=True, help="separate INDEL output directory; one shared decoding pass")
        sub.add_argument("--debug-rows", action="store_true", help="record per-row source hashes and per-column graph coordinates (needed by validate_examples.py)")
    sub.add_argument("--rows", type=positive, default=200)
    sub.add_argument("--width", type=positive, default=101)
    sub.add_argument("--gam-cache-mb", type=positive, default=1024, help="bounded GAM group cache in MiB (at least 1)")
    sub.add_argument("--batch-nodes", type=positive, default=512, help="target nodes per batch")
    sub.add_argument("--max-node-span", type=positive, default=10000, help="maximum node-ID span of one batch")
    sub.add_argument("--max-batch-alignments", type=positive, default=20000, help="fail instead of decoding a larger batch")
    sub.add_argument("--shard-size", type=positive, default=2048, help="tensors per NPY shard")
    sub.add_argument("--min-mapq", type=int, default=10, help="exclusive: alignments with MAPQ <= this are dropped")
    sub.add_argument("--min-af", type=fraction, default=0.05, help="AF threshold (single-output mode)")
    sub.add_argument("--min-variants", type=positive, default=3, help="minimum ALT-supporting records")
    sub.add_argument("--min-allele-bq", type=float, default=10, help="minimum base quality for an ALT observation")
    sub.add_argument("--max-indel-len", type=positive, default=50)
    sub.add_argument("--max-node-reads", type=nonnegative, default=800,
                     help="per target node, count support and select rows from at most N records (smallest "
                          "record SHA-256), after the prefilters; 0 = no cap")
    sub.add_argument("--chromosomes", default="all",
                     help="target nodes to keep, by the chromosome block of their node ID: all (default), autosome "
                          "(chr1-22) or a comma list of --chr-index block names (e.g. chr1,chr2,chrX); applied to "
                          "the node list before any batch")
    sub.add_argument("--chr-index", help="node ID -> chromosome block table (tensor_postprocessing chr-index .tsv); "
                                         "required unless --chromosomes all")
    sub.add_argument("--early-af-filter", action=argparse.BooleanOptionalAction, default=True,
                     help="before support counting, reject candidates whose ALT support bound / exact coverage "
                          "over all records is below the AF threshold (default on; exact without a read cap)")
    sub.add_argument("--decoder", choices=("auto", "native", "python"), default="auto",
                     help="record decoder: auto (default) = the native C++ decoder if built and self-tested (native.py "
                          "compile), else Python; identical outputs. PANSOMA_DECODER=native|python overrides auto")


def make_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("discover", "build"):
        sub = commands.add_parser(name)
        sub.add_argument("--gam", required=True, help="sorted BGZF GAM")
        sub.add_argument("--output", required=True, help="new or empty directory")
        if name == "discover":
            sub.add_argument("--min-mapq", type=int, default=5, help="exclusive MAPQ threshold")
            sub.add_argument("--node-alt", type=fraction, default=0.05, help="select nodes with imperfect fraction > this")
            sub.add_argument("--max-alignments", type=positive, help="exploratory partial scan only")
        else:
            sub.add_argument("--nodes", required=True, help="target node IDs, one per line")
            sub.add_argument("--index", help="default: GAM path + .gai")
            add_build_arguments(sub)
    return parser


def main(argv=None):
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        {"discover": discover, "build": build}[args.command](args)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()

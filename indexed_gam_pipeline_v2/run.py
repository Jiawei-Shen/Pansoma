#!/usr/bin/env python3
"""Command line for the indexed GAM tensor pipeline.

    index     build a .gai for a sorted BGZF GAM (only if vg did not produce one)
    discover  scan the GAM once and select nodes with an imperfect-mapping fraction above a threshold
    validate  check indexed retrieval for a node set against an independent full scan
    build     build candidate tensors for a node list (see build.py for the outputs)
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from indexed_gam_pipeline_v2.common import load_nodes, new_output, write_json
from indexed_gam_pipeline_v2.gam_reader import IndexedGam, build_index, scan_gam


def index_gam(args):
    output = Path(args.output)
    if output.exists():
        raise ValueError(f"Index output already exists: {output}")
    report = build_index(args.gam, output)
    print(json.dumps(dict(gam=str(Path(args.gam).resolve()), index=str(output.resolve()), **report), indent=2))


def discover(args):
    """Select nodes where more than `node_alt` of MAPQ-passing mappings carry an edit."""
    out = new_output(args.output)
    allowed = set(load_nodes(args.chr_nodes)) if args.chr_nodes else None
    stats = defaultdict(lambda: [0, 0, 0])  # perfect, imperfect, max read length
    count = 0
    for alignment in scan_gam(args.gam, args.max_alignments):
        count += 1
        if alignment.mapping_quality <= args.min_mapq:
            continue
        for mapping in alignment.path.mapping:
            nid = mapping.position.node_id
            if not nid or (allowed is not None and nid not in allowed):
                continue
            imperfect = any(e.from_length != e.to_length or e.sequence for e in mapping.edit)
            stats[nid][int(imperfect)] += 1
            stats[nid][2] = max(stats[nid][2], len(alignment.sequence))
        if count % 100000 == 0:
            print(f"Scanned {count:,} alignments; {len(stats):,} nodes", flush=True)
    selected = sorted(nid for nid, (p, n, _) in stats.items() if n >= 1 and n / (p + n) > args.node_alt)
    if args.max_nodes:
        selected = selected[:args.max_nodes]
    (out / "target_nodes.txt").write_text("".join(f"{n}\n" for n in selected))
    write_json(out / "node_stats.json", {str(n): dict(perfect=p, not_perfect=q, max_read_length=r)
                                         for n, (p, q, r) in stats.items()})
    report = dict(gam=str(Path(args.gam).resolve()), alignments_scanned=count, nodes_observed=len(stats),
                  nodes_selected=len(selected), min_mapq=args.min_mapq, node_alt=args.node_alt,
                  exploratory=bool(args.max_alignments or args.max_nodes),
                  max_alignments=args.max_alignments, max_nodes=args.max_nodes)
    write_json(out / "discovery_report.json", report)
    print(json.dumps(report, indent=2))


def validate(args):
    """Indexed fetch must return exactly the records a full sequential scan finds (as a multiset)."""
    out = new_output(args.output)
    nodes = set(load_nodes(args.nodes))
    reader = IndexedGam(args.gam, args.index)

    def digest(alignment):
        return hashlib.sha256(alignment.SerializeToString(deterministic=True)).hexdigest()

    start = time.monotonic()
    metrics = {}
    indexed = Counter(digest(a) for a in reader.fetch(nodes, metrics))
    indexed_seconds = time.monotonic() - start
    print(f"Indexed query: {metrics}", flush=True)
    start = time.monotonic()
    expected = Counter()
    scanned = 0
    for alignment in scan_gam(args.gam):
        scanned += 1
        if any(m.position.node_id in nodes for m in alignment.path.mapping):
            expected[digest(alignment)] += 1
        if scanned % 500000 == 0:
            print(f"Reference scan: {scanned:,} alignments", flush=True)
    ok = indexed == expected
    report = dict(passed=ok, gam=str(Path(args.gam).resolve()), gai_version=reader.version,
                  target_nodes=len(nodes), sequential_alignments=scanned,
                  indexed_query=metrics, expected_alignments=sum(expected.values()),
                  missing_alignments=sum((expected - indexed).values()),
                  extra_alignments=sum((indexed - expected).values()),
                  indexed_seconds=round(indexed_seconds, 3),
                  sequential_seconds=round(time.monotonic() - start, 3))
    write_json(out / "validation_report.json", report)
    print(json.dumps(report, indent=2))
    if not ok:
        raise ValueError("Indexed GAM retrieval differs from independent full scan")


def build(args):
    from indexed_gam_pipeline_v2.build import build as build_tensors
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
    sub.add_argument("--snv-min-af", type=fraction, help="AF threshold for the SNV output (split mode)")
    sub.add_argument("--indel-min-af", type=fraction, help="AF threshold for the INDEL output (split mode)")
    if outputs:
        sub.add_argument("--snv-output", help="separate SNV output directory; one shared decoding pass")
        sub.add_argument("--indel-output", help="separate INDEL output directory; one shared decoding pass")
        sub.add_argument("--debug-rows", action="store_true", help="record per-row source hashes and per-column graph coordinates (needed by validate_examples.py)")
        sub.add_argument("--max-tensors", type=positive, help="stop after this many tensors (bounded test runs)")
        sub.add_argument("--variant-type", choices=("snp", "indel", "all"), default="all")
    sub.add_argument("--rows", type=positive, default=200)
    sub.add_argument("--width", type=positive, default=101)
    sub.add_argument("--gam-cache-mb", type=int, default=1024, help="bounded GAM group cache in MiB; 0 disables reuse")
    sub.add_argument("--batch-nodes", type=positive, default=512, help="target nodes per batch")
    sub.add_argument("--max-node-span", type=positive, default=10000, help="maximum node-ID span of one batch")
    sub.add_argument("--max-batch-alignments", type=positive, default=20000, help="fail instead of decoding a larger batch")
    sub.add_argument("--shard-size", type=positive, default=2048, help="tensors per NPY shard")
    sub.add_argument("--min-mapq", type=int, default=10, help="exclusive: alignments with MAPQ <= this are dropped")
    sub.add_argument("--min-af", type=fraction, default=0.05, help="AF threshold (single-output mode)")
    sub.add_argument("--min-variants", type=positive, default=3, help="minimum ALT-supporting records")
    sub.add_argument("--min-allele-bq", type=float, default=10, help="minimum base quality for an ALT observation")
    sub.add_argument("--max-indel-len", type=positive, default=50)
    sub.add_argument("--candidate-unit", choices=("site", "allele"), default="site",
                     help="site: one tensor per (node, start, SNV|INDEL), every passing allele listed in the "
                          "summary; allele: one tensor per allele")
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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("index", "discover", "validate", "build"):
        sub = commands.add_parser(name)
        sub.add_argument("--gam", required=True, help="sorted BGZF GAM")
        sub.add_argument("--output", required=True, help="new .gai file for index; new or empty directory otherwise")
        if name == "discover":
            sub.add_argument("--min-mapq", type=int, default=5, help="exclusive MAPQ threshold")
            sub.add_argument("--node-alt", type=fraction, default=0.05, help="select nodes with imperfect fraction > this")
            sub.add_argument("--chr-nodes", help="optional node whitelist text file")
            sub.add_argument("--max-alignments", type=positive, help="exploratory partial scan only")
            sub.add_argument("--max-nodes", type=positive, help="limit selected nodes for smoke tests")
        elif name != "index":
            sub.add_argument("--nodes", required=True, help="target node IDs, one per line")
            sub.add_argument("--index", help="default: GAM path + .gai")
        if name == "build":
            add_build_arguments(sub)
    args = parser.parse_args(argv)
    try:
        {"index": index_gam, "discover": discover, "validate": validate, "build": build}[args.command](args)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()

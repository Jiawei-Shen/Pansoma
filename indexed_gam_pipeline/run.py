#!/usr/bin/env python3
"""Discover nodes, validate indexed reads, or build tensors directly from GAM."""
import argparse
from collections import Counter, defaultdict
from contextlib import closing, nullcontext
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from indexed_gam_pipeline.gam_reader import IndexedGam, build_index, scan_gam
from indexed_gam_pipeline.segments import on_chromosome, raw_segments, orient


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + "\n")


def new_output(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"Output directory must be empty: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_nodes(path):
    result = sorted({int(line.strip()) for line in Path(path).read_text().splitlines()
                     if line.strip()})
    if not result or result[0] <= 0:
        raise ValueError("Node list must contain positive node IDs")
    return result


def index_gam(args):
    output = Path(args.output)
    if output.exists():
        raise ValueError(f"Index output already exists: {output}")
    report = build_index(args.gam, output)
    print(json.dumps(dict(gam=str(Path(args.gam).resolve()), index=str(output.resolve()),
                          **report), indent=2))


def discover(args):
    out = new_output(args.output)
    allowed = set(load_nodes(args.chr_nodes)) if args.chr_nodes else None
    stats = defaultdict(lambda: [0, 0, 0, 0])
    count = 0
    for alignment in scan_gam(args.gam, args.max_alignments):
        count += 1
        if alignment.mapping_quality <= getattr(args, "min_mapq", 5) or not on_chromosome(alignment, getattr(args, "chr", "")):
            continue
        for mapping in alignment.path.mapping:
            nid = mapping.position.node_id
            if not nid or (allowed is not None and nid not in allowed):
                continue
            imperfect = any(e.from_length != e.to_length or e.sequence
                            for e in mapping.edit)
            stats[nid][int(imperfect)] += 1
            stats[nid][2] = max(stats[nid][2], sum(e.to_length for e in mapping.edit))
            cigar_length = sum(len(str(e.from_length if e.from_length == e.to_length
                                       or e.to_length == 0 else e.to_length)) + 1
                               for e in mapping.edit
                               if e.from_length == e.to_length or not e.from_length or not e.to_length)
            stats[nid][3] = max(stats[nid][3], cigar_length)
        if count % 100000 == 0:
            print(f"Scanned {count:,} alignments; {len(stats):,} nodes", flush=True)
    selected = sorted(nid for nid, (p, n, _, _) in stats.items()
                      if n >= 1 and n / (p + n) > args.node_alt)
    if args.max_nodes:
        selected = selected[:args.max_nodes]
    (out / "target_nodes.txt").write_text("".join(f"{n}\n" for n in selected))
    write_json(out / "node_stats.json", {
        str(n): {"perfect": p, "not_perfect": q,
                 "max_read_length": r, "max_cigar_length": c}
        for n, (p, q, r, c) in stats.items()})
    report = dict(gam=str(Path(args.gam).resolve()), alignments_scanned=count,
                  nodes_observed=len(stats), nodes_selected=len(selected),
                  node_alt=args.node_alt,
                  chromosome=getattr(args, "chr", ""),
                  exploratory=bool(args.max_alignments or args.max_nodes),
                  max_alignments=args.max_alignments, max_nodes=args.max_nodes)
    write_json(out / "discovery_report.json", report)
    print(json.dumps(report, indent=2))


def digest(alignment):
    # A multiset, not read-name deduplication: retain duplicate records and mates.
    return hashlib.sha256(alignment.SerializeToString(deterministic=True)).hexdigest()


def segment_digest(segment):
    return hashlib.sha256(json.dumps(segment, sort_keys=True).encode()).hexdigest()


def validate(args):
    out = new_output(args.output)
    nodes = set(load_nodes(args.nodes))
    reader = IndexedGam(args.gam, args.index)
    start = time.monotonic()
    metrics = {}
    indexed = Counter()
    indexed_segments = defaultdict(Counter)
    for alignment in reader.fetch(nodes, metrics):
        if not on_chromosome(alignment, getattr(args, "chr", "")):
            continue
        indexed[digest(alignment)] += 1
        for nid, segment in raw_segments(alignment, nodes):
            indexed_segments[nid][segment_digest(segment)] += 1
    indexed_seconds = time.monotonic() - start
    print(f"Indexed query: {metrics}", flush=True)
    expected = Counter()
    expected_segments = defaultdict(Counter)
    scanned = 0
    start = time.monotonic()
    for alignment in scan_gam(args.gam):
        scanned += 1
        if on_chromosome(alignment, getattr(args, "chr", "")) and any(
                m.position.node_id in nodes for m in alignment.path.mapping):
            expected[digest(alignment)] += 1
            for nid, segment in raw_segments(alignment, nodes):
                expected_segments[nid][segment_digest(segment)] += 1
        if scanned % 500000 == 0:
            print(f"Reference scan: {scanned:,} alignments", flush=True)
    ok = indexed == expected and indexed_segments == expected_segments
    report = dict(passed=ok, gam=str(Path(args.gam).resolve()),
                  gai_version=reader.version, target_nodes=len(nodes),
                  chromosome=getattr(args, "chr", ""),
                  sequential_alignments=scanned,
                  indexed_query=metrics,
                  expected_alignments=sum(expected.values()),
                  missing_alignments=sum((expected - indexed).values()),
                  extra_alignments=sum((indexed - expected).values()),
                  segments_equal=indexed_segments == expected_segments,
                  per_node_segments={str(n): sum(indexed_segments[n].values())
                                     for n in sorted(nodes)},
                  indexed_seconds=round(indexed_seconds, 3),
                  sequential_seconds=round(time.monotonic() - start, 3),
                  tensor_tested=False)
    write_json(out / "validation_report.json", report)
    print(json.dumps(report, indent=2))
    if not ok:
        raise ValueError("Indexed GAM retrieval differs from independent full scan")


def node_records(args, nodes, sqlite_connection=None):
    records = {}
    if args.node_json:
        opener = gzip.open if args.node_json.endswith(".gz") else open
        with opener(args.node_json, "rt") as stream:
            data = json.load(stream)
        if isinstance(data, dict):
            data = data["nodes"]
        for record in data:
            nid = int(record["node_id"])
            if nid in nodes:
                records[nid] = dict(record)
                if records[nid].get("sequence"):
                    records[nid]["sequence"] = records[nid]["sequence"].upper()
    if getattr(args, "node_sqlite", None):
        uri = Path(args.node_sqlite).resolve().as_uri() + "?mode=ro"
        manager = closing(sqlite3.connect(uri, uri=True)) if sqlite_connection is None else nullcontext(sqlite_connection)
        with manager as connection:
            # One read transaction avoids a filesystem lock cycle per node on
            # shared storage; bounded IN queries also amortize SQL execution.
            if sqlite_connection is None:
                connection.execute("BEGIN")
            ordered = sorted(nodes)
            for offset in range(0, len(ordered), 900):
                batch = ordered[offset:offset + 900]
                placeholders = ",".join("?" for _ in batch)
                for node_id, seq in connection.execute(
                        f"SELECT node_id, seq FROM nodes WHERE node_id IN ({placeholders})",
                        [str(n) for n in batch]):
                    if seq and seq != "*":
                        nid = int(node_id)
                        sequence = seq.upper()
                        rec = records.setdefault(nid, {"node_id": str(nid)})
                        if rec.get("sequence") and rec["sequence"] != sequence:
                            raise ValueError(f"SQLite/JSON sequence mismatch at node {nid}")
                        rec["sequence"] = sequence
    if args.gfa:
        opener = gzip.open if args.gfa.endswith(".gz") else open
        with opener(args.gfa, "rt") as stream:
            for line in stream:
                if line.startswith("S\t"):
                    fields = line.rstrip().split("\t")
                    nid = int(fields[1])
                    if nid in nodes and fields[2] != "*":
                        rec = records.setdefault(nid, {"node_id": str(nid)})
                        if rec.get("sequence") and rec["sequence"].upper() != fields[2].upper():
                            raise ValueError(f"GFA/JSON sequence mismatch at node {nid}")
                        rec["sequence"] = fields[2].upper()
    missing = [n for n in nodes if not records.get(n, {}).get("sequence")
               or records[n]["sequence"] == "*"]
    if missing:
        raise ValueError(f"Missing graph sequences for {len(missing)} nodes, e.g. {missing[:10]}")
    return records


def batches(nodes, size, max_span):
    batch = []
    for nid in sorted(nodes):
        if batch and (len(batch) >= size or nid - batch[0] > max_span):
            yield batch
            batch = []
        batch.append(nid)
    if batch:
        yield batch


def build_legacy(args):
    import numpy as np
    from pangenome_ml_data_generation.tensors.builders import build_tensors_from_segments

    nodes = load_nodes(args.nodes)
    records = node_records(args, set(nodes))
    reader = IndexedGam(args.gam, args.index)
    out = new_output(args.output)
    write_json(out / "candidate_nodes.json", [records[n] for n in nodes])
    (out / "target_nodes.txt").write_text("".join(f"{n}\n" for n in nodes))
    tensors, metadata = [], []
    shard = total = 0
    manifest = dict(status="running", arguments=vars(args), nodes=len(nodes),
                    gai_version=reader.version, tensor_format_version="indexed-gam-legacy-v1",
                    shape=[5, 201, 100], dtype="int8")
    write_json(out / "run_report.json", manifest)
    with (out / "variant_summary.ndjson").open("w") as summary:
        def flush():
            nonlocal shard, total
            if not tensors:
                return
            np.save(out / f"shard_{shard:05d}_data.npy", np.stack(tensors))
            for i, meta in enumerate(metadata):
                summary.write(json.dumps(dict(meta, shard_index=shard,
                                              index_within_shard=i)) + "\n")
            total += len(tensors)
            shard += 1
            tensors.clear()
            metadata.clear()

        for i, batch in enumerate(batches(nodes, args.batch_nodes, args.max_node_span), 1):
            wanted = set(batch)
            segments = defaultdict(list)
            metrics = {}
            segment_count = 0
            for alignment in reader.fetch(wanted, metrics):
                if not on_chromosome(alignment, getattr(args, "chr", "")):
                    continue
                for nid, segment in raw_segments(alignment, wanted, args.min_mapq):
                    segments[nid].append(orient(segment, len(records[nid]["sequence"])))
                    segment_count += 1
                    if segment_count > args.max_batch_segments:
                        raise ValueError("Batch segment limit exceeded; reduce --batch-nodes/--max-node-span "
                                         "or explicitly increase --max-batch-segments")
            for nid in batch:
                _, _, xs, metas = build_tensors_from_segments(
                    nid, records[nid]["sequence"], segments.pop(nid, []),
                    args.min_af, args.min_variants, args.min_allele_bq,
                    args.variant_type, args.max_indel_len)
                for tensor, meta in zip(xs, metas):
                    tensors.append(tensor)
                    metadata.append(meta)
                    if len(tensors) == args.shard_size:
                        flush()
            print(f"Batch {i}: {len(batch)} nodes, {segment_count} segments, "
                  f"{metrics['decoded_alignments']} alignments decoded", flush=True)
        flush()
    manifest.update(status="complete", tensors=total, shards=shard)
    write_json(out / "run_report.json", manifest)
    print(json.dumps(manifest, indent=2))


def build(args):
    if getattr(args, "format", "candidate-v3") == "legacy":
        return build_legacy(args)
    from indexed_gam_pipeline.build_v2 import build as candidate_build
    for key, default in (("rows", 200), ("width", 101), ("debug_rows", False)):
        if not hasattr(args, key):
            setattr(args, key, default)
    return candidate_build(args)


def positive(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def fraction(value):
    number = float(value)
    if not 0 <= number <= 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("index", "discover", "validate", "build"):
        sub = commands.add_parser(name)
        sub.add_argument("--gam", required=True)
        if name != "index":
            sub.add_argument("--chr", default="", help="Match Alignment.refpos.name, as in the original GAM pipeline")
        sub.add_argument("--output", required=True,
                         help="New .gai file for index; new or empty directory otherwise")
        if name == "discover":
            sub.add_argument("--min-mapq", type=int, default=5, help="Exclusive MAPQ threshold")
            sub.add_argument("--node-alt", type=fraction, default=0.05)
            sub.add_argument("--chr-nodes", help="Optional node whitelist text file")
            sub.add_argument("--max-alignments", type=positive, help="Exploratory partial scan ONLY")
            sub.add_argument("--max-nodes", type=positive, help="Limit selected nodes for smoke testing")
        elif name != "index":
            sub.add_argument("--nodes", required=True)
            sub.add_argument("--index", help="Default: GAM path + .gai")
        if name == "build":
            sub.add_argument("--format", choices=("candidate-v4", "candidate-v3", "candidate-v2", "legacy"), default="candidate-v4")
            sub.add_argument("--gbz", help="Matching GBZ graph for candidate-v4 occurrences")
            sub.add_argument("--gbz-query", help="Compiled gbz_node_counts helper")
            sub.add_argument("--occurrence-cache", help="Reusable GBZ occurrence SQLite cache")
            sub.add_argument("--walk-counts", help="Reusable SQLite node-to-distinct-W-count cache; required for candidate-v3")
            sub.add_argument("--rows", type=positive, default=200)
            sub.add_argument("--width", type=positive, default=101)
            sub.add_argument("--debug-rows", action="store_true", help="Include row paths and column graph coordinates")
            sub.add_argument("--gfa", help="Matching GFA with original vg node IDs")
            sub.add_argument("--node-sqlite", help="Matching graph node SQLite index (nodes.node_id, nodes.seq)")
            sub.add_argument("--node-json", help="node_id/sequence records, optionally with coordinates")
            sub.add_argument("--gam-cache-mb", type=int, default=64,
                             help="Bounded GAM group cache in MiB; 0 disables reuse (default: 64)")
            sub.add_argument("--workers", type=positive, default=1, help="Candidate compute processes; shared parent GAM cache")
            sub.add_argument("--early-alt-filter", action="store_true", help="Safely reject candidates below ALT support bound before overlap counting")
            sub.add_argument("--node-index-cache-nodes", type=int, default=0, help="FIFO node index cap across workers, 0 disables, maximum 1000")
            sub.add_argument("--node-index-cache-mb", type=int, default=64, help="Total estimated index container budget across workers")
            sub.add_argument("--batch-nodes", type=positive, default=128)
            sub.add_argument("--max-node-span", type=positive, default=10000)
            sub.add_argument("--max-batch-segments", type=positive, default=1000000)
            sub.add_argument("--shard-size", type=positive, default=4096)
            sub.add_argument("--max-tensors", type=positive,
                             help="Stop after writing this many tensors (for bounded test runs)")
            sub.add_argument("--min-mapq", type=int, default=10)
            sub.add_argument("--min-af", type=fraction, default=0.05)
            sub.add_argument("--min-variants", type=positive, default=3)
            sub.add_argument("--min-allele-bq", type=float, default=10)
            sub.add_argument("--max-indel-len", type=positive, default=50)
            sub.add_argument("--variant-type", choices=("snp", "indel", "all"), default="all")
    args = parser.parse_args()
    if args.command == "build" and not (args.gfa or args.node_json or args.node_sqlite):
        parser.error("build requires --gfa, --node-sqlite, or --node-json")
    try:
        {"index": index_gam, "discover": discover,
         "validate": validate, "build": build}[args.command](args)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()

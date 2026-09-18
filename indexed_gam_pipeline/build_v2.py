"""Bounded node batches, complete records, and versioned candidate tensor outputs."""
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from indexed_gam_pipeline.candidates import (VERSION, CHANNELS, V3_VERSION, V3_CHANNELS, BASES, OPS, ROW_ORDER,
                                             decode_alignment, overlap, make_tensor)
from indexed_gam_pipeline.gam_reader import IndexedGam
from indexed_gam_pipeline.segments import on_chromosome


def build(args):
    if getattr(args, "format", "candidate-v4") == "candidate-v4":
        from indexed_gam_pipeline.gbz_counts import GBZCounts
        if not all(getattr(args, key, None) for key in ("gbz", "gbz_query", "occurrence_cache")):
            raise ValueError("candidate-v4 requires --gbz, --gbz-query and --occurrence-cache")
        with GBZCounts(args.gbz, args.gbz_query, args.occurrence_cache) as lookup:
            return _build(args, lookup)
    if getattr(args, "format", "candidate-v3") == "candidate-v3":
        from indexed_gam_pipeline.walk_counts import WalkCounts
        if not getattr(args, "walk_counts", None):
            raise ValueError("candidate-v3 requires --walk-counts; build the cache once with walk_counts.py")
        with WalkCounts(args.walk_counts) as lookup:
            return _build(args, lookup)
    return _build(args, None)


def _build(args, walk_lookup):
    from indexed_gam_pipeline.run import load_nodes, node_records, batches, new_output, write_json
    if not 1 <= args.max_indel_len <= 50:
        raise ValueError("--max-indel-len must be in [1,50]")
    if args.width < args.max_indel_len:
        raise ValueError("--width must accommodate --max-indel-len")
    nodes = load_nodes(args.nodes)
    reader = IndexedGam(args.gam, args.index)
    out = new_output(args.output)
    (out / "target_nodes.txt").write_text("".join(f"{n}\n" for n in nodes))
    parameters = {k: getattr(args, k) for k in ("min_mapq", "min_af", "min_variants",
        "min_allele_bq", "max_indel_len", "variant_type", "rows", "width")}
    is_gbz = getattr(args, "format", "candidate-v4") == "candidate-v4"
    schema = 4 if is_gbz else (3 if walk_lookup is not None else 2)
    channels = V3_CHANNELS if walk_lookup is not None else CHANNELS
    version = V3_VERSION if walk_lookup is not None else VERSION
    if is_gbz:
        channels = CHANNELS + ["node_distinct_gbwt_path_count"]
        version = "indexed-gam-candidate-v4"
    manifest = dict(schema_version=schema, tensor_format_version=version, status="running",
        arguments=vars(args), parameters=parameters, shape=[len(channels), args.rows, args.width],
        dtype="int32" if walk_lookup is not None else "int16", channels=channels, encodings=dict(bases=BASES, padding=0,
        quality_without_read_base=-1, missing_quality=-1, operations=OPS,
        event_flags={"difference": 1, "candidate_region": 2}),
        coordinates="zero-based forward node; half-open intervals and candidate_columns",
        normalization="exact forward node/position/alleles; no repeat or cross-node normalization",
        insertion_overlap="closed boundary [mapping.start,mapping.end]; REF requires adjacent M/X columns on both sides",
        repeated_visits="one count and row per record; ALT then REF then other, earliest mapping breaks ties",
        row_selection="ALT then REF then other; descending MAPQ then record SHA256; identical records retained",
        row_order=ROW_ORDER,
        gai_version=reader.version, nodes=len(nodes), shards=0, tensors=0,
        unsupported_events=0, filtered_candidates=0, debug_rows=args.debug_rows)
    if walk_lookup is not None and not is_gbz:
        manifest["walk_count_lookup"] = dict(path=str(walk_lookup.path), **walk_lookup.metadata)
        manifest["encodings"]["node_distinct_w_record_count"] = {
            "value": "exact raw count, not normalized or clipped", "padding": 0,
            "node_absent_from_W": 0, "insertion_and_gap": "count of the row's anchor node",
            "missing_coverage": 0}
    if is_gbz:
        manifest["occurrence_lookup"] = dict(path=str(walk_lookup.path), **walk_lookup.metadata)
        manifest["encodings"]["node_distinct_gbwt_path_count"] = {
            "value": "exact distinct GBWT path count; either orientation; repeated visits counted once",
            "padding_and_missing_coverage": 0,
            "insertion_and_gap": "count of the row's anchor node",
            "missing_node_or_sequence_mismatch": "error"}
    write_json(out / "run_report.json", manifest)
    write_json(out / "manifest.json", manifest)
    tensors, metadata = [], []
    with (out / "variant_summary.ndjson").open("w") as summary, \
         (out / "unsupported_events.ndjson").open("w") as unsupported, \
         (out / "filtered_candidates.ndjson").open("w") as filtered:
        def flush():
            if not tensors:
                return
            shard = manifest["shards"]
            np.save(out / f"shard_{shard:05d}_data.npy", np.stack(tensors))
            for index, meta in enumerate(metadata):
                summary.write(json.dumps(dict(meta, schema_version=schema, channels=channels,
                    parameters=parameters, shard_index=shard, index_within_shard=index)) + "\n")
            manifest["shards"] += 1
            manifest["tensors"] += len(tensors)
            tensors.clear()
            metadata.clear()

        for bi, batch in enumerate(batches(nodes, args.batch_nodes, args.max_node_span), 1):
            wanted = set(batch)
            metrics = {}
            alignments = []
            context_nodes = set(wanted)
            for alignment in reader.fetch(wanted, metrics):
                if alignment.mapping_quality <= args.min_mapq or not on_chromosome(alignment, getattr(args, "chr", "")):
                    continue
                alignments.append(alignment)
                context_nodes.update(m.position.node_id for m in alignment.path.mapping)
                if len(alignments) > args.max_batch_segments:
                    raise ValueError("Complete-alignment batch limit exceeded; reduce --batch-nodes")
            print(f"Batch {bi}: retrieved {len(alignments)} alignments; loading {len(context_nodes)} graph nodes", flush=True)
            records = node_records(args, context_nodes)
            print(f"Batch {bi}: graph loaded; decoding original edits", flush=True)
            sequences = {n: r["sequence"] for n, r in records.items()}
            walk_counts = (walk_lookup.get_counts(context_nodes, sequences) if is_gbz else
                           walk_lookup.get_counts(context_nodes) if walk_lookup is not None else None)
            reads, candidates = [], set()
            by_node = defaultdict(list)
            for ai, alignment in enumerate(alignments):
                read, rejected = decode_alignment(alignment, sequences, args.max_indel_len)
                reads.append(read)
                for node in {v.node for v in read.visits} & wanted:
                    by_node[node].append(read)
                candidates.update(o.candidate for o in read.observations if o.candidate.node in wanted)
                for event in rejected:
                    unsupported.write(json.dumps(dict(event, record_sha256=read.digest,
                        in_target_nodes=event["node_id"] in wanted,
                        batch_index=bi, record_index=ai, read_name=read.name)) + "\n")
                    manifest["unsupported_events"] += 1
            print(f"Batch {bi}: counting {len(candidates)} candidates", flush=True)
            for candidate in sorted(candidates):
                eligible = []
                for read in by_node[candidate.node]:
                    hit = overlap(read, candidate, args.min_allele_bq)
                    if hit is not None:
                        eligible.append((read, *hit))
                counts = Counter(s for _, s, _ in eligible)
                af = counts["alt"] / len(eligible) if eligible else 0
                reasons = []
                if counts["alt"] < args.min_variants:
                    reasons.append("min_variants")
                if af < args.min_af:
                    reasons.append("min_af")
                if args.variant_type != "all" and (candidate.kind == "SNP") != (args.variant_type == "snp"):
                    reasons.append("variant_type")
                if reasons:
                    filtered.write(json.dumps(dict(candidate.metadata(), reasons=reasons,
                        coverage=len(eligible), alt_count=counts["alt"], ref_count=counts["ref"],
                        other_count=counts["other"], af=af)) + "\n")
                    manifest["filtered_candidates"] += 1
                    continue
                tensor, meta = make_tensor(candidate, eligible, args.rows, args.width, args.debug_rows, node_walk_counts=walk_counts)
                meta["tensor_format_version"] = version
                tensors.append(tensor)
                metadata.append(meta)
                if len(tensors) >= args.shard_size:
                    flush()
                if getattr(args, "max_tensors", None) is not None and manifest["tensors"] + len(tensors) >= args.max_tensors:
                    flush()
                    break
            print(f"Batch {bi}: {len(batch)} target nodes, {len(reads)} complete alignments, "
                  f"{len(context_nodes)} context nodes, {len(candidates)} candidates", flush=True)
            if getattr(args, "max_tensors", None) is not None and manifest["tensors"] >= args.max_tensors:
                break
        flush()
    manifest["status"] = "complete"
    if is_gbz:
        manifest["occurrence_lookup"] = dict(path=str(walk_lookup.path), **walk_lookup.metadata)
    write_json(out / "run_report.json", manifest)
    write_json(out / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))

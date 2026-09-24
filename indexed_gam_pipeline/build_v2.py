"""Bounded node batches, complete records, and versioned candidate tensor outputs."""
from collections import Counter, defaultdict
from contextlib import ExitStack, closing, nullcontext
import sqlite3
import json
import time
from pathlib import Path

import numpy as np
from indexed_gam_pipeline.tensor_storage import STORAGE_VERSION

from indexed_gam_pipeline.candidates import (VERSION, CHANNELS, V3_VERSION, V3_CHANNELS, BASES, OPS, ROW_ORDER, ROW_SELECTION_VERSION, WINDOW_ENCODING_VERSION,
                                             decode_alignment, overlap, make_tensor)
from indexed_gam_pipeline.gam_reader import IndexedGam
from indexed_gam_pipeline.candidate_work import DEFAULT_MAX_NODE_READS, SITE_UNIT
from indexed_gam_pipeline.segments import on_chromosome


def build(args):
    start = time.perf_counter()
    from indexed_gam_pipeline.split_outputs import split_enabled
    split_enabled(args)  # Validate paths/thresholds before creating any output.
    _dispatch(args)
    # Include initial cache hashing, helper startup and shutdown in wall time.
    out = Path(args.output)
    manifest = json.loads((out / "manifest.json").read_text())
    manifest["timing"]["total_wall_seconds"] = time.perf_counter() - start
    from indexed_gam_pipeline.run import write_json
    write_json(out / "manifest.json", manifest)
    write_json(out / "run_report.json", manifest)
    for folder in manifest.get('variant_outputs', {}).values():
        child = json.loads((Path(folder)/'manifest.json').read_text())
        child['timing']['total_wall_seconds'] = manifest['timing']['total_wall_seconds']
        write_json(Path(folder)/'manifest.json', child)
        write_json(Path(folder)/'run_report.json', child)
    print(json.dumps({"timing": manifest["timing"]}))


def _dispatch(args):
    if getattr(args, "graph_index", None):
        from indexed_gam_pipeline.graph_index import GraphIndex
        if getattr(args, "format", "candidate-v4") != "candidate-v4":
            raise ValueError("--graph-index requires candidate-v4")
        if any(getattr(args, key, None) for key in ("node_sqlite", "node_json", "gfa", "gbz", "gbz_query", "occurrence_cache", "walk_counts")):
            raise ValueError("--graph-index replaces legacy sequence and occurrence inputs")
        with GraphIndex(args.graph_index) as lookup:
            return _build(args, lookup)
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
    # Keep one read snapshot and page cache for the entire build. Reopening the
    # graph database per batch incurs expensive lock/I/O cycles on shared storage.
    with ExitStack() as stack:
        connection = None
        if getattr(args, "node_sqlite", None):
            uri = Path(args.node_sqlite).resolve().as_uri() + "?mode=ro"
            connection = stack.enter_context(closing(sqlite3.connect(uri, uri=True)))
            connection.execute("BEGIN")
        return _build_impl(args, walk_lookup, connection)


def _build_impl(args, walk_lookup, sequence_connection):
    from indexed_gam_pipeline.run import load_nodes, node_records, batches, new_output, write_json
    if not 1 <= args.max_indel_len <= 50:
        raise ValueError("--max-indel-len must be in [1,50]")
    if args.width < args.max_indel_len:
        raise ValueError("--width must accommodate --max-indel-len")
    if not 0 <= getattr(args, "node_index_cache_nodes", 0) <= 1000:
        raise ValueError("--node-index-cache-nodes must be in [0,1000]")
    if getattr(args, "node_index_cache_mb", 64) < 0:
        raise ValueError("--node-index-cache-mb must be nonnegative")
    if getattr(args, "max_node_reads", DEFAULT_MAX_NODE_READS) < 0:
        raise ValueError("--max-node-reads must be nonnegative")
    if getattr(args, "candidate_unit", "site") not in ("site", "allele"):
        raise ValueError("--candidate-unit must be site or allele")
    nodes = load_nodes(args.nodes)
    from indexed_gam_pipeline.edit_columns import select_decode_mode
    decode_policy = select_decode_mode(args.gam, getattr(args, 'decode_mode', 'auto'))
    if getattr(args, 'gam_reader', 'python') == 'vg':
        from indexed_gam_pipeline.vg_reader import VgGam
        reader = VgGam(args.gam, args.index, getattr(args, 'vg', '/scratch/jshen/bin/vg_v1.77.0'))
    else:
        reader = IndexedGam(args.gam, args.index, cache_bytes=getattr(args, "gam_cache_mb", 64) * 1024 * 1024)
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
        dtype="int8", tensor_storage_version=STORAGE_VERSION, channels=channels, encodings=dict(bases=BASES, padding=0,
        quality_without_read_base=-1, missing_quality=-1,
        base_quality="clip(raw quality, -1, 127); unavailable quality encodes as -1",
        mapping_quality="clip(raw MAPQ, -1, 127)", operations=OPS,
        event_flags={"difference": 1, "candidate_region": 2}),
        coordinates="zero-based forward node; half-open intervals and candidate_columns",
        normalization="exact forward node/position/alleles; no repeat or cross-node normalization",
        insertion_overlap="closed boundary [mapping.start,mapping.end]; REF requires adjacent M/X columns on both sides",
        repeated_visits="one count and row per record; ALT then REF then other, earliest mapping breaks ties",
        row_selection=ROW_ORDER, row_selection_version=ROW_SELECTION_VERSION,
        window_encoding_version=WINDOW_ENCODING_VERSION,
        row_order=ROW_ORDER,
        gai_version=reader.version, nodes=len(nodes), shards=0, tensors=0,
        unsupported_events=0, filtered_candidates=0, debug_rows=args.debug_rows)
    manifest['decode_policy'] = decode_policy
    if walk_lookup is not None and not is_gbz:
        manifest["walk_count_lookup"] = dict(path=str(walk_lookup.path), **walk_lookup.metadata)
        manifest["encodings"]["node_distinct_w_record_count"] = {
            "value": "min(exact raw count // 4, 127)", "padding": 0,
            "divisor": 4, "maximum_encoded_value": 127, "node_absent_from_W": 0, "insertion_and_gap": "count of the row's anchor node",
            "missing_coverage": 0}
    if is_gbz:
        manifest["occurrence_lookup"] = dict(path=str(walk_lookup.path), **walk_lookup.metadata)
        manifest["encodings"]["node_distinct_gbwt_path_count"] = {
            "value": "min(exact distinct GBWT path count // 4, 127); either orientation; revisits counted once",
            "divisor": 4, "maximum_encoded_value": 127, "saturation": "counts >=508 encode as 127",
            "padding_and_missing_coverage": 0,
            "insertion_and_gap": "count of the row's anchor node",
            "missing_node_or_sequence_mismatch": "error"}
    write_json(out / "run_report.json", manifest)
    write_json(out / "manifest.json", manifest)
    timings = defaultdict(float)
    worker_timings = Counter()
    manifest['candidate_optimization'] = dict(early_alt_filter=getattr(args,'early_alt_filter',False),
        workers=getattr(args,'workers',1), node_index_cache_nodes=getattr(args,'node_index_cache_nodes',0),
        node_index_cache_mb=getattr(args,'node_index_cache_mb',64),
        max_node_reads=getattr(args,'max_node_reads',DEFAULT_MAX_NODE_READS),
        max_node_reads_rule='records ordered by SHA-256 of the serialized GAM record; first N per target node',
        early_af_filter=getattr(args,'early_af_filter',True), early_rejected=0, early_af_rejected=0,
        memory_policy='Single parent GAM/GBZ cache; batch-scoped fork workers; total node FIFO cap shared by budget division')
    manifest['sample_unit'] = SITE_UNIT if getattr(args,'candidate_unit','site') == 'site' else 'allele'
    if manifest['sample_unit'] == SITE_UNIT:
        manifest['site_definition'] = ('one tensor per (node, start, SNV|INDEL); INS and DEL at one start share an '
            'indel site; alleles passing all filters are listed in alleles[] by ALT count; the tensor is the '
            'representative (first) allele\'s candidate tensor; no repeat normalization')
    tensors, metadata = [], []
    from indexed_gam_pipeline.split_outputs import SplitOutputs
    split = getattr(args, 'snv_output', None) is not None
    with (SplitOutputs(args, manifest, timings) if split else nullcontext()) as split_outputs, \
         (out / "variant_summary.ndjson").open("w") as summary, \
         (out / "unsupported_events.ndjson").open("w") as unsupported, \
         (out / "filtered_candidates.ndjson").open("w") as filtered:
        def flush():
            if split_outputs is not None:
                split_outputs.flush()
                manifest['timing'] = dict(timings)
                write_json(out/'manifest.json', manifest)
                write_json(out/'run_report.json', manifest)
                return
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
            summary.flush()
            manifest["timing"] = dict(timings)
            write_json(out / "manifest.json", manifest)
            write_json(out / "run_report.json", manifest)

        for bi, batch in enumerate(batches(nodes, args.batch_nodes, args.max_node_span), 1):
            batch_started = time.perf_counter()
            wanted = set(batch)
            metrics = {}
            alignments = []
            context_nodes = set(wanted)
            phase_start = time.perf_counter()
            for alignment in reader.fetch(wanted, metrics):
                if alignment.mapping_quality <= args.min_mapq or not on_chromosome(alignment, getattr(args, "chr", "")):
                    continue
                alignments.append(alignment)
                context_nodes.update(m.position.node_id for m in alignment.path.mapping)
                if len(alignments) > args.max_batch_segments:
                    raise ValueError("Complete-alignment batch limit exceeded; reduce --batch-nodes")
            timings["gam_fetch_seconds"] += time.perf_counter() - phase_start
            phase_start = time.perf_counter()
            print(f"Batch {bi}: retrieved {len(alignments)} alignments; loading {len(context_nodes)} graph nodes", flush=True)
            if getattr(args, "graph_index", None):
                records = walk_lookup.get_nodes(context_nodes)
                sequences = {n: r["sequence"] for n, r in records.items()}
                walk_counts = {n: r["distinct_path_count"] for n, r in records.items()}
                timings["graph_index_seconds"] += time.perf_counter() - phase_start
            else:
                records = node_records(args, context_nodes, sqlite_connection=sequence_connection)
                sequences = {n: r["sequence"] for n, r in records.items()}
                timings["graph_sequences_seconds"] += time.perf_counter() - phase_start
                phase_start = time.perf_counter()
                walk_counts = (walk_lookup.get_counts(context_nodes, sequences) if is_gbz else
                               walk_lookup.get_counts(context_nodes) if walk_lookup is not None else None)
                timings["occurrence_seconds"] += time.perf_counter() - phase_start
            print(f"Batch {bi}: graph loaded; decoding original edits", flush=True)
            phase_start = time.perf_counter()
            reads, candidates = [], set()
            by_node = defaultdict(list)
            for ai, alignment in enumerate(alignments):
                read, rejected = decode_alignment(alignment, sequences, args.max_indel_len,
                                                  mode=decode_policy['selected'], target_nodes=wanted)
                # All graph sequences are already loaded; the decoded Read owns
                # its context, so release this original protobuf immediately.
                alignments[ai] = None
                alignment = None
                reads.append(read)
                for node in {v.node for v in read.visits} & wanted:
                    by_node[node].append(read)
                candidates.update(o.candidate for o in read.observations if o.candidate.node in wanted)
                for event in rejected:
                    rejected_meta = dict(event, record_sha256=read.digest,
                        in_target_nodes=event["node_id"] in wanted,
                        batch_index=bi, record_index=ai, read_name=read.name)
                    unsupported.write(json.dumps(rejected_meta) + "\n")
                    if split_outputs is not None: split_outputs.record('unsupported', rejected_meta)
                    manifest["unsupported_events"] += 1
            timings["decode_edits_seconds"] += time.perf_counter() - phase_start
            phase_start = time.perf_counter()
            print(f"Batch {bi}: counting {len(candidates)} candidates", flush=True)
            from indexed_gam_pipeline.candidate_work import (alt_support_bounds, candidate_results,
                exact_coverage, af_threshold, group_sites, site_id)
            tasks=[];kept=defaultdict(list)
            t=time.perf_counter()
            alt_filter=getattr(args,'early_alt_filter',False); af_filter=getattr(args,'early_af_filter',True)
            sites=getattr(args,'candidate_unit','site')=='site'
            bounds=alt_support_bounds(reads,candidates,args.min_allele_bq) if alt_filter or af_filter else None
            # AF upper bound = ALT support bound / exact coverage, both over all records
            # (before any --max-node-reads cap). Without a cap it only removes candidates
            # process_node would reject for min_af.
            coverage=(exact_coverage([c for c in candidates if not alt_filter or bounds[c]>=args.min_variants],by_node)
                      if af_filter else None)
            for candidate in sorted(candidates):
                if alt_filter and bounds[candidate] < args.min_variants:
                    rejected_meta = dict(candidate.metadata(), reasons=['min_variants'],
                        alt_support_upper_bound=bounds[candidate], coverage_not_evaluated=True,
                        **({'site_id':site_id(candidate)} if sites else {}))
                    filtered.write(json.dumps(rejected_meta)+'\n')
                    if split_outputs is not None: split_outputs.record('filtered', rejected_meta)
                    manifest['filtered_candidates']+=1
                    manifest['candidate_optimization']['early_rejected']+=1
                elif coverage is not None and (coverage[candidate]==0 or
                        bounds[candidate]/coverage[candidate] < af_threshold(candidate,args)):
                    cov=coverage[candidate]
                    rejected_meta = dict(candidate.metadata(), reasons=['min_af'],
                        alt_support_upper_bound=bounds[candidate], coverage=cov,
                        af_upper_bound=bounds[candidate]/cov if cov else 0.0, support_not_evaluated=True,
                        **({'site_id':site_id(candidate)} if sites else {}))
                    filtered.write(json.dumps(rejected_meta)+'\n')
                    if split_outputs is not None: split_outputs.record('filtered', rejected_meta)
                    manifest['filtered_candidates']+=1
                    manifest['candidate_optimization']['early_af_rejected']+=1
                else:
                    kept[candidate.node].append(candidate)
            prefilter_seconds=time.perf_counter()-t
            for node,items in kept.items():
                units=group_sites(items) if sites else [(c,) for c in items]
                for offset in range(0,len(units),16):
                    tasks.append((node,units[offset:offset+16]))
            batch_worker_timings=Counter();cache_peaks=dict(nodes_per_worker=0,estimated_bytes_per_worker=0)
            stop=False
            with candidate_results(tasks,by_node,args,walk_counts) as completed:
                for results,compute_times,cache_stats in completed:
                    batch_worker_timings.update(compute_times)
                    cache_peaks['nodes_per_worker']=max(cache_peaks['nodes_per_worker'],cache_stats['peak_nodes'])
                    cache_peaks['estimated_bytes_per_worker']=max(cache_peaks['estimated_bytes_per_worker'],cache_stats['peak_estimated_bytes'])
                    for tensor,meta in results:
                        if tensor is None:
                            filtered.write(json.dumps(meta)+'\n')
                            if split_outputs is not None: split_outputs.record('filtered', meta)
                            manifest['filtered_candidates']+=1
                            continue
                        meta['tensor_format_version']=version
                        if split_outputs is not None:
                            split_outputs.append(tensor, meta)
                        else:
                            tensors.append(tensor);metadata.append(meta)
                            if len(tensors)>=args.shard_size:flush()
                        if getattr(args,'max_tensors',None) is not None and manifest['tensors']+(split_outputs.buffered if split_outputs is not None else len(tensors))>=args.max_tensors:
                            flush();stop=True;break
                    if stop:break
            worker_timings.update(batch_worker_timings)
            manifest['candidate_worker_summed_timing']=dict(worker_timings)
            timings["candidate_tensors_and_shard_writes_seconds"] += time.perf_counter() - phase_start
            with (out / "batch_timing.ndjson").open("a") as batch_log:
                batch_log.write(json.dumps(dict(batch=bi, first_node=batch[0], last_node=batch[-1],
                    target_nodes=len(batch), alignments=len(reads), context_nodes=len(context_nodes),
                    candidates=len(candidates), tensors_written=manifest["tensors"],
                    tensors_buffered=split_outputs.buffered if split_outputs is not None else len(tensors), elapsed_seconds=time.perf_counter()-batch_started,
                    cumulative_stage_seconds=dict(timings), gam_query=metrics,
                    prefilter_seconds=prefilter_seconds, candidate_worker_summed_timing=dict(batch_worker_timings),
                    node_index_cache_peaks=cache_peaks)) + "\n")
            print(f"Batch {bi}: {len(batch)} target nodes, {len(reads)} complete alignments, "
                  f"{len(context_nodes)} context nodes, {len(candidates)} candidates", flush=True)
            # Release decoded reads before fetching the next batch. Otherwise
            # the previous by_node/reads references survive the next GAM fetch
            # and graph lookup, adding both batches to the memory peak.
            by_node.clear()
            reads.clear()
            alignments.clear()
            read = alignment = None
            if getattr(args, "max_tensors", None) is not None and manifest["tensors"] >= args.max_tensors:
                break
        phase_start = time.perf_counter()
        flush()
        timings["candidate_tensors_and_shard_writes_seconds"] += time.perf_counter() - phase_start
    manifest["timing"] = dict(timings)
    if is_gbz:
        manifest["occurrence_performance"] = walk_lookup.performance
    manifest["gam_group_cache"] = reader.cache_stats
    manifest["status"] = "complete"
    if is_gbz:
        manifest["occurrence_lookup"] = dict(path=str(walk_lookup.path), **walk_lookup.metadata)
    if split_outputs is not None:
        split_outputs.finish(manifest)
    write_json(out / "run_report.json", manifest)
    write_json(out / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))

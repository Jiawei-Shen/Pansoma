"""The tensor builder: for each node batch, fetch -> graph -> decode -> filter -> encode -> shard.

One builder process owns one GAM reader (with its group cache) and one read-only
graph index connection for the whole run. Per batch it retains only the decoded
reads of that batch; everything is released before the next fetch.

Per batch, candidates go through two cheap prefilters over all records (ALT support
upper bound < min_variants; ALT bound / exact coverage < AF threshold), then each
target node's records are capped (--max-node-reads, deterministic by record digest),
and support is counted per unit: a site (default) or a single allele.

Outputs (per output directory):
    shard_XXXXX_data.npy        stacked (n, 8, rows, width) int8 tensors
    variant_summary.ndjson      one record per tensor: identity, counts, AF, row groups, shard location
                                (site mode: plus site_id and every passing allele in alleles[])
    filtered_candidates.ndjson  rejected candidate alleles with reasons and counts
    unsupported_events.ndjson   oversized indels and complex replacements seen while decoding
    manifest.json               format/encoding versions, parameters, provenance, counters, timing
    batch_timing.ndjson         per-batch stage seconds and counters (shared directory only)
    target_nodes.txt            the requested node list (shared directory only)

With --snv-output/--indel-output the GAM is read and decoded once; SNP tensors go
to the SNV directory and INS/DEL tensors to the INDEL directory, each with its own
AF threshold, shards, summary and manifest. --output then holds the shared manifest,
batch timings and the complete audit streams, but no shards.
"""
from collections import Counter, defaultdict
from copy import deepcopy
import json
from pathlib import Path
import time

import numpy as np

from indexed_gam_pipeline_v2.candidates import (FORMAT_VERSION, SCHEMA_VERSION, STORAGE_VERSION,
    ROW_SELECTION_VERSION, WINDOW_ENCODING_VERSION, ROW_ORDER, CHANNELS, BASES, OPS, STRAND, COUNT_LOG_SCALE,
    NodeReads, decode_alignment, alt_support_bounds, exact_coverage, make_tensor)
from indexed_gam_pipeline_v2.common import batches, load_nodes, new_output, on_chromosome, write_json
from indexed_gam_pipeline_v2.gam_reader import IndexedGam
from indexed_gam_pipeline_v2.graph_index import GraphIndex

KIND = {"SNP": "SNV", "INS": "INDEL", "DEL": "INDEL"}
SITE_UNIT = "site-v1"
DEFAULT_MAX_NODE_READS = 800
PARAMETERS = ("min_mapq", "min_af", "min_variants", "min_allele_bq", "max_indel_len",
              "variant_type", "rows", "width", "max_node_reads", "candidate_unit", "early_af_filter")
ALLELE_FIELDS = ("candidate_id", "start", "end", "ref", "alt", "event_type", "event_length",
                 "coverage", "alt_count", "ref_count", "other_count", "af")
SITE_DEFINITION = ("one tensor per (node, start, SNV|INDEL); INS and DEL at one start share an INDEL site; "
                   "every allele is filtered on its own; passing alleles are listed in alleles[] by ALT count "
                   "(ties: candidate order) and the tensor is the first (representative) allele's tensor; "
                   "no repeat/left normalization")
READ_CAP_RULE = ("records on a target node ordered by SHA-256 of the serialized GAM record; the first N are used "
                 "for support counting and rows of every candidate on that node; applied after both prefilters")


def validate_args(args):
    if not 1 <= args.max_indel_len <= 50:
        raise ValueError("--max-indel-len must be in [1, 50]")
    if args.width < args.max_indel_len:
        raise ValueError("--width must accommodate --max-indel-len")
    if args.gam_cache_mb < 0:
        raise ValueError("--gam-cache-mb must be nonnegative")
    if args.max_node_reads < 0:
        raise ValueError("--max-node-reads must be nonnegative")
    if args.candidate_unit not in ("site", "allele"):
        raise ValueError("--candidate-unit must be site or allele")
    split = [getattr(args, k, None) for k in ("snv_output", "indel_output", "snv_min_af", "indel_min_af")]
    if all(v is None for v in split):
        return False
    if any(v is None for v in split):
        raise ValueError("Split output needs --snv-output, --indel-output, --snv-min-af and --indel-min-af")
    if not all(0 <= v <= 1 for v in split[2:]):
        raise ValueError("Split AF thresholds must be in [0, 1]")
    if args.variant_type != "all":
        raise ValueError("Split output requires --variant-type all")
    paths = [Path(p).resolve() for p in (args.output, *split[:2])]
    if len(set(paths)) != 3 or any(a in b.parents for a in paths for b in paths if a != b):
        raise ValueError("Shared, SNV and INDEL output directories must be distinct and non-nested")
    return True


def site_key(candidate):
    """SNV and INDEL sites stay separate (separate outputs); INS and DEL at one start share an INDEL site."""
    return candidate.node, candidate.start, KIND[candidate.kind]


def site_id(candidate):
    return "%d:%d:%s" % site_key(candidate)


def candidate_units(candidates, site_unit):
    """Sorted candidates -> evaluation units in output order: one tuple per site, or per allele."""
    if not site_unit:
        return [(c,) for c in candidates]
    sites = defaultdict(list)
    for c in candidates:
        sites[site_key(c)].append(c)
    return [tuple(sites[k]) for k in sorted(sites)]


def capped_reads(reads, cap):
    """At most `cap` records of one node: those with the smallest record digests.

    Deterministic and independent of fetch order, so every candidate on the node and
    every rerun sees the same records. Shallower nodes are returned unchanged.
    """
    return sorted(reads, key=lambda r: r.digest)[:cap] if cap and len(reads) > cap else reads


def af_threshold(candidate, args):
    threshold = getattr(args, "snv_min_af" if candidate.kind == "SNP" else "indel_min_af", None)
    return args.min_af if threshold is None else threshold


def count_support(candidate, node_reads, args):
    """Classify every covering record once; returns (eligible, counts, failed filter reasons)."""
    eligible = node_reads.eligible(candidate, args.min_allele_bq)
    counts = Counter(s for _, s, _ in eligible)
    af = counts["alt"] / len(eligible) if eligible else 0
    reasons = []
    if counts["alt"] < args.min_variants:
        reasons.append("min_variants")
    if af < af_threshold(candidate, args):
        reasons.append("min_af")
    if args.variant_type != "all" and (candidate.kind == "SNP") != (args.variant_type == "snp"):
        reasons.append("variant_type")
    summary = dict(coverage=len(eligible), alt_count=counts["alt"], ref_count=counts["ref"],
                   other_count=counts["other"], af=af)
    return eligible, summary, reasons


def evaluate_unit(alleles, node_reads, args, path_counts):
    """Count support for every allele of one unit and encode the unit's tensor.

    Returns (rejected allele records, tensor, metadata); tensor and metadata are None
    when no allele passes. An allele unit holds one allele. A site unit holds every
    prefiltered allele of one (node, start, SNV|INDEL): each allele is filtered exactly
    as on its own, passing alleles are ranked by ALT count (ties: candidate order), and
    the first one's tensor is the site's tensor — byte-identical to that allele's own.
    """
    site = args.candidate_unit == "site"
    extra = dict(site_id=site_id(alleles[0])) if site else {}
    evaluated = [(c, *count_support(c, node_reads, args)) for c in alleles]
    rejected = [dict(c.metadata(), reasons=reasons, **summary, **extra)
                for c, _, summary, reasons in evaluated if reasons]
    passing = sorted((e for e in evaluated if not e[3]), key=lambda e: (-e[2]["alt_count"], e[0]))
    if not passing:
        return rejected, None, None
    candidate, eligible, _, _ = passing[0]
    tensor, meta = make_tensor(candidate, eligible, path_counts, args.rows, args.width, args.debug_rows)
    if site:
        listed = [{k: dict(c.metadata(), **summary)[k] for k in ALLELE_FIELDS} for c, _, summary, _ in passing]
        meta.update(sample_unit=SITE_UNIT, site_id=extra["site_id"], alleles=listed, allele_count=len(listed),
                    second_allele_af=listed[1]["af"] if len(listed) > 1 else 0.0)
    return rejected, tensor, meta


class OutputDir:
    """One output directory: manifest, NDJSON audit streams and (optionally) NPY shards."""

    def __init__(self, path, manifest, shard_size, tensors=True):
        self.path = new_output(path)
        self.manifest = manifest
        self.shard_size = shard_size
        self.tensors, self.metadata = [], []
        self.streams = {"filtered": (self.path / "filtered_candidates.ndjson").open("w"),
                        "unsupported": (self.path / "unsupported_events.ndjson").open("w")}
        if tensors:
            self.streams["summary"] = (self.path / "variant_summary.ndjson").open("w")
        self.save()

    @property
    def buffered(self):
        return len(self.tensors)

    def record(self, stream, meta):
        self.streams[stream].write(json.dumps(meta) + "\n")
        self.manifest["filtered_candidates" if stream == "filtered" else "unsupported_events"] += 1
        if stream == "filtered" and meta.get("coverage_not_evaluated"):
            self.manifest["early_rejected"] += 1
        elif stream == "filtered" and meta.get("support_not_evaluated"):
            self.manifest["early_af_rejected"] += 1

    def add(self, tensor, meta):
        self.tensors.append(tensor)
        self.metadata.append(meta)
        if len(self.tensors) >= self.shard_size:
            self.flush()

    def flush(self):
        if not self.tensors:
            return
        shard = self.manifest["shards"]
        np.save(self.path / f"shard_{shard:05d}_data.npy", np.stack(self.tensors))
        for i, meta in enumerate(self.metadata):
            self.streams["summary"].write(json.dumps(dict(meta, schema_version=SCHEMA_VERSION,
                channels=CHANNELS, parameters=self.manifest["parameters"],
                shard_index=shard, index_within_shard=i)) + "\n")
        self.streams["summary"].flush()
        self.manifest["shards"] += 1
        self.manifest["tensors"] += len(self.tensors)
        self.tensors.clear()
        self.metadata.clear()
        self.save()

    def save(self, **extra):
        self.manifest.update(extra)
        write_json(self.path / "manifest.json", self.manifest)

    def close(self):
        for stream in self.streams.values():
            stream.close()


def build(args):
    started = time.perf_counter()
    for key, default in (("index", None), ("chr", ""), ("debug_rows", False), ("max_tensors", None),
                         ("snv_output", None), ("indel_output", None), ("snv_min_af", None), ("indel_min_af", None),
                         ("max_node_reads", DEFAULT_MAX_NODE_READS), ("candidate_unit", "site"),
                         ("early_af_filter", True)):
        if not hasattr(args, key):
            setattr(args, key, default)
    split = validate_args(args)
    nodes = load_nodes(args.nodes)
    reader = IndexedGam(args.gam, args.index, cache_bytes=args.gam_cache_mb * 1024 * 1024)
    parameters = {k: getattr(args, k) for k in PARAMETERS}
    typed = {}
    with GraphIndex(args.graph_index) as graph:
        manifest = dict(status="running", schema_version=SCHEMA_VERSION,
            tensor_format_version=FORMAT_VERSION, tensor_storage_version=STORAGE_VERSION,
            shape=[len(CHANNELS), args.rows, args.width], dtype="int8", channels=CHANNELS,
            encodings=dict(bases=BASES, padding=0, quality_without_read_base=-1, missing_quality=-1,
                base_quality="clip(raw quality, -1, 127)", mapping_quality="clip(raw MAPQ, -1, 127)",
                operations=OPS,
                candidate_alt="ALT base of each candidate-region column (bases encoding; DEL columns = 6), "
                              "identical in every row; 0 outside the region or without evidence",
                node_distinct_gbwt_path_count=f"min(127, floor({COUNT_LOG_SCALE} * log2(count + 1) + 0.5)); "
                                              f"count ~ 2 ** (value / {COUNT_LOG_SCALE}) - 1; "
                                              "insertion and gap columns use the anchor node",
                strand=dict(STRAND, definition="anchor mapping orientation relative to the candidate node's "
                                               "forward strand; one value per row")),
            row_selection_version=ROW_SELECTION_VERSION, row_order=ROW_ORDER,
            window_encoding_version=WINDOW_ENCODING_VERSION,
            parameters=parameters, arguments=vars(args), gai_version=reader.version,
            graph_index=dict(path=str(graph.path), **graph.metadata),
            sample_unit=SITE_UNIT if args.candidate_unit == "site" else "allele",
            **(dict(site_definition=SITE_DEFINITION) if args.candidate_unit == "site" else {}),
            read_cap=dict(max_node_reads=args.max_node_reads, rule=READ_CAP_RULE),
            nodes=len(nodes), shards=0, tensors=0, filtered_candidates=0, early_rejected=0,
            early_af_rejected=0, unsupported_events=0, debug_rows=args.debug_rows, timing={})
        shared = OutputDir(args.output, manifest, args.shard_size, tensors=not split)
        (shared.path / "target_nodes.txt").write_text("".join(f"{n}\n" for n in nodes))
        try:
            if split:
                for kind, path, variant_type, af in (("SNV", args.snv_output, "snp", args.snv_min_af),
                                                     ("INDEL", args.indel_output, "indel", args.indel_min_af)):
                    m = deepcopy(manifest)
                    m["parameters"] = dict(parameters, variant_type=variant_type, min_af=af)
                    m["arguments"] = dict(manifest["arguments"], output=str(path), variant_type=variant_type, min_af=af)
                    m["shared_output"] = str(shared.path)
                    typed[kind] = OutputDir(path, m, args.shard_size)
                shared.save(output_layout="split", variant_outputs={k: str(v.path) for k, v in typed.items()})
            _build_batches(args, nodes, reader, graph, shared, typed, started)
        finally:
            for sink in (shared, *typed.values()):
                sink.close()
    summary = {k: shared.manifest[k] for k in ("tensors", "shards", "filtered_candidates", "unsupported_events")}
    if typed:
        summary["tensors_by_type"] = shared.manifest["tensors_by_type"]
    print(json.dumps(dict(summary, timing=shared.manifest["timing"])), flush=True)
    return shared.manifest


def _build_batches(args, nodes, reader, graph, shared, typed, started):
    """The batch loop. `typed` is {'SNV': OutputDir, 'INDEL': OutputDir} in split mode, else empty."""
    sinks = typed or {"ALL": shared}

    def record(stream, meta):
        """Audit records go to the shared directory and, in split mode, to the typed one too."""
        shared.record(stream, meta)
        if typed and meta.get("event_type") in KIND:
            typed[KIND[meta["event_type"]]].record(stream, meta)

    def written():
        return sum(s.manifest["tensors"] + s.buffered for s in sinks.values())

    timings = defaultdict(float)
    for bi, batch in enumerate(batches(nodes, args.batch_nodes, args.max_node_span), 1):
        batch_started = time.perf_counter()
        wanted = set(batch)
        # 1. Fetch every complete alignment touching the batch; collect all visited nodes.
        metrics = {}
        alignments = []
        context_nodes = set(wanted)
        for alignment in reader.fetch(wanted, metrics):
            if alignment.mapping_quality <= args.min_mapq or not on_chromosome(alignment, args.chr):
                continue
            alignments.append(alignment)
            context_nodes.update(m.position.node_id for m in alignment.path.mapping)
            if len(alignments) > args.max_batch_alignments:
                raise ValueError("Batch alignment limit exceeded; reduce --batch-nodes or raise --max-batch-alignments")
        timings["gam_fetch_seconds"] += time.perf_counter() - batch_started
        # 2. Sequences and path counts for the targets and all context nodes.
        t = time.perf_counter()
        records = graph.get_nodes(context_nodes)
        sequences = {n: r["sequence"] for n, r in records.items()}
        path_counts = {n: r["distinct_path_count"] for n, r in records.items()}
        timings["graph_index_seconds"] += time.perf_counter() - t
        # 3. Decode edits; keep candidate observations on target nodes only.
        t = time.perf_counter()
        reads, candidates, by_node = [], set(), defaultdict(list)
        for ai in range(len(alignments)):
            read, rejected = decode_alignment(alignments[ai], sequences, args.max_indel_len, target_nodes=wanted)
            alignments[ai] = None  # the decoded Read owns everything we still need
            reads.append(read)
            for node in {v.node for v in read.visits} & wanted:
                by_node[node].append(read)
            candidates.update(o.candidate for o in read.observations if o.candidate.node in wanted)
            for event in rejected:
                record("unsupported", dict(event, record_sha256=read.digest,
                    in_target_nodes=event["node_id"] in wanted, batch_index=bi,
                    record_index=ai, read_name=read.name))
        timings["decode_edits_seconds"] += time.perf_counter() - t
        # 4. Prefilters over ALL records of each node, before any read cap, with no overlap scan:
        #    ALT support upper bound < min_variants, then ALT bound / exact coverage < AF threshold.
        #    Both only drop candidates that support counting would reject for the same reason.
        t = time.perf_counter()
        site = args.candidate_unit == "site"
        bounds = alt_support_bounds(reads, candidates, args.min_allele_bq)
        ordered = sorted(candidates)
        coverage = (exact_coverage([c for c in ordered if bounds[c] >= args.min_variants], reads)
                    if args.early_af_filter else {})
        kept = []
        for candidate in ordered:
            extra = dict(site_id=site_id(candidate)) if site else {}
            if bounds[candidate] < args.min_variants:
                record("filtered", dict(candidate.metadata(), reasons=["min_variants"],
                    alt_support_upper_bound=bounds[candidate], coverage_not_evaluated=True, **extra))
            elif args.early_af_filter and (coverage[candidate] == 0 or
                                           bounds[candidate] / coverage[candidate] < af_threshold(candidate, args)):
                depth = coverage[candidate]
                record("filtered", dict(candidate.metadata(), reasons=["min_af"],
                    alt_support_upper_bound=bounds[candidate], coverage=depth,
                    af_upper_bound=bounds[candidate] / depth if depth else 0.0,
                    support_not_evaluated=True, **extra))
            else:
                kept.append(candidate)
        # 5. Per unit (site or allele), in sorted order: cap the node's records, count
        #    support for every allele, encode the representative. Units are sorted by node,
        #    so one NodeReads (the capped records and their oriented visit views) serves every
        #    unit of a node and is dropped when the next node starts.
        node_reads = None
        stop = False
        for unit in candidate_units(kept, site):
            node = unit[0].node
            if node_reads is None or node_reads.node != node:
                node_reads = NodeReads(node, capped_reads(by_node[node], args.max_node_reads), args.width)
            rejected_alleles, tensor, meta = evaluate_unit(unit, node_reads, args, path_counts)
            for rejection in rejected_alleles:
                record("filtered", rejection)
            if tensor is None:
                continue
            sinks[KIND[meta["event_type"]] if typed else "ALL"].add(tensor, meta)
            if args.max_tensors is not None and written() >= args.max_tensors:
                stop = True
                break
        timings["candidate_tensors_and_shard_writes_seconds"] += time.perf_counter() - t
        shared.save(timing=dict(timings))
        with (shared.path / "batch_timing.ndjson").open("a") as log:
            log.write(json.dumps(dict(batch=bi, first_node=batch[0], last_node=batch[-1],
                target_nodes=len(batch), alignments=len(reads), context_nodes=len(context_nodes),
                candidates=len(candidates), tensors_written=written(),
                elapsed_seconds=time.perf_counter() - batch_started,
                cumulative_stage_seconds=dict(timings), gam_query=metrics)) + "\n")
        print(f"Batch {bi}: {len(batch)} target nodes, {len(reads)} alignments, "
              f"{len(context_nodes)} context nodes, {len(candidates)} candidates, "
              f"{written()} tensors so far", flush=True)
        # 6. Release this batch before the next fetch so two batches never coexist in memory.
        by_node.clear()
        node_reads = None
        reads.clear()
        alignments.clear()
        if stop:
            break
    t = time.perf_counter()
    for sink in sinks.values():
        sink.flush()
    timings["candidate_tensors_and_shard_writes_seconds"] += time.perf_counter() - t
    timings["total_wall_seconds"] = time.perf_counter() - started
    final = dict(timing=dict(timings), gam_group_cache=reader.cache_stats,
                 graph_index_performance=graph.performance, status="complete")
    for sink in typed.values():
        sink.save(**final)
    if typed:
        shared.save(tensors=sum(s.manifest["tensors"] for s in typed.values()),
                    shards=sum(s.manifest["shards"] for s in typed.values()),
                    tensors_by_type={k: s.manifest["tensors"] for k, s in typed.items()})
    shared.save(**final)

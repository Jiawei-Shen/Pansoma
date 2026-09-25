"""The tensor builder: for each node batch, fetch -> graph -> decode -> filter -> encode -> shard.

One builder process owns one GAM reader (with its group cache) and one read-only
graph index connection for the whole run. Per batch it retains only the decoded
reads of that batch; everything is released before the next fetch.

Per batch, candidates go through two cheap prefilters over all records (ALT support
upper bound < min_variants; ALT bound / exact coverage < AF threshold), then each
target node's records are capped (--max-node-reads, deterministic by record digest),
and support is counted per site: every allele of one (node, start, SNV|INDEL), each
allele filtered on its own, all passing alleles in one site tensor.

The GAM is read and decoded once. SNP sites go to the SNV directory and INS/DEL sites
to the INDEL directory, each with its own AF threshold (--snv-min-af, --indel-min-af):

    SNV and INDEL directories (--snv-output, --indel-output)
        shard_XXXXX_data.npy        stacked (n, 8, rows, width) int8 tensors
        variant_summary.ndjson      one record per site tensor: identity, counts, AF, row groups, site_id,
                                    every passing allele in alleles[], shard location
        filtered_candidates.ndjson  rejected candidate alleles of that type, with reasons and counts
        unsupported_events.ndjson   oversized indels of that type seen while decoding
        manifest.json               format/encoding versions, parameters, provenance, counters, timing
    shared directory (--output; no shards)
        manifest.json               the same plus output_layout, variant_outputs and tensors_by_type
        filtered_candidates.ndjson  every audit record (the typed directories get those whose
        unsupported_events.ndjson   event_type is SNP, INS or DEL)
        batch_timing.ndjson         per-batch stage seconds and counters
        target_nodes.txt            the requested node list after --chromosomes
        displaced_nodes.tsv         node, records: nodes outside the batch that left-normalization moved
                                    an indel of a target node onto (orchestrate.py builds them in
                                    supplement rounds; see the README "supplement run")

Every stream is closed before any manifest is saved with status "complete".

Records are decoded by the native decoder when it is built and passes its checks
(--decoder auto, the default; native.py), else by candidates.decode_alignment; the
outputs are identical and the manifest's `decoder` records which one ran and why.
"""
from collections import Counter, defaultdict, deque
from copy import deepcopy
import json
from pathlib import Path
import time

import numpy as np

from .candidates import (FORMAT_VERSION, SCHEMA_VERSION, STORAGE_VERSION,
    ROW_SELECTION_VERSION, WINDOW_ENCODING_VERSION, ROW_ORDER, CHANNELS, BASES, OPS, STRAND, COUNT_LINEAR_MAX,
    NodeReads, decode_alignment, alt_support_bounds, exact_coverage, make_site_tensor)
from .common import batches, load_nodes, new_output, write_json
from .tensor_postprocessing.chr_index import select_nodes
from .gam_reader import IndexedGam
from .graph_index import GraphIndex
from .native import select_decoder

KIND ={"SNP": "SNV", "INS": "INDEL", "DEL": "INDEL"}
SITE_UNIT = "site-v2"
PARAMETERS = ("min_mapq", "min_af", "min_variants", "min_allele_bq", "max_indel_len",
              "variant_type", "rows", "width", "max_node_reads", "candidate_unit", "early_af_filter")
ALLELE_FIELDS = ("candidate_id", "start", "end", "ref", "alt", "event_type", "event_length", "path",
                 "coverage", "alt_count", "ref_count", "other_count", "af")
SITE_DEFINITION = ("one tensor per (node, start, SNV|INDEL) after indel left-normalization; INS and DEL at one "
                   "start share an INDEL site; every allele is filtered on its own; passing alleles are listed in "
                   "alleles[] by ALT count (ties: candidate order) as A1, A2, ...; every record covering any passing "
                   "allele is a row, labeled with the allele it carries (or REF/OTHER), and channel 2 spells that "
                   "allele over the site layout (longest insertion's slots + longest deletion's span); top-level "
                   "counts/AF are A1's (the representative)")
READ_CAP_RULE = ("records on a target node ordered by SHA-256 of the serialized GAM record; the first N are used "
                 "for support counting and rows of every candidate on that node; applied after both prefilters")
# manifest.arguments: v2's `run build` argparse namespace, key for key in v2's order. Options v3 removed
# are recorded with the only value production ever used (LEGACY_ARGUMENTS), so manifests stay byte-identical.
ARGUMENTS = ("command", "gam", "output", "nodes", "index", "graph_index", "snv_min_af", "indel_min_af", "snv_output",
             "indel_output", "debug_rows", "max_tensors", "variant_type", "rows", "width", "gam_cache_mb",
             "batch_nodes", "max_node_span", "max_batch_alignments", "shard_size", "min_mapq", "min_af",
             "min_variants", "min_allele_bq", "max_indel_len", "candidate_unit", "max_node_reads", "chromosomes",
             "chr_index", "early_af_filter", "decoder")
LEGACY_ARGUMENTS = {"max_tensors": None, "variant_type": "all", "candidate_unit": "site"}


def recorded(args, keys):
    """{key: value} for the manifest, in the order of `keys`; removed options take their fixed v2 value."""
    return {k: LEGACY_ARGUMENTS[k] if k in LEGACY_ARGUMENTS else getattr(args, k) for k in keys}


AUTO_SIZES = (512, 1024, 2048)  # --batch-nodes auto: the batch sizes tried at each position
AUTO_GAIN = 0.2  # a larger batch must cut the estimated GAM bytes per target node by at least this fraction


class BatchTooLarge(ValueError):
    pass


def fetch_estimate(reader, batch):
    """Compressed GAM bytes (BGZF block distance) the fetch of `batch` reads, from the GAI alone."""
    return sum(max(1, (end >> 16) - (start >> 16)) for start, end in reader.ranges(batch))


def adaptive_batches(nodes, reader, max_span):
    """--batch-nodes auto: yields (batch, plan). At each position the next 512, 1024 and 2048 target
    nodes (node-ID span limit scaled with the size) are priced by the GAM bytes their fetch would read;
    a larger batch is taken only if it reads at least AUTO_GAIN fewer bytes per target node than the
    smaller choice. Where neighbouring batches share the same GAM groups (long reads, coarse GAI bins)
    one large batch reads them once instead of several times; where they do not, batches stay small."""
    nodes = sorted(nodes)
    i = 0
    while i < len(nodes):
        best, estimates = None, {}
        for size in AUTO_SIZES:
            span = max_span * size // AUTO_SIZES[0]
            j = i + 1
            while j < len(nodes) and j - i < size and nodes[j] - nodes[i] <= span:
                j += 1
            per_node = fetch_estimate(reader, nodes[i:j]) / (j - i)
            estimates[str(size)] = round(per_node, 1)
            if best is None or per_node <= (1 - AUTO_GAIN) * best[1]:
                best = (j, per_node, size)
        yield nodes[i:best[0]], dict(size=best[2], bytes_per_node=estimates)
        i = best[0]


def validate_args(args):
    if args.batch_nodes != "auto" and not (isinstance(args.batch_nodes, int) and args.batch_nodes >= 1):
        raise ValueError("--batch-nodes must be a positive integer or auto")
    if not 1 <= args.max_indel_len <= 50:
        raise ValueError("--max-indel-len must be in [1, 50]")
    if args.width < args.max_indel_len:
        raise ValueError("--width must accommodate --max-indel-len")
    if args.gam_cache_mb < 1:
        raise ValueError("--gam-cache-mb must be positive")
    if args.max_node_reads < 0:
        raise ValueError("--max-node-reads must be nonnegative")
    split = [getattr(args, k, None) for k in ("snv_output", "indel_output", "snv_min_af", "indel_min_af")]
    if any(v is None for v in split):
        raise ValueError("Split output needs --snv-output, --indel-output, --snv-min-af and --indel-min-af")
    if not all(0 <= v <= 1 for v in split[2:]):
        raise ValueError("Split AF thresholds must be in [0, 1]")
    paths = [Path(p).resolve() for p in (args.output, *split[:2])]
    if len(set(paths)) != 3 or any(a in b.parents for a in paths for b in paths if a != b):
        raise ValueError("Shared, SNV and INDEL output directories must be distinct and non-nested")


def site_key(candidate):
    """SNV and INDEL sites stay separate (separate outputs); INS and DEL at one start share an INDEL site."""
    return candidate.node, candidate.start, KIND[candidate.kind]


def site_id(candidate):
    return "%d:%d:%s" % site_key(candidate)


def candidate_units(candidates):
    """Sorted candidates -> sites in output order: one tuple per site key, candidate order inside."""
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
    return args.snv_min_af if candidate.kind == "SNP" else args.indel_min_af


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
    summary = dict(coverage=len(eligible), alt_count=counts["alt"], ref_count=counts["ref"],
                   other_count=counts["other"], af=af)
    return eligible, summary, reasons


def evaluate_unit(alleles, node_reads, args, path_counts):
    """Count support for every allele of one site and encode the site's tensor.

    Returns (rejected allele records, tensor, metadata); tensor and metadata are None
    when no allele passes. A site holds every prefiltered allele of one (node, start,
    SNV|INDEL): each allele is filtered exactly as on its own, passing alleles are ranked by
    ALT count (ties: candidate order) and all of them go into one site tensor
    (make_site_tensor), the first as representative.
    """
    extra = dict(site_id=site_id(alleles[0]))
    evaluated = [(c, *count_support(c, node_reads, args)) for c in alleles]
    rejected = [dict(c.metadata(), reasons=reasons, **summary, **extra)
                for c, _, summary, reasons in evaluated if reasons]
    passing = sorted((e for e in evaluated if not e[3]), key=lambda e: (-e[2]["alt_count"], e[0]))
    if not passing:
        return rejected, None, None
    tensor, meta = make_site_tensor([c for c, _, _, _ in passing], [e for _, e, _, _ in passing], path_counts,
                                    args.rows, args.width, args.debug_rows)
    listed = [dict({k: dict(c.metadata(), **summary)[k] for k in ALLELE_FIELDS}, label=f"A{i + 1}")
              for i, (c, _, summary, _) in enumerate(passing)]
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
    validate_args(args)
    # The Python decoder is looked up at call time (tests patch it); it is also the native fallback.
    decode, decoder_info = select_decoder(args.decoder, python=lambda *a, **kw: decode_alignment(*a, **kw))
    print(f"Decoder: {decoder_info['used']}" + (f" ({decoder_info['reason']})" if "reason" in decoder_info else ""),
          flush=True)
    # Target nodes outside the chosen chromosome blocks are dropped before any batch (their reads are
    # still context of the remaining targets).
    nodes, selection = select_nodes(load_nodes(args.nodes), args.chromosomes, args.chr_index)
    nodes = nodes.tolist()
    if not nodes:
        raise ValueError(f"No target nodes left after --chromosomes {args.chromosomes}")
    reader = IndexedGam(args.gam, args.index, cache_bytes=args.gam_cache_mb * 1024 * 1024)
    parameters = recorded(args, PARAMETERS)
    typed = {}
    with GraphIndex(args.graph_index) as graph:
        manifest = dict(status="running", schema_version=SCHEMA_VERSION,
            tensor_format_version=FORMAT_VERSION, tensor_storage_version=STORAGE_VERSION,
            shape=[len(CHANNELS), args.rows, args.width], dtype="int8", channels=CHANNELS,
            encodings=dict(bases=BASES, padding=0, quality_without_read_base=-1, missing_quality=-1,
                base_quality="clip(raw quality, -1, 127)", mapping_quality="clip(raw MAPQ, -1, 127)",
                operations=OPS,
                site_allele="per row: the site allele the record carries (A1..Ak or REF) spelled over the "
                            "site layout columns (bases encoding, gaps = 6); 0 for OTHER records, outside "
                            "the site and without evidence",
                node_distinct_gbwt_path_count=f"count if count <= {COUNT_LINEAR_MAX}, else "
                                              f"min(127, {COUNT_LINEAR_MAX} + ceil(log2(count - {COUNT_LINEAR_MAX - 1}))); "
                                              "insertion and gap columns use the anchor node",
                strand=dict(STRAND, definition="anchor mapping orientation relative to the candidate node's "
                                               "forward strand; one value per row")),
            row_selection_version=ROW_SELECTION_VERSION, row_order=ROW_ORDER,
            window_encoding_version=WINDOW_ENCODING_VERSION,
            parameters=parameters, arguments=recorded(args, ARGUMENTS), gai_version=reader.version,
            graph_index=dict(path=str(graph.path), **graph.metadata),
            sample_unit=SITE_UNIT, site_definition=SITE_DEFINITION,
            read_cap=dict(max_node_reads=args.max_node_reads, rule=READ_CAP_RULE), decoder=decoder_info,
            nodes=len(nodes), chromosome_selection=selection, shards=0, tensors=0, filtered_candidates=0, early_rejected=0,
            early_af_rejected=0, unsupported_events=0, debug_rows=args.debug_rows, timing={})
        shared = OutputDir(args.output, manifest, args.shard_size, tensors=False)
        (shared.path / "target_nodes.txt").write_text("".join(f"{n}\n" for n in nodes))
        try:
            for kind, path, variant_type, af in (("SNV", args.snv_output, "snp", args.snv_min_af),
                                                 ("INDEL", args.indel_output, "indel", args.indel_min_af)):
                m = deepcopy(manifest)
                m["parameters"] = dict(parameters, variant_type=variant_type, min_af=af)
                m["arguments"] = dict(manifest["arguments"], output=str(path), variant_type=variant_type, min_af=af)
                m["shared_output"] = str(shared.path)
                typed[kind] = OutputDir(path, m, args.shard_size)
            shared.save(output_layout="split", variant_outputs={k: str(v.path) for k, v in typed.items()})
            _build_batches(args, nodes, reader, graph, shared, typed, started, decode)
        finally:
            for sink in (shared, *typed.values()):
                sink.close()
    summary = {k: shared.manifest[k] for k in ("tensors", "shards", "filtered_candidates", "unsupported_events")}
    summary["tensors_by_type"] = shared.manifest["tensors_by_type"]
    print(json.dumps(dict(summary, timing=shared.manifest["timing"])), flush=True)
    return shared.manifest


def _build_batches(args, nodes, reader, graph, shared, typed, started, decode):
    """The batch loop. `typed` is {'SNV': OutputDir, 'INDEL': OutputDir}; `decode` is
    decode_alignment or a native.NativeDecoder (select_decoder)."""
    sinks = typed

    def record(stream, meta):
        """Audit records go to the shared directory; SNP/INS/DEL records also to their typed one."""
        shared.record(stream, meta)
        if meta.get("event_type") in KIND:
            typed[KIND[meta["event_type"]]].record(stream, meta)

    def written():
        return sum(s.manifest["tensors"] + s.buffered for s in sinks.values())

    timings = defaultdict(float)
    displaced = Counter()  # node outside a batch -> records whose target-node indel normalization moved there
    auto = args.batch_nodes == "auto"
    planned = (adaptive_batches(nodes, reader, args.max_node_span) if auto
              else ((b, None) for b in batches(nodes, args.batch_nodes, args.max_node_span)))
    queue = deque()  # auto: halves of a batch that exceeded --max-batch-alignments
    bi = 0
    while True:
        if queue:
            batch, plan = queue.popleft()
        else:
            item = next(planned, None)
            if item is None:
                break
            batch, plan = item
        batch_started = time.perf_counter()
        wanted = set(batch)
        # 1. Fetch every complete alignment touching the batch; collect all visited nodes.
        metrics = {}
        alignments = []
        context_nodes = set(wanted)
        try:
            for alignment in reader.fetch(wanted, metrics):
                if alignment.mapping_quality <= args.min_mapq:
                    continue
                alignments.append(alignment)
                context_nodes.update(m.position.node_id for m in alignment.path.mapping)
                if len(alignments) > args.max_batch_alignments:
                    raise BatchTooLarge("Batch alignment limit exceeded; reduce --batch-nodes or raise "
                                        "--max-batch-alignments")
        except BatchTooLarge:
            if not auto or len(batch) < 2:
                raise
            half = len(batch) // 2  # auto: retry as two halves, in order
            queue.extendleft([(batch[half:], dict(plan, split=True)), (batch[:half], dict(plan, split=True))])
            alignments.clear()
            continue
        bi += 1
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
            read, rejected = decode(alignments[ai], sequences, args.max_indel_len, target_nodes=wanted)
            alignments[ai] = None  # the decoded Read owns everything we still need
            reads.append(read)
            for source, destination in read.moves:
                if source in wanted and destination not in wanted:
                    displaced[destination] += 1
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
        bounds = alt_support_bounds(reads, candidates, args.min_allele_bq)
        ordered = sorted(candidates)
        coverage = (exact_coverage([c for c in ordered if bounds[c] >= args.min_variants], reads)
                    if args.early_af_filter else {})
        kept = []
        for candidate in ordered:
            extra = dict(site_id=site_id(candidate))
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
        # 5. Per site, in sorted order: cap the node's records, count support for every
        #    allele, encode the site. Sites are sorted by node, so one NodeReads (the capped
        #    records and their oriented visit views) serves every site of a node and is
        #    dropped when the next node starts.
        node_reads = None
        for unit in candidate_units(kept):
            node = unit[0].node
            if node_reads is None or node_reads.node != node:
                node_reads = NodeReads(node, capped_reads(by_node[node], args.max_node_reads), args.width)
            rejected_alleles, tensor, meta = evaluate_unit(unit, node_reads, args, path_counts)
            for rejection in rejected_alleles:
                record("filtered", rejection)
            if tensor is None:
                continue
            sinks[KIND[meta["event_type"]]].add(tensor, meta)
        timings["candidate_tensors_and_shard_writes_seconds"] += time.perf_counter() - t
        shared.save(timing=dict(timings))
        with (shared.path / "batch_timing.ndjson").open("a") as log:
            log.write(json.dumps(dict(batch=bi, first_node=batch[0], last_node=batch[-1],
                target_nodes=len(batch), alignments=len(reads), context_nodes=len(context_nodes),
                candidates=len(candidates), tensors_written=written(),
                elapsed_seconds=time.perf_counter() - batch_started,
                cumulative_stage_seconds=dict(timings), gam_query=metrics,
                **({"batch_plan": plan} if plan is not None else {}))) + "\n")
        print(f"Batch {bi}: {len(batch)} target nodes, {len(reads)} alignments, "
              f"{len(context_nodes)} context nodes, {len(candidates)} candidates, "
              f"{written()} tensors so far", flush=True)
        # 6. Release this batch before the next fetch so two batches never coexist in memory.
        by_node.clear()
        node_reads = None
        reads.clear()
        alignments.clear()
    (shared.path / "displaced_nodes.tsv").write_text("".join(f"{n}\t{c}\n" for n, c in sorted(displaced.items())))
    shared.manifest["displaced_nodes"] = len(displaced)
    t = time.perf_counter()
    for sink in sinks.values():
        sink.flush()
    timings["candidate_tensors_and_shard_writes_seconds"] += time.perf_counter() - t
    timings["total_wall_seconds"] = time.perf_counter() - started
    decoder = dict(shared.manifest["decoder"])
    if hasattr(decode, "fallbacks"):  # records the native decoder raised on and Python decoded
        decoder["native_record_fallbacks"] = decode.fallbacks
    final = dict(timing=dict(timings), gam_group_cache=reader.cache_stats,
                 graph_index_performance=graph.performance, decoder=decoder, status="complete")
    for sink in (shared, *typed.values()):  # every stream is complete before any manifest says so
        sink.close()
    for sink in typed.values():
        sink.save(**final)
    shared.save(tensors=sum(s.manifest["tensors"] for s in typed.values()),
                shards=sum(s.manifest["shards"] for s in typed.values()),
                tensors_by_type={k: s.manifest["tensors"] for k, s in typed.items()})
    shared.save(**final)

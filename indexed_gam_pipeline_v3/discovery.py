"""Discovery: select target nodes from one full scan of a sorted GAM.

A node is selected when more than `node_alt` of its MAPQ-passing mappings carry an edit
(the node_stats.json counts `perfect` / `not_perfect` per node). Two rules:

    raw         (--raw) vg's edits as written (the rule of every run before 2026-09-24)
    normalized  (default) each record is first decoded and its indels left-normalized
                exactly like the builder does, so an indel counts on the node the builder will
                see it on. Target lists from this rule already contain the nodes that
                normalization moves indels onto, which otherwise need supplement rounds.

Normalized is the default (`--raw` for the old rule). The native module scans the GAM in parallel: the file is cut at group starts taken from the
GAI into contiguous segments, every record is read exactly once by one worker, and the
per-node counts are merged in order of first appearance, so node_stats.json and
target_nodes.txt are byte-identical to the sequential scan. Without the native module (or
for an exploratory --max-alignments scan) the sequential pure-Python scan runs (raw rule).
"""
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
import os
from pathlib import Path
import sys

import numpy as np
import pysam

from .common import new_output, write_json
from .gam_reader import IndexedGam, group, scan_gam

SEGMENTS_PER_PROCESS = 4
CHUNK = 1_000_000


def default_processes():
    return int(os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 1)


# --- sequential pure-Python scan (raw rule) -------------------------------------------

def python_counts(gam, min_mapq, max_alignments=None):
    """(nodes, perfect, not_perfect, max_read_length, counters) in first-appearance order."""
    stats = defaultdict(lambda: [0, 0, 0])  # perfect, imperfect, max read length
    count = used = 0
    for alignment in scan_gam(gam, max_alignments):
        count += 1
        if alignment.mapping_quality <= min_mapq:
            continue
        used += 1
        for mapping in alignment.path.mapping:
            nid = mapping.position.node_id
            if not nid:
                continue
            imperfect = any(e.from_length != e.to_length or e.sequence for e in mapping.edit)
            stats[nid][int(imperfect)] += 1
            stats[nid][2] = max(stats[nid][2], len(alignment.sequence))
        if count % 100000 == 0:
            print(f"Scanned {count:,} alignments; {len(stats):,} nodes", flush=True)
    columns = list(zip(*((n, p, q, r) for n, (p, q, r) in stats.items()))) or [[], [], [], []]
    arrays = [np.array(c, dtype=np.int64) for c in columns]
    return (*arrays, dict(alignments=count, used=used, fallbacks=0))


# --- parallel native scan -------------------------------------------------------------

def segments(gam, index, count):
    """[(start, end)] virtual offsets of `count` (or fewer) contiguous segments that tile the
    GAM at group starts (GAI run starts); the first starts at 0, the last ends at EOF (None)."""
    reader = IndexedGam(gam, index, cache_bytes=1)
    starts = sorted({start for _, _, runs in reader.bins for start, _ in runs} - {0})
    size = Path(gam).stat().st_size
    cuts, k = [0], 0
    for i in range(1, count):
        target = size * i // count
        while k < len(starts) and (starts[k] >> 16) < target:
            k += 1
        if k < len(starts) and starts[k] > cuts[-1]:
            cuts.append(starts[k])
    return list(zip(cuts, cuts[1:] + [None]))


def _scan_segment(job):
    gam, (start, end), normalized, min_mapq, max_indel, graph_index = job
    from . import native
    module, info = native.load()
    if module is None:
        raise RuntimeError(f"native module unavailable in a discovery worker: {info.get('reason')}")
    counter = module.Discovery()
    graph = None
    if normalized:
        from .graph_index import GraphIndex
        graph = GraphIndex(graph_index)
    groups = 0
    with pysam.BGZFile(str(gam), "rb") as stream:
        stream.seek(start)
        while end is None or stream.tell() < end:
            messages = group(stream)
            if messages is None:
                break
            groups += 1
            if not messages:
                continue
            if normalized:
                nodes = [n for n in module.group_nodes(messages).tolist() if n]
                sequences = {n: r["sequence"] for n, r in graph.get_nodes(nodes).items()}
                counter.add_normalized(messages, sequences, min_mapq, max_indel)
            else:
                counter.add_raw(messages, min_mapq)
        if end is not None and stream.tell() != end:
            raise ValueError(f"discovery segment {start} did not end at the next segment start {end}")
    if graph is not None:
        graph.db.close()
    nodes, perfect, not_perfect, max_len, counters = counter.result()
    return nodes, perfect, not_perfect, max_len, dict(counters, groups=groups)


def merge(results):
    """Per-segment counts (in segment order) -> one set in order of first appearance."""
    nodes = np.concatenate([r[0] for r in results]) if results else np.zeros(0, dtype=np.int64)
    if not len(nodes):
        empty = np.zeros(0, dtype=np.int64)
        return empty, empty, empty, empty
    unique, first, inverse = np.unique(nodes, return_index=True, return_inverse=True)
    perfect = np.zeros(len(unique), dtype=np.int64)
    not_perfect = np.zeros(len(unique), dtype=np.int64)
    max_len = np.zeros(len(unique), dtype=np.int64)
    np.add.at(perfect, inverse, np.concatenate([r[1] for r in results]))
    np.add.at(not_perfect, inverse, np.concatenate([r[2] for r in results]))
    np.maximum.at(max_len, inverse, np.concatenate([r[3] for r in results]))
    order = np.argsort(first, kind="stable")
    return unique[order], perfect[order], not_perfect[order], max_len[order]


def native_counts(gam, index, processes, normalized=False, min_mapq=5, max_indel=50, graph_index=None):
    jobs = [(str(gam), segment, normalized, min_mapq, max_indel, graph_index)
            for segment in segments(gam, index, max(1, processes) * SEGMENTS_PER_PROCESS)]
    if processes <= 1:
        results = [_scan_segment(job) for job in jobs]
    else:
        with ProcessPoolExecutor(processes, mp_context=multiprocessing.get_context("fork")) as pool:
            results = list(pool.map(_scan_segment, jobs))
    counters = {k: sum(r[4][k] for r in results) for k in ("alignments", "used", "fallbacks", "groups")}
    return (*merge(results), dict(counters, segments=len(jobs)))


# --- outputs ---------------------------------------------------------------------------

def write_node_stats(path, nodes, perfect, not_perfect, max_len):
    """node_stats.json exactly as write_json(dict) would write it, streamed (tens of millions of nodes)."""
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w") as stream:
        if not len(nodes):
            stream.write("{}\n")
        else:
            stream.write("{\n")
            for start in range(0, len(nodes), CHUNK):
                stop = start + CHUNK
                part = ",\n".join(
                    f'  "{n}": {{\n    "perfect": {p},\n    "not_perfect": {q},\n    "max_read_length": {r}\n  }}'
                    for n, p, q, r in zip(nodes[start:stop].tolist(), perfect[start:stop].tolist(),
                                          not_perfect[start:stop].tolist(), max_len[start:stop].tolist()))
                stream.write((",\n" if start else "") + part)
            stream.write("\n}\n")
    temporary.replace(path)


def discover(args):
    """Select nodes where more than `node_alt` of MAPQ-passing mappings carry an edit."""
    from . import native
    normalized = not getattr(args, "raw", False)
    max_alignments = getattr(args, "max_alignments", None)
    processes = getattr(args, "processes", None) or default_processes()
    max_indel = getattr(args, "max_indel_len", 50)
    graph_index = getattr(args, "graph_index", None)
    if normalized and not graph_index:
        raise ValueError("discover needs --graph-index (node sequences for normalization), or --raw")
    if normalized and max_alignments:
        raise ValueError("--max-alignments is for exploratory scans with --raw")
    module, info = native.load()
    out = new_output(args.output)
    if max_alignments or module is None:
        if normalized:
            raise ValueError(f"normalized discovery needs the native module ({info.get('reason')}); use --raw")
        if module is None:
            print(f"Discovery in pure Python, one process ({info.get('reason')})", flush=True)
        nodes, perfect, not_perfect, max_len, counters = python_counts(args.gam, args.min_mapq, max_alignments)
        engine, processes, counters["segments"] = "python", 1, 1
    else:
        index = getattr(args, "index", None)
        nodes, perfect, not_perfect, max_len, counters = native_counts(
            args.gam, index, processes, normalized, args.min_mapq, max_indel, graph_index)
        engine = "native"
    with np.errstate(divide="ignore", invalid="ignore"):
        chosen = (not_perfect >= 1) & (not_perfect / (perfect + not_perfect) > args.node_alt)
    selected = np.sort(nodes[chosen])
    (out / "target_nodes.txt").write_text("".join(f"{n}\n" for n in selected.tolist()))
    write_node_stats(out / "node_stats.json", nodes, perfect, not_perfect, max_len)
    report = dict(gam=str(Path(args.gam).resolve()), alignments_scanned=counters["alignments"],
                  nodes_observed=int(len(nodes)), nodes_selected=int(len(selected)), min_mapq=args.min_mapq,
                  node_alt=args.node_alt, exploratory=bool(max_alignments), max_alignments=max_alignments,
                  rule="normalized" if normalized else "raw", engine=engine, processes=processes,
                  segments=counters["segments"], alignments_passing_mapq=counters["used"])
    if normalized:
        report.update(max_indel_len=max_indel, graph_index=str(Path(graph_index).resolve()),
                      normalization_fallbacks=counters["fallbacks"])
    write_json(out / "discovery_report.json", report)
    print(json.dumps(report, indent=2))
    sys.stdout.flush()

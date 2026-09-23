#!/usr/bin/env python3
"""Independent audit of a --debug-rows tensor folder against the source GAM and graph index.

Re-queries the original records with the GAI, recounts coverage from raw mapping
intervals (without using candidates.overlap), and checks every selected row's
graph bases, path-count cells, candidate-ALT stripe, strand, grouping and the
recorded row sampling. When the build capped records per node (max_node_reads),
the same cap is re-applied independently from the record digests. Requires
tensors built with --debug-rows.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from indexed_gam_pipeline_v2.candidates import (BASES, FORMAT_VERSION, ROW_SELECTION_VERSION, STORAGE_VERSION,
                                                STRAND, WINDOW_ENCODING_VERSION, Candidate, candidate_alt_codes,
                                                encode_count, rc)
from indexed_gam_pipeline_v2.common import on_chromosome, write_json
from indexed_gam_pipeline_v2.gam_reader import IndexedGam
from indexed_gam_pipeline_v2.graph_index import GraphIndex


def check(condition, message):
    if not condition:
        raise ValueError(message)


def independent_coverage(gam, index, metadata, sequences, min_mapq, chromosome=""):
    """From raw mapping intervals of every fetched record:
    ({candidate_id: Counter(record sha256)}, {node: [sha256 of every record mapped to it]}, query)."""
    covered = {m["candidate_id"]: Counter() for m in metadata}
    on_node = defaultdict(list)
    targets = {m["node_id"] for m in metadata}
    query = {}
    for alignment in IndexedGam(gam, index).fetch(targets, query):
        if alignment.mapping_quality <= min_mapq or not on_chromosome(alignment, chromosome):
            continue
        digest = hashlib.sha256(alignment.SerializeToString(deterministic=True)).hexdigest()
        for node in {m.position.node_id for m in alignment.path.mapping} & targets:
            on_node[node].append(digest)
        for meta in metadata:
            hit = False
            for mapping in alignment.path.mapping:
                if mapping.position.node_id != meta["node_id"]:
                    continue
                if not any(e.from_length or e.to_length for e in mapping.edit):
                    continue
                length = len(sequences[meta["node_id"]])
                span = sum(e.from_length for e in mapping.edit)
                start = mapping.position.offset
                if mapping.position.is_reverse:
                    start = length - start - span
                end = start + span
                hit = (start <= meta["start"] <= end if meta["event_type"] == "INS"
                       else start < meta["end"] and end > meta["start"])
                if hit:
                    break
            if hit:
                covered[meta["candidate_id"]][digest] += 1
    return covered, on_node, query


def cap_records(digests, cap):
    """The builder's per-node read cap, re-applied: the `cap` smallest digests (a multiset), or None."""
    return Counter(sorted(digests)[:cap]) if cap and len(digests) > cap else None


def validate(folder, gam, index, graph_index):
    folder = Path(folder)
    manifest = json.loads((folder / "manifest.json").read_text())
    metadata = [json.loads(s) for s in (folder / "variant_summary.ndjson").read_text().splitlines()]
    check(manifest.get("debug_rows"), "Validation requires tensors built with --debug-rows")
    check(manifest.get("tensor_format_version") == FORMAT_VERSION, "Unexpected tensor format version")
    check(manifest.get("tensor_storage_version") == STORAGE_VERSION, "Unexpected storage version")
    graph_nodes = {m["node_id"] for m in metadata} | {
        c["node_id"] for m in metadata for row in m["rows"] for c in row["columns"] if c}
    with GraphIndex(graph_index) as lookup:
        records = lookup.get_nodes(graph_nodes)
    sequences = {n: r["sequence"] for n, r in records.items()}
    counts = {n: r["distinct_path_count"] for n, r in records.items()}
    parameters = manifest["parameters"]
    covered, on_node, query = independent_coverage(gam, index, metadata, sequences, parameters["min_mapq"],
                                                   manifest["arguments"].get("chr", ""))
    cap = parameters.get("max_node_reads", 0)
    kept = {node: cap_records(digests, cap) for node, digests in on_node.items()}
    capped_nodes = sum(k is not None for k in kept.values())
    results = []
    for meta in metadata:
        x = np.load(folder / f"shard_{meta['shard_index']:05d}_data.npy", mmap_mode="r")[meta["index_within_shard"]]
        label = meta["candidate_id"]
        n = meta["selected_alignments"]
        if kept.get(meta["node_id"]) is not None:  # records beyond the cap were never considered
            covered[label] = covered[label] & kept[meta["node_id"]]
        check(meta["window_encoding_version"] == WINDOW_ENCODING_VERSION, f"{label}: window encoding")
        check(meta["row_selection_version"] == ROW_SELECTION_VERSION, f"{label}: row selection version")
        check(list(x.shape) == manifest["shape"] and str(x.dtype) == manifest["dtype"], f"{label}: tensor shape/dtype")
        check(len(meta["rows"]) == n, f"{label}: debug row count")
        check(not x[:, n:].any(), f"{label}: unused row padding")
        check(meta["coverage"] == sum(meta[f"{c}_count"] for c in ("alt", "ref", "other")), f"{label}: counts")
        check(meta["af"] == meta["alt_count"] / meta["coverage"], f"{label}: AF")
        check(sum(covered[label].values()) == meta["coverage"], f"{label}: independent source coverage")
        selected = Counter(row["record_sha256"] for row in meta["rows"])
        check(not selected - covered[label], f"{label}: selected source records")
        if n == meta["coverage"]:
            check(selected == covered[label], f"{label}: full source record multiset")
        check(Counter(row["support"] for row in meta["rows"]) == Counter(meta["selected_counts"]),
              f"{label}: selected support counts")
        check(bool(np.all(x[1, :n][x[0, :n] == 6] == -1)), f"{label}: gap quality")
        lo, hi = meta["candidate_columns"]
        alt_codes = candidate_alt_codes(Candidate(meta["node_id"], meta["start"], meta["ref"], meta["alt"], meta["event_type"]))
        check(lo == meta["anchor_column"] and hi - lo == min(len(alt_codes), manifest["shape"][2] - lo),
              f"{label}: candidate columns")
        actual_paths, checked = [], 0
        for ri, row in enumerate(meta["rows"]):
            strand = STRAND["reverse"] if row["reversed_for_candidate"] else STRAND["forward"]
            path, previous = [], None
            for ci, col in enumerate(row["columns"]):
                if col is None:
                    check(not x[:, ri, ci].any(), f"{label}: all-zero missing evidence")
                    continue
                check(int(x[2, ri, ci]) == (alt_codes[ci - lo] if lo <= ci < hi else 0),
                      f"{label}: candidate ALT channel at row {ri}, column {ci}")
                check(int(x[7, ri, ci]) == strand, f"{label}: strand channel at row {ri}, column {ci}")
                if ci == meta["anchor_column"] and "offset" in col:
                    check(col["mapping_index"] == row["anchor_mapping_index"] and col["offset"] == meta["start"],
                          f"{label}: anchor coordinate")
                check(int(x[6, ri, ci]) == encode_count(counts[col["node_id"]]),
                      f"{label}: path count at row {ri}, column {ci}")
                if col["mapping_index"] != previous:
                    path.append((col["node_id"], col["reverse"]))
                    previous = col["mapping_index"]
                if "offset" in col:
                    base = sequences[col["node_id"]][col["offset"]]
                    if col["reverse"]:
                        base = rc(base)
                    check(int(x[5, ri, ci]) == BASES.get(base, 5), f"{label}: graph base at row {ri}, column {ci}")
                    checked += 1
            actual_paths.append(tuple(path))
        # Recompute the ranking, grouping and uniform sampling from the recorded audit.
        audit = meta["selection_audit"]
        check(Counter(a["record_sha256"] for a in audit) == covered[label], f"{label}: all preselection records")
        ordered_groups = {}
        for a in sorted(audit, key=lambda a: (-a["window_mismatch_bp"], -a["mapping_quality"],
                                              a["record_sha256"], a["anchor_mapping_index"])):
            ordered_groups.setdefault(tuple(tuple(p) for p in a["path"]), []).append(a)
        ranked = [a for group in ordered_groups.values() for a in group]
        k = min(manifest["shape"][1], len(ranked))
        indices = ([len(ranked) // 2] if k == 1 else
                   [i * (len(ranked) - 1) // (k - 1) for i in range(k)] if k else [])
        check(meta["selected_grouped_ranks"] == indices, f"{label}: uniform sampling ranks")
        check(n == k, f"{label}: sampled depth")
        for ri, idx in enumerate(indices):
            a, row = ranked[idx], meta["rows"][ri]
            check(row["record_sha256"] == a["record_sha256"] and row["anchor_mapping_index"] == a["anchor_mapping_index"]
                  and row["support"] == a["support"], f"{label}: sampled record ordering")
            check(actual_paths[ri] == tuple(tuple(p) for p in a["path"]), f"{label}: sampled path")
            mismatch = int(np.count_nonzero((x[0, ri] != 0) & (x[5, ri] != 0) & (x[0, ri] != x[5, ri]) & (x[4, ri] != 6)))
            check(mismatch == a["window_mismatch_bp"] == row["window_mismatch_bp"] == meta["window_mismatch_bp"][ri],
                  f"{label}: visible edit bp count")
        next_row = 0
        for group in meta["row_groups"]:
            check(group["start_row"] == next_row and group["end_row"] > next_row, f"{label}: group boundaries")
            key = tuple((v["node_id"], v["reverse"]) for v in group["path"])
            check(all(p == key for p in actual_paths[next_row:group["end_row"]]), f"{label}: group membership")
            next_row = group["end_row"]
        check(next_row == n, f"{label}: all rows grouped")
        results.append(dict(candidate_id=label, event_type=meta["event_type"], coverage=meta["coverage"],
                            selected_alignments=n, node_path_groups=len(meta["row_groups"]),
                            reference_columns_checked=checked, passed=True))
    check(len(results) == manifest["tensors"], "Manifest tensor count")
    return dict(passed=True, examples=len(results), event_types=dict(Counter(r["event_type"] for r in results)),
                reference_columns_checked=sum(r["reference_columns_checked"] for r in results),
                selected_rows_checked=sum(r["selected_alignments"] for r in results),
                indexed_query=query, source_gam=str(gam), source_index=index, graph_index=str(graph_index),
                max_node_reads=cap, capped_nodes=capped_nodes, candidates=results)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", help="tensor output directory built with --debug-rows")
    parser.add_argument("--gam", required=True)
    parser.add_argument("--index", help="default: GAM path + .gai")
    parser.add_argument("--graph-index", required=True)
    parser.add_argument("--output", required=True, help="audit report JSON")
    args = parser.parse_args(argv)
    report = validate(args.folder, args.gam, args.index, args.graph_index)
    write_json(args.output, report)
    print(json.dumps({k: v for k, v in report.items() if k != "candidates"}, indent=2))


if __name__ == "__main__":
    main()

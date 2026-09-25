#!/usr/bin/env python3
"""Independent audit of a --debug-rows tensor folder against the source GAM and graph index.

Re-queries the original records with the GAI, recounts every site allele's coverage from
raw mapping intervals (without using NodeReads), and checks every selected row's
graph bases, path-count cells, site-allele channel (re-derived from the recorded site
layout, graph sequence and row labels), strand, allele blocks and the recorded row
sampling. When the build capped records per node (max_node_reads), the same cap is
re-applied independently from the record digests. Requires tensors built with --debug-rows.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np

from ..candidates import (BASES, FORMAT_VERSION, ROW_SELECTION_VERSION, STORAGE_VERSION,
                                                STRAND, WINDOW_ENCODING_VERSION, encode_count, rc)
from ..common import write_json
from ..gam_reader import IndexedGam
from ..graph_index import GraphIndex


def check(condition, message):
    if not condition:
        raise ValueError(message)


def independent_coverage(gam, index, metadata, sequences, min_mapq):
    """From raw mapping intervals of every fetched record:
    ({candidate_id: Counter(record sha256)}, {node: [sha256 of every record mapped to it]}, query)."""
    covered = {m["candidate_id"]: Counter() for m in metadata}
    on_node = defaultdict(list)
    targets = {m["node_id"] for m in metadata}
    query = {}
    for alignment in IndexedGam(gam, index).fetch(targets, query):
        if alignment.mapping_quality <= min_mapq:
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


GAP = BASES["-"]


def expected_site_codes(alleles, reference, slots, label):
    """Channel-2 codes over the site columns for a row carrying `label` (independent of SiteLayout)."""
    if label == "OTHER":
        return None
    span = [BASES.get(b, 5) for b in reference]
    inserted = [GAP] * slots
    if label != "REF":
        allele = alleles[label]
        if allele["event_type"] == "INS":
            inserted = [BASES.get(b, 5) for b in allele["alt"]] + [GAP] * (slots - len(allele["alt"]))
        elif allele["event_type"] == "DEL":
            span = [GAP] * len(allele["ref"]) + span[len(allele["ref"]):]
        else:
            span = [BASES.get(b, 5) for b in allele["alt"]]
    return inserted + span


def validate(folder, gam, index, graph_index):
    folder = Path(folder)
    manifest = json.loads((folder / "manifest.json").read_text())
    metadata = [json.loads(s) for s in (folder / "variant_summary.ndjson").read_text().splitlines()]
    check(manifest.get("debug_rows"), "Validation requires tensors built with --debug-rows")
    check(manifest.get("tensor_format_version") == FORMAT_VERSION, "Unexpected tensor format version")
    check(manifest.get("tensor_storage_version") == STORAGE_VERSION, "Unexpected storage version")
    width = manifest["shape"][2]
    anchor = width // 2
    graph_nodes = {m["node_id"] for m in metadata} | {
        n for m in metadata for a in (m.get("alleles") or [m]) for n, _, _ in a.get("path") or []} | {
        c["node_id"] for m in metadata for row in m["rows"] for c in row["columns"] if c}
    with GraphIndex(graph_index) as lookup:
        records = lookup.get_nodes(graph_nodes)
    sequences = {n: r["sequence"] for n, r in records.items()}
    counts = {n: r["distinct_path_count"] for n, r in records.items()}
    parameters = manifest["parameters"]
    site_alleles = {}
    for m in metadata:  # allele units have no alleles[]: the record itself is its only allele
        listed = m.get("alleles") or [dict(m, label="A1")]
        site_alleles[m["candidate_id"]] = {a["label"]: dict(a, node_id=m["node_id"]) for a in listed}
    everything = [a for alleles in site_alleles.values() for a in alleles.values()]
    covered, on_node, query = independent_coverage(gam, index, everything, sequences, parameters["min_mapq"])
    cap = parameters.get("max_node_reads", 0)
    kept = {node: cap_records(digests, cap) for node, digests in on_node.items()}
    capped_nodes = sum(k is not None for k in kept.values())
    results = []
    for meta in metadata:
        x = np.load(folder / f"shard_{meta['shard_index']:05d}_data.npy", mmap_mode="r")[meta["index_within_shard"]]
        label = meta["candidate_id"]
        alleles = site_alleles[label]
        n = meta["selected_alignments"]
        for a in alleles.values():  # records beyond the cap were never considered
            if kept.get(meta["node_id"]) is not None:
                covered[a["candidate_id"]] = covered[a["candidate_id"]] & kept[meta["node_id"]]
        site = Counter()
        for a in alleles.values():
            site |= covered[a["candidate_id"]]
        check(meta["window_encoding_version"] == WINDOW_ENCODING_VERSION, f"{label}: window encoding")
        check(meta["row_selection_version"] == ROW_SELECTION_VERSION, f"{label}: row selection version")
        check(list(x.shape) == manifest["shape"] and str(x.dtype) == manifest["dtype"], f"{label}: tensor shape/dtype")
        check(len(meta["rows"]) == n, f"{label}: debug row count")
        check(not x[:, n:].any(), f"{label}: unused row padding")
        check(meta["coverage"] == sum(meta[f"{c}_count"] for c in ("alt", "ref", "other")), f"{label}: counts")
        check(meta["af"] == meta["alt_count"] / meta["coverage"], f"{label}: AF")
        check(sum(covered[label].values()) == meta["coverage"], f"{label}: independent representative coverage")
        check(sum(site.values()) == meta["site_coverage"], f"{label}: independent site coverage")
        selected = Counter(row["record_sha256"] for row in meta["rows"])
        check(not selected - site, f"{label}: selected source records")
        if n == meta["site_coverage"]:
            check(selected == site, f"{label}: full source record multiset")
        check(Counter(row["site_allele"] for row in meta["rows"]) == Counter(meta["selected_counts"]),
              f"{label}: selected allele counts")
        check(bool(np.all(x[1, :n][x[0, :n] == 6] == -1)), f"{label}: gap quality")
        # Site layout from the alleles and the graph sequence alone.
        slots = max((len(a["alt"]) for a in alleles.values() if a["event_type"] == "INS"), default=0)
        spanned = [a for a in alleles.values() if a["event_type"] != "INS"]
        span = max((len(a["ref"]) for a in spanned), default=0)
        def graph_ref(a):  # the graph bases a SNV/DEL allele covers, over every node of its path
            path = a.get("path") or []
            first = len(a["ref"]) - sum(e - s for _, s, e in path)
            return (sequences[meta["node_id"]][a["start"]:a["start"] + first]
                    + "".join(sequences[n][s:e] for n, s, e in path))
        check(all(a["ref"] == graph_ref(a) for a in spanned), f"{label}: allele REF vs graph")
        reference = max(spanned, key=lambda a: len(a["ref"]))["ref"] if spanned else ""
        check(meta["site_layout"] == dict(insertion_slots=slots, span=span, reference=reference), f"{label}: layout")
        lo, hi = meta["candidate_columns"]
        check(lo == anchor == meta["anchor_column"] and hi == min(width, lo + slots + span), f"{label}: site columns")
        checked = 0
        for ri, row in enumerate(meta["rows"]):
            strand = STRAND["reverse"] if row["reversed_for_candidate"] else STRAND["forward"]
            codes = expected_site_codes(alleles, reference, slots, row["site_allele"])
            for ci, col in enumerate(row["columns"]):
                if col is None:
                    check(not x[:, ri, ci].any(), f"{label}: all-zero missing evidence")
                    continue
                expected = codes[ci - lo] if codes is not None and lo <= ci < hi else 0
                check(int(x[2, ri, ci]) == expected, f"{label}: site allele channel at row {ri}, column {ci}")
                check(int(x[7, ri, ci]) == strand, f"{label}: strand channel at row {ri}, column {ci}")
                if ci == anchor and not slots and "offset" in col:
                    check(col["mapping_index"] == row["anchor_mapping_index"] and col["offset"] == meta["start"],
                          f"{label}: anchor coordinate")
                check(int(x[6, ri, ci]) == encode_count(counts[col["node_id"]]),
                      f"{label}: path count at row {ri}, column {ci}")
                if "offset" in col:
                    base = sequences[col["node_id"]][col["offset"]]
                    if col["reverse"]:
                        base = rc(base)
                    check(int(x[5, ri, ci]) == BASES.get(base, 5), f"{label}: graph base at row {ri}, column {ci}")
                    checked += 1
            mismatch = int(np.count_nonzero((x[0, ri] != 0) & (x[5, ri] != 0) & (x[0, ri] != x[5, ri]) & (x[4, ri] != 6)))
            check(mismatch == row["window_mismatch_bp"] == meta["window_mismatch_bp"][ri], f"{label}: visible edit bp")
        # Blocks A1..Ak, REF, OTHER; uniform sampling over the blocks ordered by record hash.
        blocks = list(alleles) + ["REF", "OTHER"]
        audit = meta["selection_audit"]
        check(Counter(a["record_sha256"] for a in audit) == site, f"{label}: all preselection records")
        ranked = sorted(audit, key=lambda a: (blocks.index(a["site_allele"]), a["record_sha256"],
                                              a["anchor_mapping_index"]))
        k = min(manifest["shape"][1], len(ranked))
        indices = ([len(ranked) // 2] if k == 1 else
                   [i * (len(ranked) - 1) // (k - 1) for i in range(k)] if k else [])
        check(sorted(meta["selected_ranks"]) == indices and n == k, f"{label}: uniform sampling ranks")
        for ri, (rank, row) in enumerate(zip(meta["selected_ranks"], meta["rows"])):
            a = ranked[rank]
            check(row["record_sha256"] == a["record_sha256"] and row["anchor_mapping_index"] == a["anchor_mapping_index"]
                  and row["site_allele"] == a["site_allele"], f"{label}: sampled record at row {ri}")
        next_row = 0
        for group in meta["row_groups"]:
            check(group["start_row"] == next_row and group["end_row"] > next_row, f"{label}: block boundaries")
            check(all(r["site_allele"] == group["allele"] for r in meta["rows"][next_row:group["end_row"]]),
                  f"{label}: block membership")
            next_row = group["end_row"]
        check(next_row == n and [g["allele"] for g in meta["row_groups"]] == [
            b for b in blocks if meta["selected_counts"].get(b)], f"{label}: block order")
        results.append(dict(candidate_id=label, event_type=meta["event_type"], coverage=meta["coverage"],
                            site_coverage=meta["site_coverage"], alleles=len(alleles), selected_alignments=n,
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

"""Byte and normalized comparison of two run roots, task directories or build directories.

Standard library only (imports no pipeline module), so it fingerprints v2 and v3 outputs alike:

    python -m indexed_gam_pipeline_v3.tools.compare_runs A B [--mask DOTTED.KEY ...] [--report FILE]

exits 1 on any difference. A file present in only one tree is a difference. File classes:

    ignored       anything under source/, logs/ or incomplete/; memory.ndjson, queue_status*.json,
                  finalize_report.json, run.sh, *.tmp
    batch_timing.ndjson
                  every row without elapsed_seconds and cumulative_stage_seconds, row order kept
    manifest.json, labels.manifest.json
                  the whole JSON, key order kept, the anchor path replaced by <ROOT>, without timing,
                  graph_index_performance, decoder, created, arguments.decoder, sources.outputs_json_sha256
    outputs.json, outputs.pre_merge.json
                  <ROOT>-normalized, without merge.created, merge.copy_seconds, merge.verify_seconds
    status.json   projected to status, tasks, processes, tensors, tensors_by_type, merged, merge_layout, labeled
    config.json   projected to tensors, chromosome_selection, tasks, processes, schedule, parts,
                  supplement.rounds, variant_outputs, postprocess and builder without candidate_unit
                  (<ROOT>-normalized)
    everything else (shards, labels, summaries, audit streams, displaced/target/part node lists,
                  validation reports, truth tables, recall reports, ...): raw bytes

--mask removes further dotted keys from every JSON file of the manifest, outputs, status and config
classes (e.g. --mask arguments.output shared_output).
"""
import argparse
import fnmatch
import hashlib
import json
from pathlib import Path
import sys

IGNORED_DIRS = ("source", "logs", "incomplete")
IGNORED_NAMES = ("memory.ndjson", "queue_status*.json", "finalize_report.json", "run.sh", "*.tmp")
TIMING_ROW = ("elapsed_seconds", "cumulative_stage_seconds")
MANIFEST_MASKS = ("timing", "graph_index_performance", "decoder", "created", "arguments.decoder",
                  "sources.outputs_json_sha256")
OUTPUTS_MASKS = ("merge.created", "merge.copy_seconds", "merge.verify_seconds")
STATUS_KEYS = ("status", "tasks", "processes", "tensors", "tensors_by_type", "merged", "merge_layout", "labeled")
CONFIG_KEYS = ("tensors", "chromosome_selection", "tasks", "processes", "schedule", "parts", "supplement",
               "variant_outputs", "postprocess", "builder")
ANCHOR = "<ROOT>"


def file_class(relative):
    """'ignored', 'timing', 'manifest', 'outputs', 'status', 'config' or 'raw' for a tree-relative path."""
    parts = Path(relative).parts
    name = parts[-1]
    if any(p in IGNORED_DIRS for p in parts[:-1]) or any(fnmatch.fnmatch(name, g) for g in IGNORED_NAMES):
        return "ignored"
    return {"batch_timing.ndjson": "timing", "manifest.json": "manifest", "labels.manifest.json": "manifest",
            "outputs.json": "outputs", "outputs.pre_merge.json": "outputs", "status.json": "status",
            "config.json": "config"}.get(name, "raw")


def _unmask(value, dotted):
    head, _, rest = dotted.partition(".")
    if isinstance(value, dict) and head in value:
        if rest:
            _unmask(value[head], rest)
        else:
            del value[head]


def _anchored(value, anchors):
    if isinstance(value, str):
        for anchor in anchors:
            value = value.replace(anchor, ANCHOR)
        return value
    if isinstance(value, list):
        return [_anchored(v, anchors) for v in value]
    if isinstance(value, dict):
        return {k: _anchored(v, anchors) for k, v in value.items()}
    return value


def normalized_json(path, kind, anchors=(), masks=()):
    """The comparable JSON value of a manifest/outputs/status/config file (key order kept)."""
    value = json.loads(Path(path).read_text())
    if kind == "status":
        value = {k: v for k, v in value.items() if k in STATUS_KEYS}
    elif kind == "config":
        value = {k: v for k, v in value.items() if k in CONFIG_KEYS}
        if isinstance(value.get("supplement"), dict):
            value["supplement"] = {"rounds": value["supplement"].get("rounds")}
        if isinstance(value.get("builder"), dict):
            value["builder"] = {k: v for k, v in value["builder"].items() if k != "candidate_unit"}
    defaults = dict(manifest=MANIFEST_MASKS, outputs=OUTPUTS_MASKS).get(kind, ())
    for dotted in (*defaults, *masks):
        _unmask(value, dotted)
    return _anchored(value, anchors)


def normalized_bytes(path, kind, anchors=(), masks=()):
    if kind == "raw":
        return None
    if kind == "timing":
        rows = []
        for line in Path(path).read_text().splitlines():
            row = json.loads(line)
            rows.append(json.dumps({k: v for k, v in row.items() if k not in TIMING_ROW}))
        return ("\n".join(rows) + "\n").encode()
    return json.dumps(normalized_json(path, kind, anchors, masks)).encode()


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _anchors(tree, anchor):
    """The anchor as given and resolved (longest first), so either spelling becomes <ROOT>."""
    anchor = Path(anchor if anchor is not None else tree)
    return tuple(sorted({str(anchor), str(anchor.resolve())}, key=len, reverse=True))


def fingerprint(tree, anchor=None, masks=()):
    """{relative path: sha256 of the normalized bytes} of every compared file under `tree`.

    `anchor` (default: the tree itself) is the directory whose path is replaced by <ROOT> in JSON files.
    """
    tree = Path(tree)
    anchors = _anchors(tree, anchor)
    result = {}
    for path in sorted(tree.rglob("*")):
        relative = path.relative_to(tree).as_posix()
        kind = file_class(relative)
        if kind == "ignored" or not path.is_file():
            continue
        data = normalized_bytes(path, kind, anchors, masks)
        result[relative] = _sha256_file(path) if data is None else hashlib.sha256(data).hexdigest()
    return result


def differences(a, b):
    """Sorted descriptions of the differences between two fingerprints."""
    found = [f"only in A: {k}" for k in sorted(a.keys() - b.keys())]
    found += [f"only in B: {k}" for k in sorted(b.keys() - a.keys())]
    found += [f"differs: {k}" for k in sorted(a.keys() & b.keys()) if a[k] != b[k]]
    return found


def _detail(a, b, relative, anchors_a, anchors_b, masks):
    """For a differing JSON file: its top-level keys whose normalized values differ."""
    kind = file_class(relative)
    if kind not in ("manifest", "outputs", "status", "config"):
        return ""
    x = normalized_json(Path(a) / relative, kind, anchors_a, masks)
    y = normalized_json(Path(b) / relative, kind, anchors_b, masks)
    keys = [k for k in dict.fromkeys([*x, *y]) if x.get(k, KeyError) != y.get(k, KeyError)]
    if not keys and list(x) != list(y):
        return " (key order)"
    return " (keys: " + ", ".join(keys) + ")" if keys else ""


def compare(a, b, masks=(), anchor_a=None, anchor_b=None, fingerprints=None):
    """Differences between trees `a` and `b` (empty list = identical after normalization).

    `fingerprints` = (fingerprint of a, fingerprint of b) when already computed.
    """
    fa, fb = fingerprints or (fingerprint(a, anchor_a, masks), fingerprint(b, anchor_b, masks))
    anchors_a, anchors_b = _anchors(a, anchor_a), _anchors(b, anchor_b)
    return [f + _detail(a, b, f.split(": ", 1)[1], anchors_a, anchors_b, masks) if f.startswith("differs") else f
            for f in differences(fa, fb)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("a")
    parser.add_argument("b")
    parser.add_argument("--mask", nargs="*", default=[], help="further dotted JSON keys to ignore")
    parser.add_argument("--report", help="also write the result as JSON")
    args = parser.parse_args(argv)
    fingerprints = fingerprint(args.a, masks=args.mask), fingerprint(args.b, masks=args.mask)
    compared = len(fingerprints[0])
    found = compare(args.a, args.b, args.mask, fingerprints=fingerprints)
    for line in found:
        print(line)
    print(f"{len(found)} differences over {compared} compared files of A")
    if args.report:
        Path(args.report).write_text(json.dumps(dict(a=str(Path(args.a).resolve()), b=str(Path(args.b).resolve()),
                                                     masks=args.mask, compared_files=compared,
                                                     differences=found), indent=2) + "\n")
    sys.exit(1 if found else 0)


if __name__ == "__main__":
    main()

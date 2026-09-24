"""Merge a finished run's task outputs into per-chromosome shards of a fixed size.

Input: <root>/outputs.json of a completed `orchestrate run` (every typed task directory
task_NNNN/ with shard_XXXXX_data.npy, variant_summary.ndjson, manifest.json and the audit
streams). Output, per tensor kind (SNV, INDEL, or ALL):

    <tensors>/<kind>/<chrom>_shard_NNNNN_data.npy      chr1..chr22; `shard_size` tensors, last one shorter
    <tensors>/<kind>/<chrom>_variant_summary.ndjson    source records in source order; shard_index /
                                                       index_within_shard rewritten, plus chrom, shard_file,
                                                       source_task, source_shard_index, source_index_within_shard,
                                                       grch38 (reference_path.linear of the candidate; null off a
                                                       unique GRCh38 node)
    <tensors>/<kind>/filtered_candidates.ndjson        task streams concatenated in task order (all chromosomes)
    <tensors>/<kind>/unsupported_events.ndjson
    <tensors>/<kind>/manifest.json                     layout, per-chromosome shards with SHA-256, provenance
    <tensors>/<kind>/validation_report.json
    <tensors>/non_autosomal/<kind>/...                 same layout for chrX, chrY, chrM, chrEBV, unplaced
    <root>/batch_timing.ndjson                         shared batch timings of every task (split runs)

Records keep their source order (task, shard, index), so shards stay node-sorted. The
chromosome of a tensor is the block of its node (chr_index). Everything is written to a
hidden `.merging/` directory first and published only after verification:
  * every output shard's data bytes hash (SHA-256, re-read from disk) equal the hash of the
    source tensors that were copied into it, in order;
  * every output summary record equals its source record apart from the position fields;
  * counts per chromosome add up to outputs.json; shard headers, lengths and file sizes agree;
  * an independent spot check reloads random tensors from the source task shards.
Sources are deleted only when `keep_sources` is False and every kind verified.
"""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import re
import shutil
import time

import numpy as np

from indexed_gam_pipeline_v2.common import read_json, sha256_file, write_json
from indexed_gam_pipeline_v2.tensor_postprocessing.chr_index import ChrIndex

LAYOUT = "chromosome-shards-v1"
DEFAULT_SHARD_SIZE = 32768
POSITION = ("shard_index", "index_within_shard")
ADDED = ("chrom", "shard_file", "source_task", "source_shard_index", "source_index_within_shard", "grch38")
# Must agree across all task manifests of one kind; copied into the merged manifest.
SHARED_KEYS = ("schema_version", "tensor_format_version", "tensor_storage_version", "shape", "dtype", "channels",
               "encodings", "row_selection_version", "row_order", "window_encoding_version", "parameters",
               "sample_unit", "site_definition", "read_cap", "debug_rows")
COPIED_KEYS = SHARED_KEYS + ("graph_index", "gai_version")
NODE = re.compile(rb'"node_id": (\d+)')
AUDIT_STREAMS = ("filtered_candidates", "unsupported_events")


def shard_name(group, index):
    return f"{group}_shard_{index:05d}_data.npy"


def summary_name(group):
    return f"{group}_variant_summary.ndjson"


def stripped(record, extra=()):
    return json.dumps({k: v for k, v in record.items() if k not in POSITION and k not in extra}, sort_keys=True)


def data_sha256(path):
    """SHA-256 of an .npy file's data region (header excluded), read from disk."""
    array = np.load(path, mmap_mode="r")
    offset, nbytes = array.offset, array.nbytes
    del array
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        stream.seek(offset)
        remaining = nbytes
        while remaining:
            chunk = stream.read(min(64 << 20, remaining))
            if not chunk:
                raise ValueError(f"Truncated shard: {path}")
            digest.update(chunk)
            remaining -= len(chunk)
    return digest.hexdigest()


class GroupWriter:
    """Fills <group>_shard_*.npy files of exactly `shard_size` tensors (last one shorter) and the summary."""

    def __init__(self, directory, group, total, shard_size, shape, linear=None):
        self.directory, self.group, self.shape, self.linear = Path(directory), group, tuple(shape), linear
        self.sizes = [shard_size] * (total // shard_size) + ([total % shard_size] if total % shard_size else [])
        self.shard, self.position, self.array, self.hasher = -1, 0, None, None
        self.summary = (self.directory / summary_name(group)).open("w")
        self.summary_hash = hashlib.sha256()
        self.shards = []

    def _next(self):
        self._close_shard()
        self.shard += 1
        if self.shard >= len(self.sizes):
            raise ValueError(f"More tensors than planned for {self.group}")
        self.array = np.lib.format.open_memmap(self.directory / shard_name(self.group, self.shard), mode="w+",
                                               dtype=np.int8, shape=(self.sizes[self.shard],) + self.shape)
        self.position, self.hasher = 0, hashlib.sha256()

    def _close_shard(self):
        if self.array is not None:
            if self.position != len(self.array):
                raise ValueError(f"Shard {self.group}/{self.shard} was not filled")
            self.array.flush()
            self.array = None
            self.shards.append(dict(file=shard_name(self.group, self.shard), tensors=self.sizes[self.shard],
                                    sha256=self.hasher.hexdigest()))

    def write(self, block, records, source_task, source_shard, first_index):
        """Append source tensors block[i] with their summary records (dicts) in order."""
        done = 0
        while done < len(block):
            if self.array is None or self.position == len(self.array):
                self._next()
            take = min(len(block) - done, len(self.array) - self.position)
            chunk = np.ascontiguousarray(block[done:done + take])
            self.array[self.position:self.position + take] = chunk
            self.hasher.update(memoryview(chunk).cast("B"))
            for k in range(take):
                source = records[done + k]
                self.summary_hash.update((stripped(source) + "\n").encode())
                record = dict(source, chrom=self.group, shard_file=shard_name(self.group, self.shard),
                              shard_index=self.shard, index_within_shard=self.position + k,
                              source_task=source_task, source_shard_index=source_shard,
                              source_index_within_shard=first_index + done + k)
                if self.linear is not None:
                    record["grch38"] = self.linear(source["node_id"], source["start"], source["ref"], source["alt"],
                                                   source["event_type"], source.get("path"))
                self.summary.write(json.dumps(record) + "\n")
            self.position += take
            done += take

    def close(self):
        self._close_shard()
        self.summary.close()
        if len(self.shards) != len(self.sizes):
            raise ValueError(f"{self.group}: wrote {len(self.shards)} of {len(self.sizes)} planned shards")
        return dict(tensors=sum(self.sizes), shards=self.shards, summary=summary_name(self.group),
                    summary_sha256_stripped=self.summary_hash.hexdigest())


def plan_kind(sources, index):
    """Reference manifest, per-source group arrays and per-group totals; checks manifests agree."""
    reference, groups, totals = None, [], Counter()
    for source in sources:
        folder = Path(source["path"])
        manifest = read_json(folder / "manifest.json")
        if manifest.get("status") != "complete":
            raise ValueError(f"Incomplete task output: {folder}")
        if reference is None:
            reference = manifest
        mismatch = [k for k in SHARED_KEYS if manifest.get(k) != reference.get(k)]
        if mismatch:
            raise ValueError(f"{folder}: manifest differs from task {sources[0]['task']} in {mismatch}")
        if manifest["tensors"] != source["tensors"] or manifest["shards"] != source["shards"]:
            raise ValueError(f"{folder}: manifest counts disagree with outputs.json")
        with (folder / "variant_summary.ndjson").open("rb") as stream:
            nodes = np.array([int(NODE.search(line).group(1)) for line in stream], dtype=np.int64)
        if len(nodes) != source["tensors"]:
            raise ValueError(f"{folder}: summary has {len(nodes)} records, manifest {source['tensors']}")
        k = index.lookup(nodes)
        if (k < 0).any():
            raise ValueError(f"{folder}: nodes outside every chromosome block, e.g. {nodes[k < 0][:5].tolist()}")
        groups.append(k)
        totals.update(Counter(k.tolist()))
    if reference is None:
        raise ValueError("No task outputs to merge")
    if reference["dtype"] != "int8":
        raise ValueError("Only int8 tensors are supported")
    return reference, groups, totals


def write_kind(kind, sources, index, shard_size, directories, reference_path=None):
    """Copy every source tensor of one kind into its group's writer; returns per-group results."""
    reference, groups, totals = plan_kind(sources, index)
    shape = reference["shape"]
    linear = None
    if reference_path is not None:
        from indexed_gam_pipeline_v2.tensor_postprocessing.reference_path import ReferencePath
        linear = ReferencePath(reference_path).linear
    writers = {k: GroupWriter(directories[index.dataset[k]], index.names[k], n, shard_size, shape, linear)
               for k, n in sorted(totals.items())}
    for source, ks in zip(sources, groups):
        folder = Path(source["path"])
        with (folder / "variant_summary.ndjson").open() as stream:
            records = [json.loads(line) for line in stream]
        position = 0
        for shard in range(source["shards"]):
            data = np.load(folder / f"shard_{shard:05d}_data.npy", mmap_mode="r")
            if data.dtype != np.int8 or list(data.shape[1:]) != shape:
                raise ValueError(f"{folder}: shard {shard} shape/dtype")
            n = len(data)
            chunk_records, chunk_groups = records[position:position + n], ks[position:position + n]
            if len(chunk_records) != n or any(r["shard_index"] != shard or r["index_within_shard"] != i
                                              for i, r in enumerate(chunk_records)):
                raise ValueError(f"{folder}: summary positions disagree with shard {shard}")
            cuts = [0] + (np.flatnonzero(np.diff(chunk_groups)) + 1).tolist() + [n]
            for a, b in zip(cuts, cuts[1:]):
                writers[int(chunk_groups[a])].write(data[a:b], chunk_records[a:b], source["task"], shard, a)
            position += n
            del data
        if position != len(records):
            raise ValueError(f"{folder}: summary has records beyond the last shard")
    results = {index.names[k]: dict(w.close(), dataset=index.dataset[k]) for k, w in writers.items()}
    for stream in AUDIT_STREAMS:
        with (directories["autosome"] / f"{stream}.ndjson").open("wb") as target:
            for source in sources:
                path = Path(source["path"]) / f"{stream}.ndjson"
                if path.exists():
                    with path.open("rb") as handle:
                        shutil.copyfileobj(handle, target, 16 << 20)
    return reference, results


def verify_summary(directory, group, result):
    """Re-read one output summary: positions, file names, chromosome, and the stripped-record hash."""
    digest, counts = hashlib.sha256(), Counter()
    expected = {s["file"]: s["tensors"] for s in result["shards"]}
    last = (-1, -1)
    with (Path(directory) / result["summary"]).open() as stream:
        for line in stream:
            r = json.loads(line)
            position = (r["shard_index"], r["index_within_shard"])
            follows = position == (last[0], last[1] + 1) or position == (last[0] + 1, 0)
            if (not follows or r["chrom"] != group or r["shard_file"] != shard_name(group, r["shard_index"])
                    or r["shard_file"] not in expected):
                raise ValueError(f"{group}: summary position/name mismatch at {position}")
            last = position
            counts[r["shard_file"]] += 1
            digest.update((stripped(r, ADDED) + "\n").encode())
    if dict(counts) != expected:
        raise ValueError(f"{group}: summary counts {dict(counts)} != shards {expected}")
    if digest.hexdigest() != result["summary_sha256_stripped"]:
        raise ValueError(f"{group}: summary records differ from the source records")


def verify_shard(path, tensors, shape, sha256):
    array = np.load(path, mmap_mode="r")
    ok = (array.dtype == np.int8 and array.shape == (tensors,) + tuple(shape)
          and Path(path).stat().st_size == array.offset + array.nbytes)
    del array
    if not ok:
        raise ValueError(f"Shard header/size mismatch: {path}")
    if data_sha256(path) != sha256:
        raise ValueError(f"Shard bytes differ from their source tensors: {path}")
    return str(path)


def merge(root, chr_index, shard_size=DEFAULT_SHARD_SIZE, keep_sources=False, workers=8, spots=200, seed=20260923,
          reference_path=None):
    """Merge every kind listed in <root>/outputs.json; returns the run-level merge report."""
    started = time.perf_counter()
    root = Path(root)
    outputs = read_json(root / "outputs.json")
    if "merge" in outputs:
        raise ValueError(f"{root} was already merged ({outputs['merge']['layout']})")
    index = ChrIndex(chr_index)
    tensors = Path(outputs["tensors_dir"])
    sources_of = {kind: sorted(items, key=lambda s: s["task"]) for kind, items in outputs["outputs"].items()}
    plans = {}
    for kind in sources_of:
        final = dict(autosome=tensors / kind, non_autosomal=tensors / "non_autosomal" / kind)
        work = {d: p / ".merging" for d, p in final.items()}
        for d, p in final.items():
            if (p / "manifest.json").exists() and read_json(p / "manifest.json").get("layout") == LAYOUT:
                raise ValueError(f"{p} already holds a merged layout")
            clashes = list(p.glob("*_shard_*_data.npy")) + list(p.glob("*_variant_summary.ndjson"))
            if clashes:
                raise ValueError(f"{p} already has merged-style files, e.g. {clashes[0]}")
            if work[d].exists():
                shutil.rmtree(work[d])  # leftovers of an interrupted merge are ours to discard
            work[d].mkdir(parents=True)
        plans[kind] = (final, work)
    # 1. Copy (kinds in parallel: they read and write disjoint files).
    with ProcessPoolExecutor(max_workers=max(1, min(workers, len(sources_of)))) as pool:
        futures = {kind: pool.submit(write_kind, kind, sources_of[kind], index, shard_size, plans[kind][1],
                                     reference_path)
                   for kind in sources_of}
        written = {kind: f.result() for kind, f in futures.items()}
    copy_seconds = time.perf_counter() - started
    # 2. Verify every shard from disk, every summary, the totals, and a spot check.
    t = time.perf_counter()
    jobs = []
    for kind, (reference, results) in written.items():
        work = plans[kind][1]
        for group, result in results.items():
            verify_summary(work[result["dataset"]], group, result)
            jobs += [(work[result["dataset"]] / s["file"], s["tensors"], reference["shape"], s["sha256"])
                     for s in result["shards"]]
        total = sum(r["tensors"] for r in results.values())
        if total != outputs["tensors_by_type"][kind] or total != sum(s["tensors"] for s in sources_of[kind]):
            raise ValueError(f"{kind}: merged {total} tensors, outputs.json lists {outputs['tensors_by_type'][kind]}")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        list(pool.map(verify_shard, *zip(*jobs)))
    spots_checked = {}
    for kind, (reference, results) in written.items():
        work = plans[kind][1]
        source_paths = {s["task"]: s["path"] for s in sources_of[kind]}
        spots_checked[kind] = spot_check_sources(work, results, source_paths, spots, seed)
    verify_seconds = time.perf_counter() - t
    # 3. Publish: move verified files into place, manifests last.
    created = datetime.now(timezone.utc).isoformat()
    report = dict(layout=LAYOUT, shard_size=shard_size, created=created, chr_index=dict(path=str(index.path),
                  sha256=index.sha256), kinds={}, copy_seconds=copy_seconds, verify_seconds=verify_seconds,
                  spot_checked=spots_checked, sources_kept=keep_sources,
                  reference_path=str(Path(reference_path).resolve()) if reference_path else None)
    for kind, (reference, results) in written.items():
        final, work = plans[kind]
        report["kinds"][kind] = {}
        for dataset in ("autosome", "non_autosomal"):
            groups = {g: r for g, r in results.items() if r["dataset"] == dataset}
            names = [s["file"] for r in groups.values() for s in r["shards"]] + [r["summary"] for r in groups.values()]
            if dataset == "autosome":
                names += [f"{s}.ndjson" for s in AUDIT_STREAMS]
            for name in names:
                target = final[dataset] / name
                if target.exists():
                    raise ValueError(f"Refusing to overwrite {target}")
                (work[dataset] / name).rename(target)
            manifest = dict(status="complete", layout=LAYOUT, kind=kind, dataset=dataset, shard_size=shard_size,
                            **{k: reference[k] for k in COPIED_KEYS if k in reference},
                            tensors=sum(r["tensors"] for r in groups.values()),
                            chromosomes={g: {k: v for k, v in r.items() if k != "dataset"}
                                         for g, r in sorted(groups.items())},
                            chr_index=report["chr_index"], created=created,
                            sources=dict(run_root=str(root), outputs_json_sha256=sha256_file(root / "outputs.json"),
                                         tasks=sources_of[kind]),
                            summary_fields_added=list(ADDED), summary_fields_rewritten=list(POSITION))
            if dataset == "autosome":
                manifest["audit_streams"] = {s: f"{s}.ndjson" for s in AUDIT_STREAMS}
            write_json(final[dataset] / "validation_report.json",
                       dict(passed=True, tensors=manifest["tensors"], shards=sum(len(r["shards"]) for r in groups.values()),
                            chromosomes={g: r["tensors"] for g, r in sorted(groups.items())},
                            checks=["shard header/length/size", "shard data SHA-256 from disk == copied source bytes",
                                    "summary positions, names, chromosome", "summary records == source records",
                                    "totals == outputs.json", "random tensors reloaded from source shards"]))
            write_json(final[dataset] / "manifest.json", manifest)
            work[dataset].rmdir()
            report["kinds"][kind][dataset] = dict(directory=str(final[dataset]), tensors=manifest["tensors"],
                                                  chromosomes={g: r["tensors"] for g, r in sorted(groups.items())})
    report["batch_timing"] = collect_batch_timing(root, tensors, sources_of)
    # 4. Run bookkeeping; sources are deleted last and only on request.
    backup = root / "outputs.pre_merge.json"
    if not backup.exists():
        shutil.copy2(root / "outputs.json", backup)
    if not keep_sources:
        delete_sources(tensors, sources_of)
    outputs["merge"] = report
    write_json(root / "outputs.json", outputs)
    status = read_json(root / "status.json")
    status.update(merged=True, merge_layout=LAYOUT, status="finalized" if not keep_sources else status["status"])
    write_json(root / "status.json", status)
    report["total_seconds"] = time.perf_counter() - started
    return report


def spot_check_sources(directories, results, source_paths, spots, seed):
    """Reload random merged tensors and compare them with the tensor at their recorded source position."""
    rng = random.Random(seed)
    checked = 0
    for group, result in sorted(results.items()):
        directory = directories[result["dataset"]]
        with (directory / result["summary"]).open() as stream:
            lines = stream.readlines()
        for line in rng.sample(lines, min(spots, len(lines))):
            r = json.loads(line)
            out = np.load(directory / r["shard_file"], mmap_mode="r")[r["index_within_shard"]]
            folder = Path(source_paths[r["source_task"]])
            source = np.load(folder / f"shard_{r['source_shard_index']:05d}_data.npy",
                             mmap_mode="r")[r["source_index_within_shard"]]
            if not np.array_equal(out, source):
                raise ValueError(f"Spot check failed: {group} {r['shard_file']}[{r['index_within_shard']}]")
            checked += 1
    return checked


def collect_batch_timing(root, tensors, sources_of):
    """Concatenate shared/task_*/batch_timing.ndjson (split runs) into <root>/batch_timing.ndjson."""
    tasks = sorted({s["task"] for items in sources_of.values() for s in items})
    target, lines = Path(root) / "batch_timing.ndjson", 0
    with target.open("w") as out:
        for task in tasks:
            for kind in ("shared", *sources_of):
                path = Path(tensors) / kind / f"task_{task:04d}" / "batch_timing.ndjson"
                if path.exists():
                    with path.open() as stream:
                        for line in stream:
                            out.write(json.dumps(dict(task=task, **json.loads(line))) + "\n")
                            lines += 1
                    break
    return dict(path=str(target), batches=lines)


def delete_sources(tensors, sources_of):
    """Remove the per-task directories (typed and shared) after a verified merge."""
    tasks = sorted({s["task"] for items in sources_of.values() for s in items})
    for kind, items in sources_of.items():
        for source in items:
            shutil.rmtree(source["path"])
    for task in tasks:
        shared = Path(tensors) / "shared" / f"task_{task:04d}"
        if shared.exists():
            shutil.rmtree(shared)

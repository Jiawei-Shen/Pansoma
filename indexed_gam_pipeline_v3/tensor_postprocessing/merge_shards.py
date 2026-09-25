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

The copy runs in parallel (`workers` processes): every (kind, chromosome) group is one job that
reads its records' rows with plain sequential file reads and appends them to its shards, and the
audit streams are copied in slices straight to their offsets in the concatenated file. The
bytes are those of a sequential copy.
"""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import random
import re
import shutil
import time

import numpy as np

from ..common import read_json, sha256_file, write_json
from .chr_index import ChrIndex

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


def npy_header(shape):
    """The .npy header bytes of an int8 C-order array (what np.lib.format.open_memmap writes)."""
    buffer = io.BytesIO()
    np.lib.format.write_array_header_1_0(buffer, dict(descr=np.lib.format.dtype_to_descr(np.dtype(np.int8)),
                                                      fortran_order=False, shape=tuple(shape)))
    return buffer.getvalue()


def read_rows(path, start, stop, shape):
    """Rows [start, stop) of an int8 .npy shard with one sequential read (no memory map)."""
    header = np.load(path, mmap_mode="r")  # parses the header only
    offset, row = header.offset, header.itemsize * int(np.prod(header.shape[1:], dtype=np.int64))
    if header.dtype != np.int8 or list(header.shape[1:]) != list(shape) or not 0 <= start <= stop <= len(header):
        raise ValueError(f"{path}: shape/dtype")
    del header
    rows = np.empty((stop - start,) + tuple(shape), dtype=np.int8)
    if not rows.size:
        return rows
    with open(path, "rb") as stream:
        stream.seek(offset + start * row)
        if stream.readinto(memoryview(rows).cast("B")) != rows.nbytes:
            raise ValueError(f"Truncated shard: {path}")
    return rows


class GroupWriter:
    """Fills <group>_shard_*.npy files of exactly `shard_size` tensors (last one shorter) and the summary."""

    def __init__(self, directory, group, total, shard_size, shape, linear=None):
        self.directory, self.group, self.shape, self.linear = Path(directory), group, tuple(shape), linear
        self.sizes = [shard_size] * (total // shard_size) + ([total % shard_size] if total % shard_size else [])
        self.shard, self.position, self.stream, self.hasher = -1, 0, None, None
        self.summary = (self.directory / summary_name(group)).open("w")
        self.summary_hash = hashlib.sha256()
        self.shards = []

    def _next(self):
        self._close_shard()
        self.shard += 1
        if self.shard >= len(self.sizes):
            raise ValueError(f"More tensors than planned for {self.group}")
        self.stream = open(self.directory / shard_name(self.group, self.shard), "wb")
        self.stream.write(npy_header((self.sizes[self.shard],) + self.shape))
        self.position, self.hasher = 0, hashlib.sha256()

    def _close_shard(self):
        if self.stream is not None:
            if self.position != self.sizes[self.shard]:
                raise ValueError(f"Shard {self.group}/{self.shard} was not filled")
            self.stream.close()
            self.stream = None
            self.shards.append(dict(file=shard_name(self.group, self.shard), tensors=self.sizes[self.shard],
                                    sha256=self.hasher.hexdigest()))

    def write(self, block, records, source_task, source_shard, first_index):
        """Append source tensors block[i] with their summary records (dicts) in order."""
        done = 0
        while done < len(block):
            if self.stream is None or self.position == self.sizes[self.shard]:
                self._next()
            take = min(len(block) - done, self.sizes[self.shard] - self.position)
            chunk = memoryview(np.ascontiguousarray(block[done:done + take])).cast("B")
            self.stream.write(chunk)
            self.hasher.update(chunk)
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


def scan_source(source, reference, first_task, index):
    """One task directory: complete, manifest agreeing with the reference; the group of every record."""
    folder = Path(source["path"])
    manifest = read_json(folder / "manifest.json")
    if manifest.get("status") != "complete":
        raise ValueError(f"Incomplete task output: {folder}")
    mismatch = [k for k in SHARED_KEYS if manifest.get(k) != reference.get(k)]
    if mismatch:
        raise ValueError(f"{folder}: manifest differs from task {first_task} in {mismatch}")
    if manifest["tensors"] != source["tensors"] or manifest["shards"] != source["shards"]:
        raise ValueError(f"{folder}: manifest counts disagree with outputs.json")
    with (folder / "variant_summary.ndjson").open("rb") as stream:
        nodes = np.array([int(NODE.search(line).group(1)) for line in stream], dtype=np.int64)
    if len(nodes) != source["tensors"]:
        raise ValueError(f"{folder}: summary has {len(nodes)} records, manifest {source['tensors']}")
    k = index.lookup(nodes)
    if (k < 0).any():
        raise ValueError(f"{folder}: nodes outside every chromosome block, e.g. {nodes[k < 0][:5].tolist()}")
    return k


def plan_kind(sources, index, pool=None):
    """Reference manifest, per-source group arrays and per-group totals; checks manifests agree."""
    if not sources:
        raise ValueError("No task outputs to merge")
    first = Path(sources[0]["path"])
    reference = read_json(first / "manifest.json")
    if reference.get("status") != "complete":
        raise ValueError(f"Incomplete task output: {first}")
    jobs = [(s, reference, sources[0]["task"], index) for s in sources]
    groups = list(pool.map(scan_source, *zip(*jobs), chunksize=max(1, len(jobs) // 64)) if pool
                  else map(scan_source, *zip(*jobs)))
    totals = Counter()
    for k in groups:
        totals.update(Counter(k.tolist()))
    if reference["dtype"] != "int8":
        raise ValueError("Only int8 tensors are supported")
    return reference, groups, totals


def write_group(directory, group, total, shard_size, shape, reference_path, parts):
    """Copy one group's records, `parts` = [(source, sorted record positions of the group)] in task
    order, into its shards and summary; the group's result."""
    linear = None
    if reference_path is not None:
        from .reference_path import ReferencePath
        linear = ReferencePath(reference_path).linear
    writer = GroupWriter(directory, group, total, shard_size, shape, linear)
    for source, positions in parts:
        folder = Path(source["path"])
        with (folder / "variant_summary.ndjson").open() as stream:
            lines = stream.readlines()
        start = 0
        for shard in range(source["shards"]):
            path = folder / f"shard_{shard:05d}_data.npy"
            header = np.load(path, mmap_mode="r")
            if header.dtype != np.int8 or list(header.shape[1:]) != list(shape):
                raise ValueError(f"{folder}: shard {shard} shape/dtype")
            n = len(header)
            del header
            lo, hi = np.searchsorted(positions, [start, start + n])
            local = positions[lo:hi] - start
            cuts = [0] + (np.flatnonzero(np.diff(local) != 1) + 1).tolist() + [len(local)]
            for a, b in zip(cuts, cuts[1:]):
                if a == b:
                    continue
                first, last = int(local[a]), int(local[b - 1]) + 1
                records = [json.loads(line) for line in lines[start + first:start + last]]
                if len(records) != last - first or any(r["shard_index"] != shard or r["index_within_shard"] != i
                                                      for i, r in zip(range(first, last), records)):
                    raise ValueError(f"{folder}: summary positions disagree with shard {shard}")
                writer.write(read_rows(path, first, last, shape), records, source["task"], shard, first)
            start += n
        if start != len(lines):
            raise ValueError(f"{folder}: summary has records beyond the last shard")
    return writer.close()


def copy_slice(target, offset, paths):
    """Write the files `paths`, concatenated, into `target` from byte `offset` on."""
    with open(target, "r+b") as out:
        out.seek(offset)
        for path in paths:
            with open(path, "rb") as handle:
                shutil.copyfileobj(handle, out, 16 << 20)


def stream_slices(target, paths, parts):
    """copy_slice jobs that together write the concatenation of `paths` into a presized `target`."""
    paths = [Path(p) for p in paths if Path(p).exists()]
    sizes = [p.stat().st_size for p in paths]
    with open(target, "wb") as out:
        out.truncate(sum(sizes))
    jobs, offset, step = [], 0, max(1, -(-sum(sizes) // max(1, parts)))
    batch, batch_offset, batch_bytes = [], 0, 0
    for path, size in zip(paths, sizes):
        if not batch:
            batch_offset = offset
        batch.append(path)
        batch_bytes += size
        offset += size
        if batch_bytes >= step:
            jobs.append((target, batch_offset, batch))
            batch, batch_bytes = [], 0
    if batch:
        jobs.append((target, batch_offset, batch))
    return jobs


def write_kinds(sources_of, index, shard_size, work, reference_path, pool, workers):
    """Every (kind, group) and every audit-stream slice as one pool job, largest first;
    {kind: (reference manifest, {group name: result})}."""
    planned = {kind: plan_kind(sources, index, pool) for kind, sources in sources_of.items()}
    jobs = []  # (bytes, kind, name, function, arguments)
    for kind, (reference, groups, totals) in planned.items():
        shape, row = reference["shape"], int(np.prod(reference["shape"]))
        for k, total in sorted(totals.items()):
            parts = [(s, np.flatnonzero(g == k)) for s, g in zip(sources_of[kind], groups) if (g == k).any()]
            jobs.append((total * row, kind, k, write_group,
                         (work[kind][index.dataset[k]], index.names[k], total, shard_size, shape, reference_path, parts)))
        for stream in AUDIT_STREAMS:
            paths = [Path(s["path"]) / f"{stream}.ndjson" for s in sources_of[kind]]
            for target, offset, batch in stream_slices(work[kind]["autosome"] / f"{stream}.ndjson", paths, workers):
                jobs.append((sum(p.stat().st_size for p in batch), kind, None, copy_slice, (target, offset, batch)))
    futures = [(kind, k, pool.submit(function, *arguments))
               for _, kind, k, function, arguments in sorted(jobs, key=lambda j: -j[0])]
    done = {(kind, k): future.result() for kind, k, future in futures}
    return {kind: (reference, {index.names[k]: dict(done[kind, k], dataset=index.dataset[k])
                               for k in sorted(totals)})
            for kind, (reference, _, totals) in planned.items()}


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
    # 1. Copy: every (kind, chromosome) group and every audit-stream slice is one parallel job.
    with ProcessPoolExecutor(max_workers=max(1, workers)) as pool:
        written = write_kinds(sources_of, index, shard_size, {kind: plans[kind][1] for kind in sources_of},
                              reference_path, pool, workers)
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

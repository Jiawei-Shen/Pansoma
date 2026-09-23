#!/usr/bin/env python3
"""Whole-genome runs on one Slurm node: many small tasks, one disposable builder process each.

    prepare   freeze a copy of this package, split the node list into N contiguous tasks,
              record input fingerprints and write config.json + run.sh under --root;
              with --node-stats (discovery's node_stats.json) also predict each task's cost
    run       execute the tasks (at most --processes at once, fail fast), most expensive
              predicted first when costs were recorded, else in index order; --resume skips
              tasks that already completed and validated
    task      run one task: builder subprocess, then validate its shards (used by `run`)

Layout under --root (bookkeeping) and --tensors (outputs; default <root>/tensors):
    <root>/source/                frozen package copy (hashes recorded in config.json)
    <root>/parts/nodes_NNNN.txt   node list of each task
    <root>/logs/task_NNNN.log     builder stdout/stderr; task_NNNN.resources.txt from /usr/bin/time -v
    <root>/queue_status.json, status.json, memory.ndjson, outputs.json
    <tensors>/shared/task_NNNN/   per-task shared manifest, batch timings, audit streams (split mode)
    <tensors>/SNV/task_NNNN/      SNV shards + summary + manifest                        (split mode)
    <tensors>/INDEL/task_NNNN/    INDEL shards + summary + manifest                      (split mode)
    <tensors>/ALL/task_NNNN/      everything in one directory                    (single-output mode)
"""
import argparse
from array import array
from collections import Counter
import json
import mmap
import os
import re
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from indexed_gam_pipeline_v2.build import KIND, SITE_UNIT
from indexed_gam_pipeline_v2.candidates import FORMAT_VERSION, STORAGE_VERSION
from indexed_gam_pipeline_v2.common import read_json, sha256_file, stamp, write_json
from indexed_gam_pipeline_v2.graph_index import GraphIndex

PACKAGE = Path(__file__).resolve().parent.name
# Every builder option is frozen into config.json and passed explicitly, so a prepared
# run never depends on the CLI defaults of the code that later executes it.
BUILDER_OPTIONS = ("rows", "width", "gam_cache_mb", "batch_nodes", "max_node_span", "max_batch_alignments",
                   "shard_size", "min_mapq", "min_af", "min_variants", "min_allele_bq", "max_indel_len", "chr",
                   "candidate_unit", "max_node_reads", "early_af_filter")


# --- prepare ------------------------------------------------------------------

def partition_nodes(source, folder, tasks):
    """Split a sorted unique node list into `tasks` contiguous, near-equal files."""
    nodes = np.fromfile(str(source), dtype=np.int64, sep="\n")
    if not len(nodes) or nodes[0] <= 0 or np.any(np.diff(nodes) <= 0):
        raise ValueError("Node list must be sorted, unique and positive")
    if not 1 <= tasks <= len(nodes):
        raise ValueError("Need at least one node per task")
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    parts = []
    for i, chunk in enumerate(np.array_split(nodes, tasks)):
        path = folder / f"nodes_{i:04d}.txt"
        path.write_text("".join(f"{n}\n" for n in chunk.tolist()))
        parts.append(dict(nodes_file=str(path), nodes=len(chunk), first_node=int(chunk[0]),
                          last_node=int(chunk[-1]), sha256=sha256_file(path)))
    return parts


NODE_STATS_RECORD = re.compile(rb'"(\d+)":\s*\{\s*"perfect":\s*(\d+),\s*"not_perfect":\s*(\d+)')


def node_costs(stats_path, nodes):
    """Discovery's not_perfect count of every node in `nodes` (sorted int64 array).

    The number of edited MAPQ-passing mappings on a node tracks how much decoding and
    candidate work the node costs. node_stats.json covers every observed node (several GB
    for a genome), so it is scanned as bytes instead of parsed into one dict.
    """
    found, counts = array("q"), array("q")
    with open(stats_path, "rb") as stream, mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ) as data:
        for match in NODE_STATS_RECORD.finditer(data):
            found.append(int(match.group(1)))
            counts.append(int(match.group(3)))
    found, counts = np.frombuffer(found, dtype=np.int64), np.frombuffer(counts, dtype=np.int64)
    order = np.argsort(found, kind="stable")
    found, counts = found[order], counts[order]
    at = np.minimum(np.searchsorted(found, nodes), max(len(found) - 1, 0))
    if not len(found) or np.any(found[at] != nodes):
        raise ValueError("node_stats.json lacks some target nodes: " + str(stats_path))
    return counts[at]


def prepare(args):
    root = Path(args.root).resolve()
    root.mkdir(parents=True, exist_ok=False)
    if (args.snv_min_af is None) != (args.indel_min_af is None):
        raise ValueError("Give both --snv-min-af and --indel-min-af, or neither")
    if not 1 <= args.processes <= args.tasks:
        raise ValueError("Require tasks >= processes >= 1")
    with GraphIndex(args.graph_index) as graph:
        graph_metadata = graph.metadata
    tensors = Path(args.tensors).resolve() if args.tensors else root / "tensors"
    if tensors.exists() and any(tensors.iterdir()):
        raise ValueError(f"Tensor output directory must be new or empty: {tensors}")
    package = Path(__file__).resolve().parent
    shutil.copytree(package, root / "source" / PACKAGE,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "tests"))
    parts = partition_nodes(args.nodes, root / "parts", args.tasks)
    schedule = dict(order="task index")
    if args.node_stats:
        nodes = np.fromfile(str(args.nodes), dtype=np.int64, sep="\n")
        costs = np.cumsum(node_costs(args.node_stats, nodes))
        end = 0
        for part in parts:
            begin, end = end, end + part["nodes"]
            part["predicted_cost"] = int(costs[end - 1] - (costs[begin - 1] if begin else 0))
        schedule = dict(order="predicted cost descending, then task index",
                        predicted_cost="sum of discovery not_perfect over the task's nodes",
                        node_stats=stamp(args.node_stats))
    index = args.index or str(args.gam) + ".gai"
    config = dict(
        created=time.strftime("%Y-%m-%dT%H:%M:%S"),
        python=sys.executable,
        tensors=str(tensors),
        inputs=dict(gam=stamp(args.gam), index=stamp(index),
                    graph_index=dict(stamp(args.graph_index), metadata=graph_metadata),
                    nodes=dict(stamp(args.nodes), sha256=sha256_file(args.nodes), count=sum(p["nodes"] for p in parts))),
        source_sha256={str(p.relative_to(root)): sha256_file(p) for p in sorted((root / "source").rglob("*")) if p.is_file()},
        tasks=args.tasks, processes=args.processes, schedule=schedule, parts=parts,
        builder={k: getattr(args, k) for k in BUILDER_OPTIONS},
        variant_outputs=(dict(SNV=args.snv_min_af, INDEL=args.indel_min_af) if args.snv_min_af is not None else None))
    write_json(root / "config.json", config)
    (root / "run.sh").write_text("\n".join([
        "#!/bin/bash", "set -euo pipefail",
        "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1",
        "export PYTHONUNBUFFERED=1",
        "cd " + shlex.quote(str(root / "source")),
        "exec " + shlex.quote(sys.executable) + f" -m {PACKAGE}.orchestrate run --root " + shlex.quote(str(root)) + ' "$@"', ""]))
    (root / "run.sh").chmod(0o755)
    print(json.dumps(dict(root=str(root), tensors=str(tensors), tasks=args.tasks, processes=args.processes,
                          nodes=config["inputs"]["nodes"]["count"], variant_outputs=config["variant_outputs"]), indent=2))
    return config


def verify(root, config):
    """Inputs, frozen source and partitions must be exactly what `prepare` recorded."""
    for key, expected in config["inputs"].items():
        current = stamp(expected["path"])
        if any(current[k] != expected[k] for k in ("size", "mtime_ns")):
            raise ValueError("Input changed since prepare: " + key)
    for relative, expected in config["source_sha256"].items():
        if sha256_file(root / relative) != expected:
            raise ValueError("Frozen source changed: " + relative)
    for part in config["parts"]:
        if sha256_file(part["nodes_file"]) != part["sha256"]:
            raise ValueError("Partition changed: " + part["nodes_file"])


# --- task -----------------------------------------------------------------------

def task_outputs(root, config, index):
    """{output type: directory} for one task, all under config['tensors'].

    Split mode: 'shared' (coordinator files) plus one directory per variant type.
    Single-output mode: just 'ALL', which also holds the coordinator files.
    """
    base, name = Path(config["tensors"]), f"task_{index:04d}"
    if config["variant_outputs"]:
        return {"shared": base / "shared" / name, **{kind: base / kind / name for kind in config["variant_outputs"]}}
    return {"ALL": base / "ALL" / name}


def build_command(root, config, index):
    b = config["builder"]
    outputs = task_outputs(root, config, index)
    command = [config["python"], "-m", f"{PACKAGE}.run", "build",
               "--gam", config["inputs"]["gam"]["path"], "--index", config["inputs"]["index"]["path"],
               "--graph-index", config["inputs"]["graph_index"]["path"],
               "--nodes", config["parts"][index]["nodes_file"],
               "--output", str(outputs.get("shared", outputs.get("ALL"))),
               "--variant-type", "all"]
    for key in BUILDER_OPTIONS:
        if key == "chr":
            if b["chr"]:
                command += ["--chr", b["chr"]]
        elif key == "early_af_filter":
            command += ["--early-af-filter" if b[key] else "--no-early-af-filter"]
        else:
            command += ["--" + key.replace("_", "-"), str(b[key])]
    if config["variant_outputs"]:
        for kind, af in config["variant_outputs"].items():
            flag = kind.lower()
            command += [f"--{flag}-output", str(outputs[kind]), f"--{flag}-min-af", str(af)]
    return command


def validate_shards(folder, shard_size):
    """Every shard file, its header and the summary must agree with the manifest."""
    folder = Path(folder)
    manifest = read_json(folder / "manifest.json")
    if manifest["status"] != "complete" or manifest["tensor_format_version"] != FORMAT_VERSION:
        raise ValueError("Incomplete or unexpected tensor format: " + str(folder))
    if manifest["dtype"] != "int8" or manifest["tensor_storage_version"] != STORAGE_VERSION:
        raise ValueError("Unexpected storage encoding: " + str(folder))
    shape = tuple(manifest["shape"])
    expected = manifest["shards"]
    if len(list(folder.glob("shard_*_data.npy"))) != expected:
        raise ValueError("Shard file count mismatch: " + str(folder))
    sizes = []
    for i in range(expected):
        path = folder / f"shard_{i:05d}_data.npy"
        x = np.load(path, mmap_mode="r", allow_pickle=False)
        if x.ndim != 4 or x.shape[1:] != shape or x.dtype != np.int8:
            raise ValueError(f"Unexpected shard shape/dtype: {path}")
        if not 0 < len(x) <= shard_size or (i < expected - 1 and len(x) != shard_size):
            raise ValueError(f"Invalid shard length: {path}")
        if path.stat().st_size != x.offset + x.nbytes:
            raise ValueError(f"Truncated or overlong shard: {path}")
        sizes.append(len(x))
        del x
    seen = [0] * expected
    events = Counter()
    site_mode = manifest.get("sample_unit") == SITE_UNIT
    sites = set()
    with (folder / "variant_summary.ndjson").open() as stream:
        for line in stream:
            m = json.loads(line)
            shard = m["shard_index"]
            if not 0 <= shard < expected or m["index_within_shard"] != seen[shard]:
                raise ValueError("Summary shard/index ordering mismatch: " + str(folder))
            seen[shard] += 1
            if m["coverage"] != sum(m[k] for k in ("alt_count", "ref_count", "other_count")):
                raise ValueError("Candidate coverage mismatch: " + str(folder))
            if m["af"] != m["alt_count"] / m["coverage"] or m["selected_alignments"] != min(shape[1], m["coverage"]):
                raise ValueError("Candidate AF/row count mismatch: " + str(folder))
            if (m.get("sample_unit") == SITE_UNIT) != site_mode:
                raise ValueError("Summary record and manifest disagree on the sample unit: " + str(folder))
            if site_mode:
                # One record per site; the representative (tensor) allele comes first, then
                # the other passing alleles by non-increasing ALT count, all at the same start.
                alleles, kind = m["alleles"], KIND[m["event_type"]]
                if (m["site_id"] != f"{m['node_id']}:{m['start']}:{kind}" or not alleles
                        or m["allele_count"] != len(alleles) or alleles[0]["candidate_id"] != m["candidate_id"]
                        or any(a["alt_count"] < b["alt_count"] for a, b in zip(alleles, alleles[1:]))
                        or any(a["start"] != m["start"] or KIND[a["event_type"]] != kind for a in alleles)):
                    raise ValueError(f"Site allele list mismatch at {m['site_id']}: {folder}")
                if m["site_id"] in sites:
                    raise ValueError(f"Duplicate site {m['site_id']}: {folder}")
                sites.add(m["site_id"])
            events[m["event_type"]] += 1
    if seen != sizes or sum(sizes) != manifest["tensors"]:
        raise ValueError("Summary, shard and manifest tensor counts differ: " + str(folder))
    report = dict(passed=True, tensors=sum(sizes), shards=expected, shard_size=shard_size,
                  last_shard_tensors=sizes[-1] if sizes else 0, shape_per_tensor=list(shape),
                  dtype="int8", event_types=dict(events), sites=len(sites) if site_mode else None)
    write_json(folder / "validation_report.json", report)
    return report


def validate_task(root, config, index):
    """Validate every tensor output of a task; returns {type: report}."""
    allowed = {"SNV": {"SNP"}, "INDEL": {"INS", "DEL"}, "ALL": {"SNP", "INS", "DEL"}}
    reports = {}
    for kind, folder in task_outputs(root, config, index).items():
        if kind == "shared":
            continue
        report = validate_shards(folder, config["builder"]["shard_size"])
        if set(report["event_types"]) - allowed[kind]:
            raise ValueError(f"Unexpected event type in {folder}")
        reports[kind] = report
    return reports


def task(root, index):
    config = read_json(root / "config.json")
    logs = root / "logs"
    logs.mkdir(exist_ok=True)
    command = build_command(root, config, index)
    if Path("/usr/bin/time").exists():
        command = ["/usr/bin/time", "-v", "-o", str(logs / f"task_{index:04d}.resources.txt"), *command]
    subprocess.run(command, check=True, cwd=root / "source")
    validate_task(root, config, index)


# --- run ------------------------------------------------------------------------

class MemoryRecorder:
    """Sample the controller's process tree RSS every `interval` seconds into memory.ndjson."""

    def __init__(self, root, interval=30):
        self.path = Path(root) / "memory.ndjson"
        self.interval = interval
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self.loop, daemon=True)
        self.peak_rss_kib = 0

    def sample(self):
        pending, rows, seen = [os.getpid()], [], set()
        while pending:
            pid = pending.pop()
            if pid in seen:
                continue
            seen.add(pid)
            try:
                status = {line.split(":", 1)[0]: line.split(":", 1)[1].strip()
                          for line in Path(f"/proc/{pid}/status").read_text().splitlines() if ":" in line}
                rows.append(dict(pid=pid, name=status.get("Name"),
                                 rss_kib=int(status.get("VmRSS", "0 kB").split()[0]),
                                 hwm_kib=int(status.get("VmHWM", "0 kB").split()[0])))
                pending.extend(map(int, Path(f"/proc/{pid}/task/{pid}/children").read_text().split()))
            except (FileNotFoundError, ProcessLookupError, PermissionError):
                pass
        total = sum(r["rss_kib"] for r in rows)
        self.peak_rss_kib = max(self.peak_rss_kib, total)
        with self.path.open("a") as out:
            out.write(json.dumps(dict(unix=time.time(), sum_rss_kib=total, processes=rows)) + "\n")

    def loop(self):
        while not self.stop.wait(self.interval):
            self.sample()

    def __enter__(self):
        self.sample()
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join()
        self.sample()


def execute_queue(commands, ledger_path, processes, poll_seconds=0.5, log_dir=None, cwd=None, order=None):
    """Run {task index: command} with at most `processes` concurrent; any failure stops everything.

    Tasks start in `order` (default: ascending index). Each task's stdout/stderr goes to
    `log_dir/task_NNNN.log` when a log directory is given.
    """
    queue = list(order) if order is not None else sorted(commands)
    if sorted(queue) != sorted(commands):
        raise ValueError("Queue order must list every task exactly once")
    ledger = dict(status="running", job_id=os.environ.get("SLURM_JOB_ID"), processes=processes,
                  total_tasks=len(commands), completed_tasks=0, started_unix=time.time(),
                  tasks={str(i): dict(status="pending") for i in commands})
    active, logs = {}, {}
    if log_dir is not None:
        Path(log_dir).mkdir(parents=True, exist_ok=True)

    def save():
        ledger["elapsed_seconds"] = time.time() - ledger["started_unix"]
        write_json(ledger_path, ledger)

    def interrupted(signum, frame):
        raise RuntimeError("Controller signal " + str(signum))

    previous = {s: signal.signal(s, interrupted) for s in (signal.SIGTERM, signal.SIGINT)}
    save()
    try:
        while queue or active:
            for i, proc in list(active.items()):
                code = proc.poll()
                if code is None:
                    continue
                del active[i]
                if i in logs:
                    logs.pop(i).close()
                item = ledger["tasks"][str(i)]
                item.update(status="complete" if code == 0 else "failed", returncode=code,
                            wall_seconds=time.time() - item["started_unix"])
                if code:
                    raise RuntimeError(f"Task {i} exited with code {code}")
                ledger["completed_tasks"] += 1
                print(f"COMPLETE task {i}: {ledger['completed_tasks']}/{len(commands)}", flush=True)
            while queue and len(active) < processes:
                i = queue.pop(0)
                if log_dir is not None:
                    logs[i] = (Path(log_dir) / f"task_{i:04d}.log").open("w")
                proc = subprocess.Popen(commands[i], start_new_session=True, cwd=cwd,
                                        stdout=logs.get(i), stderr=subprocess.STDOUT if i in logs else None)
                active[i] = proc
                ledger["tasks"][str(i)].update(status="running", pid=proc.pid, started_unix=time.time())
            save()
            if active:
                time.sleep(poll_seconds)
        ledger["status"] = "complete"
    except BaseException as error:
        ledger.update(status="failed", error=str(error))
        raise
    finally:
        for proc in active.values():
            if proc.poll() is None:
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        deadline = time.monotonic() + 10
        for i, proc in active.items():
            try:
                proc.wait(timeout=max(0.01, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                proc.wait()
            ledger["tasks"][str(i)].update(status="cancelled", returncode=proc.returncode)
        for log in logs.values():
            log.close()
        save()
        for s, handler in previous.items():
            signal.signal(s, handler)
    return ledger


def completed_tasks(root, config):
    """Task indices whose outputs exist, are complete and validate."""
    done = set()
    for i in range(config["tasks"]):
        outputs = task_outputs(root, config, i)
        if all((folder / "manifest.json").exists() for folder in outputs.values()):
            try:
                validate_task(root, config, i)
                done.add(i)
            except (ValueError, KeyError, OSError, json.JSONDecodeError):
                pass
    return done


def set_aside_partial(root, config, pending):
    """Move any existing output of a pending task out of the way before rerunning it."""
    moved = []
    for i in pending:
        for kind, folder in task_outputs(root, config, i).items():
            if folder.exists():
                target = root / "incomplete" / time.strftime("%Y%m%dT%H%M%S") / kind / folder.name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(folder), str(target))
                moved.append(str(target))
    return moved


def run(root, resume=False):
    root = Path(root).resolve()
    config = read_json(root / "config.json")
    if int(os.environ.get("SLURM_CPUS_PER_TASK", config["processes"])) < config["processes"]:
        raise ValueError("Not enough allocated CPUs for the configured process count")
    verify(root, config)
    done = completed_tasks(root, config) if resume else set()
    pending = [i for i in range(config["tasks"]) if i not in done]
    if not resume and any(folder.exists() for i in pending for folder in task_outputs(root, config, i).values()):
        raise ValueError("Task outputs already exist; use --resume to continue an interrupted run")
    moved = set_aside_partial(root, config, pending) if resume else []
    status = dict(status="running", job_id=os.environ.get("SLURM_JOB_ID"), started_unix=time.time(),
                  tasks=config["tasks"], processes=config["processes"], resumed=resume,
                  previously_completed=len(done), pending=len(pending), set_aside=moved)
    write_json(root / "status.json", status)
    commands = {i: [config["python"], "-m", f"{PACKAGE}.orchestrate", "task", "--root", str(root), "--index", str(i)]
                for i in pending}
    # Long tasks first (longest-processing-time order), so no expensive task starts late
    # and runs alone at the end; without recorded costs this is plain index order.
    order = sorted(pending, key=lambda i: (-config["parts"][i].get("predicted_cost", 0), i))
    memory = MemoryRecorder(root)
    try:
        with memory:
            if commands:
                execute_queue(commands, root / "queue_status.json", config["processes"],
                              log_dir=root / "logs", cwd=root / "source", order=order)
            verify(root, config)
            catalog = catalog_outputs(root, config)
            status.update(status="complete", tensors=catalog["tensors"], tensors_by_type=catalog["tensors_by_type"])
    except BaseException as error:
        status.update(status="failed", error=str(error))
        raise
    finally:
        status["elapsed_seconds"] = time.time() - status["started_unix"]
        status["peak_sampled_rss_kib"] = memory.peak_rss_kib
        write_json(root / "status.json", status)


def catalog_outputs(root, config):
    """outputs.json: every tensor directory of every task with its tensor count."""
    kinds = list(config["variant_outputs"] or ["ALL"])
    entries = {kind: [] for kind in kinds}
    for i in range(config["tasks"]):
        outputs = task_outputs(root, config, i)
        for kind in kinds:
            manifest = read_json(outputs[kind] / "manifest.json")
            entries[kind].append(dict(task=i, path=str(outputs[kind]), tensors=manifest["tensors"],
                                      shards=manifest["shards"]))
    catalog = dict(root=str(root), tensors_dir=config["tensors"],
                   tensors_by_type={k: sum(e["tensors"] for e in v) for k, v in entries.items()}, outputs=entries)
    catalog["tensors"] = sum(catalog["tensors_by_type"].values())
    write_json(root / "outputs.json", catalog)
    return catalog


def main(argv=None):
    from indexed_gam_pipeline_v2.run import add_build_arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="action", required=True)
    p = commands.add_parser("prepare", help="freeze source, partition nodes, write config.json and run.sh")
    p.add_argument("--root", required=True, help="new run directory (config, frozen source, partitions, logs, status)")
    p.add_argument("--tensors", help="tensor output directory (default: <root>/tensors)")
    p.add_argument("--gam", required=True)
    p.add_argument("--index", help="default: GAM path + .gai")
    p.add_argument("--nodes", required=True, help="sorted target node list, e.g. discovery output")
    p.add_argument("--node-stats", help="discovery node_stats.json: predict task costs and run expensive tasks first")
    p.add_argument("--tasks", type=int, default=512)
    p.add_argument("--processes", type=int, default=32)
    p.add_argument("--chr", default="")
    add_build_arguments(p, outputs=False)
    p.set_defaults(gam_cache_mb=8192)
    r = commands.add_parser("run", help="execute all pending tasks")
    r.add_argument("--root", required=True)
    r.add_argument("--resume", action="store_true", help="skip validated tasks; set aside partial outputs and rerun them")
    t = commands.add_parser("task", help="run and validate one task (used by `run`)")
    t.add_argument("--root", required=True)
    t.add_argument("--index", type=int, required=True)
    args = parser.parse_args(argv)
    if args.action == "prepare":
        prepare(args)
    elif args.action == "task":
        task(Path(args.root).resolve(), args.index)
    else:
        run(args.root, args.resume)


if __name__ == "__main__":
    main()

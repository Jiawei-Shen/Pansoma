#!/usr/bin/env python3
"""Whole-genome runs on one Slurm node: many small tasks, one disposable builder process each.

    prepare   freeze a copy of this package, split the node list into N contiguous tasks,
              record input fingerprints and write config.json + run.sh under --root;
              with --node-stats (discovery's node_stats.json) also predict each task's cost
    run       execute the tasks (at most --processes at once, fail fast), most expensive
              predicted first when costs were recorded, else in index order; --resume skips
              tasks that already completed and validated
    task      run one task: builder subprocess, then validate its shards (used by `run`)
    finalize  merge per chromosome, then label, as frozen at prepare (what `run` does at its
              end); each step is skipped when already done

Layout under --root (bookkeeping) and --tensors (outputs; default <root>/tensors):
    <root>/source/                frozen package copy without tests/ and tools/ (hashes recorded in config.json)
    <root>/parts/nodes_NNNN.txt   node list of each task (supplement rounds: parts/supplement_NN/)
    <root>/logs/task_NNNN.log     builder stdout/stderr; task_NNNN.resources.txt from /usr/bin/time -v
    <root>/queue_status.json, status.json, memory.ndjson, outputs.json
    <tensors>/shared/task_NNNN/   per-task shared manifest, batch timings, audit streams
    <tensors>/SNV/task_NNNN/      SNV shards + summary + manifest
    <tensors>/INDEL/task_NNNN/    INDEL shards + summary + manifest
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

from .build import KIND, SITE_UNIT
from .candidates import FORMAT_VERSION, STORAGE_VERSION
from .common import read_json, sha256_file, stamp, write_json
from .graph_index import GraphIndex
from .tensor_postprocessing.chr_index import select_nodes

PACKAGE = Path(__file__).resolve().parent.name
# Every builder option is frozen into config.json and passed explicitly, so a prepared
# run never depends on the CLI defaults of the code that later executes it.
BUILDER_OPTIONS = ("rows", "width", "gam_cache_mb", "batch_nodes", "max_node_span", "max_batch_alignments",
                   "shard_size", "min_mapq", "min_af", "min_variants", "min_allele_bq", "max_indel_len",
                   "chromosomes", "max_node_reads", "early_af_filter", "decoder")
LABEL_INPUTS = ("somatic_vcf", "somatic_bed", "germline_vcf", "germline_bed", "reference_fasta")


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


def node_counts(stats_path, nodes):
    """Discovery's (not_perfect, perfect + not_perfect) counts of every node in `nodes` (sorted int64 array).

    node_stats.json covers every observed node (several GB for a genome), so it is scanned as
    bytes instead of parsed into one dict.
    """
    found, perfect, edited = array("q"), array("q"), array("q")
    with open(stats_path, "rb") as stream, mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ) as data:
        for match in NODE_STATS_RECORD.finditer(data):
            found.append(int(match.group(1)))
            perfect.append(int(match.group(2)))
            edited.append(int(match.group(3)))
    found, perfect, edited = (np.frombuffer(a, dtype=np.int64) for a in (found, perfect, edited))
    order = np.argsort(found, kind="stable")
    found = found[order]
    at = np.minimum(np.searchsorted(found, nodes), max(len(found) - 1, 0))
    if not len(found) or np.any(found[at] != nodes):
        raise ValueError("node_stats.json lacks some target nodes: " + str(stats_path))
    edited = edited[order][at]
    return edited, perfect[order][at] + edited


def node_costs(stats_path, nodes):
    """Discovery's not_perfect count of every node in `nodes`: the number of edited MAPQ-passing
    mappings on a node tracks how much decoding and candidate work the node costs."""
    return node_counts(stats_path, nodes)[0]


def postprocess_options(args):
    """The finalize step (merge, then labels) frozen into config.json; its input files are fingerprinted."""
    labels = {k: getattr(args, k) for k in LABEL_INPUTS}
    if any(labels.values()) and not all(labels.values()):
        raise ValueError("Labels need all of --" + ", --".join(k.replace("_", "-") for k in LABEL_INPUTS))
    if args.merge_shard_size < 0:
        raise ValueError("--merge-shard-size must be nonnegative (0 = no merge)")
    if (args.merge_shard_size or args.chromosomes != "all") and not args.chr_index:
        raise ValueError("--chr-index is needed to merge per chromosome (or use --merge-shard-size 0) "
                         "and for --chromosomes other than all")
    if all(labels.values()) and not (args.merge_shard_size and args.reference_path and args.truth_dir):
        raise ValueError("Labels need the merge (--merge-shard-size > 0), --reference-path and --truth-dir")
    inputs = {}
    if args.reference_path:
        inputs["reference_path"] = stamp(Path(args.reference_path) / "meta.json")
    if all(labels.values()):
        inputs.update({k: stamp(v) for k, v in labels.items()})
    return dict(merge_shard_size=args.merge_shard_size, keep_sources=args.keep_sources,
                reference_path=str(Path(args.reference_path).resolve()) if args.reference_path else None,
                labels=dict({k: str(Path(v).resolve()) for k, v in labels.items()},
                            truth_dir=str(Path(args.truth_dir).resolve())) if all(labels.values()) else None,
                inputs=inputs)


def native_decoder_status(decoder):
    """{available, reason} of the native decoder the tasks would use (native.load() of this package).

    --decoder native with an unusable module is refused; under auto a warning says that the tasks
    will decode in Python. --decoder python does not check.
    """
    if decoder == "python":
        return dict(available=None, reason="not checked (--decoder python)")
    from . import native
    module, info = native.load()
    if module is not None:
        return dict(available=True, reason=None)
    if decoder == "native":
        raise ValueError("--decoder native, but the native decoder is unavailable: " + info["reason"])
    print("WARNING: the native decoder is unavailable, so the tasks will decode in Python (much slower): "
          + info["reason"], file=sys.stderr, flush=True)
    return dict(available=False, reason=info["reason"])


def read_config(root):
    """config.json of a run root prepared by this package; roots of any other package are refused."""
    config = read_json(Path(root) / "config.json")
    if config.get("package") != PACKAGE:
        raise ValueError(f"{root} was prepared by {config.get('package') or 'another package (no config.package)'}, "
                         f"not {PACKAGE}; resume, run or finalize it with the package that prepared it")
    return config


def prepare(args):
    postprocess = postprocess_options(args)
    native_decoder = native_decoder_status(args.decoder)
    root = Path(args.root).resolve()
    root.mkdir(parents=True, exist_ok=False)
    if not 1 <= args.processes <= args.tasks:
        raise ValueError("Require tasks >= processes >= 1")
    with GraphIndex(args.graph_index) as graph:
        graph_metadata = graph.metadata
    tensors = Path(args.tensors).resolve() if args.tensors else root / "tensors"
    if tensors.exists() and any(tensors.iterdir()):
        raise ValueError(f"Tensor output directory must be new or empty: {tensors}")
    package = Path(__file__).resolve().parent
    shutil.copytree(package, root / "source" / PACKAGE,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "tests", "tools", ".fastdecode-build-*"))
    # --chromosomes: drop target nodes outside the chosen blocks before partitioning, so tasks
    # stay balanced and cost predictions cover only what is built.
    nodes, selection = select_nodes(np.fromfile(str(args.nodes), dtype=np.int64, sep="\n"),
                                    args.chromosomes, args.chr_index)
    nodes_file = Path(args.nodes)
    if args.chromosomes != "all":
        nodes_file = root / "nodes_selected.txt"
        nodes_file.write_text("".join(f"{n}\n" for n in nodes.tolist()))
        selection.update(nodes_file=str(nodes_file), sha256=sha256_file(nodes_file))
    parts = partition_nodes(nodes_file, root / "parts", args.tasks)
    schedule = dict(order="task index")
    # Deep nodes (more MAPQ>5 mappings in discovery than --downsample-reads) are built alone from a
    # fixed sample of that many records; without node stats only the builder's single-node limit applies.
    downsample = dict(reads=args.downsample_reads, nodes=0,
                      rule="discovery perfect + not_perfect > reads" if args.node_stats else "no node stats")
    if args.node_stats:
        edited, mappings = node_counts(args.node_stats, nodes)
        deep = mappings > args.downsample_reads
        if deep.any():
            table = root / "downsample_nodes.tsv"
            table.write_text("node\tmappings\n" + "".join(
                f"{n}\t{m}\n" for n, m in zip(nodes[deep].tolist(), mappings[deep].tolist())))
            downsample.update(nodes=int(deep.sum()), mappings=int(mappings[deep].sum()),
                              nodes_file=str(table), sha256=sha256_file(table))
        costs = np.cumsum(edited)
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
        package=PACKAGE,
        python=sys.executable,
        tensors=str(tensors),
        inputs=dict(gam=stamp(args.gam), index=stamp(index),
                    graph_index=dict(stamp(args.graph_index), metadata=graph_metadata),
                    nodes=dict(stamp(args.nodes), sha256=sha256_file(args.nodes), count=int(selection["nodes_in"])),
                    **({"chr_index": stamp(args.chr_index)} if args.chr_index else {}),
                    **postprocess.pop("inputs")),
        chromosome_selection=selection,
        source_sha256={str(p.relative_to(root)): sha256_file(p) for p in sorted((root / "source").rglob("*")) if p.is_file()},
        tasks=args.tasks, processes=args.processes, schedule=schedule, parts=parts,
        # max_tasks is frozen here (not derived from processes at run time), so a hand-lowered
        # processes after an out-of-memory failure does not move supplement task boundaries.
        supplement=dict(max_rounds=args.supplement_rounds, min_records=args.supplement_min_records,
                        max_tasks=4 * args.processes, rounds=[]),
        builder={k: getattr(args, k) for k in BUILDER_OPTIONS},
        downsample=downsample,
        native_decoder=native_decoder,
        variant_outputs=dict(SNV=args.snv_min_af, INDEL=args.indel_min_af),
        postprocess=postprocess)
    write_json(root / "config.json", config)
    (root / "run.sh").write_text("\n".join([
        "#!/bin/bash", "set -euo pipefail",
        "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1",
        "export PYTHONUNBUFFERED=1",
        "cd " + shlex.quote(str(root / "source")),
        "exec " + shlex.quote(sys.executable) + f" -m {PACKAGE}.orchestrate run --root " + shlex.quote(str(root)) + ' "$@"', ""]))
    (root / "run.sh").chmod(0o755)
    print(json.dumps(dict(root=str(root), tensors=str(tensors), tasks=args.tasks, processes=args.processes,
                          nodes=int(len(nodes)), chromosome_selection={k: v for k, v in selection.items() if k != "sha256"},
                          downsample={k: v for k, v in downsample.items() if k != "sha256"},
                          variant_outputs=config["variant_outputs"], postprocess=postprocess,
                          native_decoder=native_decoder), indent=2))
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
    table = config.get("downsample", {}).get("nodes_file")
    if table and sha256_file(table) != config["downsample"]["sha256"]:
        raise ValueError("Downsample table changed: " + table)


# --- task -----------------------------------------------------------------------

def task_outputs(root, config, index):
    """{output type: directory} for one task, all under config['tensors']: 'shared' (coordinator
    files) plus one directory per variant type (SNV, INDEL)."""
    base, name = Path(config["tensors"]), f"task_{index:04d}"
    return {"shared": base / "shared" / name, **{kind: base / kind / name for kind in config["variant_outputs"]}}


def build_command(root, config, index):
    b = config["builder"]
    outputs = task_outputs(root, config, index)
    command = [config["python"], "-m", f"{PACKAGE}.run", "build",
               "--gam", config["inputs"]["gam"]["path"], "--index", config["inputs"]["index"]["path"],
               "--graph-index", config["inputs"]["graph_index"]["path"],
               "--nodes", config["parts"][index]["nodes_file"],
               "--output", str(outputs["shared"])]
    if "chr_index" in config["inputs"]:
        command += ["--chr-index", config["inputs"]["chr_index"]["path"]]
    for key in BUILDER_OPTIONS:
        if key == "early_af_filter":
            command += ["--early-af-filter" if b[key] else "--no-early-af-filter"]
        else:
            command += ["--" + key.replace("_", "-"), str(b[key])]
    if "downsample" in config:  # roots prepared before downsampling keep the builder's default
        command += ["--downsample-reads", str(config["downsample"]["reads"])]
        if config["downsample"].get("nodes_file"):
            command += ["--downsample-nodes", config["downsample"]["nodes_file"]]
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
    if manifest.get("sample_unit") != SITE_UNIT:
        raise ValueError("Unexpected sample unit (site tensors expected): " + str(folder))
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
            if m["af"] != m["alt_count"] / m["coverage"] or m["selected_alignments"] != min(shape[1], m["site_coverage"]):
                raise ValueError("Candidate AF/row count mismatch: " + str(folder))
            # Rows: blocks A1..Ak, REF, OTHER in that order, contiguous, matching the site counts.
            blocks = list(m["allele_labels"]) + ["REF", "OTHER"]
            groups, n = m["row_groups"], m["selected_alignments"]
            if (sum(m["site_counts"].values()) != m["site_coverage"] or sum(m["selected_counts"].values()) != n
                    or [g["allele"] for g in groups] != sorted((g["allele"] for g in groups), key=blocks.index)
                    or [(g["start_row"], g["end_row"]) for g in groups] != [
                        (sum(g2["end_row"] - g2["start_row"] for g2 in groups[:i]),
                         sum(g2["end_row"] - g2["start_row"] for g2 in groups[:i + 1])) for i in range(len(groups))]
                    or (groups and groups[-1]["end_row"] != n)
                    or any(m["selected_counts"].get(g["allele"]) != g["end_row"] - g["start_row"] for g in groups)):
                raise ValueError("Site row blocks and counts disagree: " + str(folder))
            if m.get("sample_unit") != SITE_UNIT:
                raise ValueError("Summary record and manifest disagree on the sample unit: " + str(folder))
            # One record per site; the representative (tensor) allele comes first, then
            # the other passing alleles by non-increasing ALT count, all at the same start.
            alleles, kind = m["alleles"], KIND[m["event_type"]]
            if (m["site_id"] != f"{m['node_id']}:{m['start']}:{kind}" or not alleles
                    or m["allele_count"] != len(alleles) or alleles[0]["candidate_id"] != m["candidate_id"]
                    or any(a["alt_count"] < b["alt_count"] for a, b in zip(alleles, alleles[1:]))
                    or any(a["start"] != m["start"] or KIND[a["event_type"]] != kind for a in alleles)):
                raise ValueError(f"Site allele list mismatch at {m['site_id']}: {folder}")
            if [a["label"] for a in alleles] != list(m["allele_labels"]) or any(
                    m["allele_labels"][a["label"]] != a["candidate_id"] for a in alleles):
                raise ValueError(f"Site allele labels mismatch at {m['site_id']}: {folder}")
            if m["site_id"] in sites:
                raise ValueError(f"Duplicate site {m['site_id']}: {folder}")
            sites.add(m["site_id"])
            events[m["event_type"]] += 1
    if seen != sizes or sum(sizes) != manifest["tensors"]:
        raise ValueError("Summary, shard and manifest tensor counts differ: " + str(folder))
    report = dict(passed=True, tensors=sum(sizes), shards=expected, shard_size=shard_size,
                  last_shard_tensors=sizes[-1] if sizes else 0, shape_per_tensor=list(shape),
                  dtype="int8", event_types=dict(events), sites=len(sites))
    write_json(folder / "validation_report.json", report)
    return report


def validate_task(root, config, index):
    """Validate the SNV and INDEL outputs of a task; returns {type: report}."""
    reports = {}
    for kind, folder in task_outputs(root, config, index).items():
        if kind == "shared":
            continue
        reports[kind] = validate_shards(folder, config["builder"]["shard_size"])
    return reports


def task(root, index):
    config = read_config(root)
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
    `log_dir/task_NNNN.log` when a log directory is given. The ledger is written at the start,
    whenever a task starts or ends, and at the end.
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
            changed = False
            for i, proc in list(active.items()):
                code = proc.poll()
                if code is None:
                    continue
                changed = True
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
                changed = True
            if changed:
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
    config = read_config(root)
    if (root / "outputs.json").exists() and "merge" in read_json(root / "outputs.json"):
        raise ValueError("This run was already merged (task outputs may be gone); use `orchestrate finalize` "
                         "to finish labeling instead of run/--resume")
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
            config = run_supplement(root, config)
            verify(root, config)
            catalog = catalog_outputs(root, config)
            status.update(status="complete", tensors=catalog["tensors"], tensors_by_type=catalog["tensors_by_type"])
            write_json(root / "status.json", status)
            finalize(root, config)  # a no-op with merge_shard_size 0
            status.update({k: v for k, v in read_json(root / "status.json").items() if k not in status or
                           k in ("status", "merged", "merge_layout", "labeled")})
    except BaseException as error:
        status.update(status="failed", error=str(error))
        raise
    finally:
        status["elapsed_seconds"] = time.time() - status["started_unix"]
        status["peak_sampled_rss_kib"] = memory.peak_rss_kib
        write_json(root / "status.json", status)


def run_supplement(root, config):
    """Supplement rounds: nodes that left-normalization moved target-node indels onto and no
    task covers become extra tasks of this run (appended to config.json), built and validated
    like the others, so `finalize` merges and labels them together with the main tasks.

    Each round lists the displaced nodes of every task so far minus every node already in a
    task (and outside the chromosome selection), splits them into up to supplement.max_tasks
    tasks (4 x processes, frozen at prepare) and runs them; it stops at an empty list or after
    `max_rounds`. Rounds are recorded in the config, so --resume simply continues. Returns the
    (possibly extended) config.
    """
    settings = config["supplement"]
    while len(settings["rounds"]) < settings["max_rounds"] and not (
            settings["rounds"] and settings["rounds"][-1]["nodes"] == 0):
        number = len(settings["rounds"]) + 1
        folder = root / "parts" / f"supplement_{number:02d}"
        folder.mkdir(parents=True, exist_ok=True)
        summary = displaced_nodes(root, folder / "nodes.txt", [p["nodes_file"] for p in config["parts"]],
                                  settings["min_records"])
        first = config["tasks"]
        if summary["nodes"]:
            count = min(summary["nodes"], settings["max_tasks"])
            parts = partition_nodes(folder / "nodes.txt", folder, count)
            for part in parts:
                part["supplement_round"] = number
            config["parts"] += parts
            config["tasks"] += count
        settings["rounds"].append(dict(summary, round=number, first_task=first, tasks=config["tasks"] - first))
        config["supplement"] = settings
        write_json(root / "config.json", config)
        print(f"SUPPLEMENT round {number}: {summary['nodes']} nodes -> {config['tasks'] - first} tasks", flush=True)
        pending = list(range(first, config["tasks"]))  # new tasks; an interrupted round resumes via run --resume
        if pending:
            commands = {i: [config["python"], "-m", f"{PACKAGE}.orchestrate", "task", "--root", str(root),
                            "--index", str(i)] for i in pending}
            execute_queue(commands, root / f"queue_status_supplement_{number:02d}.json", config["processes"],
                          log_dir=root / "logs", cwd=root / "source", order=pending)
    return config


def finalize(root, config=None):
    """After every task validated: merge per chromosome, then label (as configured at prepare).

    Each step is skipped when already done (outputs.json has `merge`; every merged directory has
    labels.manifest.json), so an interrupted finalize can simply be repeated. Returns what was
    done ({} when nothing was left, or when merge_shard_size is 0).
    """
    from .tensor_postprocessing.merge_shards import merge
    from .tensor_postprocessing.truth_labels import label_run
    root = Path(root).resolve()
    config = config or read_config(root)
    post = config["postprocess"]
    report = {}
    if post["merge_shard_size"]:
        if "merge" not in read_json(root / "outputs.json"):
            report["merge"] = merge(root, config["inputs"]["chr_index"]["path"], post["merge_shard_size"],
                                    keep_sources=post["keep_sources"], workers=min(8, config["processes"]),
                                    reference_path=post["reference_path"])
        labels = post["labels"]
        kinds = list(config["variant_outputs"])
        tensors = Path(config["tensors"])
        if labels and not all((tensors / kind / "labels.manifest.json").exists() for kind in kinds):
            report["labels"] = label_run(tensors, kinds, post["reference_path"], labels["reference_fasta"],
                                         labels["somatic_vcf"], labels["somatic_bed"], labels["germline_vcf"],
                                         labels["germline_bed"], labels["truth_dir"])
            status = read_json(root / "status.json")
            status.update(labeled=True)
            write_json(root / "status.json", status)
    return report


def displaced_nodes(root, output, exclude=(), min_records=3):
    """Write the node list of a supplement run; returns a summary.

    Left-normalization can move an indel from a target node onto a node that is not a
    target of any task (tensors are only built on target nodes), so the site would be
    missing. Every builder lists such nodes in displaced_nodes.tsv; this collects them over
    all tasks, drops the run's own nodes and every node in `exclude` (e.g. earlier supplement
    lists) and nodes seen in fewer than `min_records` records, and writes the rest sorted.
    A supplement run over them (same options) builds those sites with their full read sets;
    its own displaced_nodes can feed a further round until the list is empty.
    """
    root = Path(root).resolve()
    config = read_config(root)
    counts = Counter()
    for i in range(config["tasks"]):
        table = task_outputs(root, config, i)["shared"] / "displaced_nodes.tsv"
        for line in table.read_text().splitlines():
            node, records = line.split("\t")
            counts[int(node)] += int(records)
    selection = config["chromosome_selection"]
    covered = [np.fromfile(str(selection.get("nodes_file", config["inputs"]["nodes"]["path"])), dtype=np.int64, sep="\n")]
    covered += [np.fromfile(str(path), dtype=np.int64, sep="\n") for path in exclude]
    covered = np.unique(np.concatenate(covered))
    nodes = np.array(sorted(n for n, c in counts.items() if c >= min_records), dtype=np.int64)
    at = np.minimum(np.searchsorted(covered, nodes), max(len(covered) - 1, 0))
    fresh = nodes[covered[at] != nodes] if len(covered) else nodes
    fresh, _ = select_nodes(fresh, selection["selection"], config["inputs"].get("chr_index", {}).get("path"))
    Path(output).write_text("".join(f"{n}\n" for n in fresh.tolist()))
    summary = dict(output=str(Path(output).resolve()), nodes=int(len(fresh)), displaced_nodes=len(counts),
                   already_covered=int(len(nodes) - len(fresh)), below_min_records=len(counts) - int(len(nodes)),
                   records=int(sum(counts[n] for n in fresh.tolist())))
    return summary


def catalog_outputs(root, config):
    """outputs.json: every tensor directory of every task with its tensor count. Sampled deep nodes
    of all tasks go to <root>/downsampled_nodes.tsv (task directories may be deleted by the merge)."""
    kinds = list(config["variant_outputs"])
    entries = {kind: [] for kind in kinds}
    sampled = []
    for i in range(config["tasks"]):
        outputs = task_outputs(root, config, i)
        for kind in kinds:
            manifest = read_json(outputs[kind] / "manifest.json")
            entries[kind].append(dict(task=i, path=str(outputs[kind]), tensors=manifest["tensors"],
                                      shards=manifest["shards"]))
        table = outputs["shared"] / "downsampled_nodes.tsv"
        if table.exists():
            sampled += [f"{i}\t{line}\n" for line in table.read_text().splitlines()[1:]]
    catalog = dict(root=str(root), tensors_dir=config["tensors"],
                   tensors_by_type={k: sum(e["tensors"] for e in v) for k, v in entries.items()}, outputs=entries)
    catalog["tensors"] = sum(catalog["tensors_by_type"].values())
    if sampled:  # only then, so outputs.json of other runs is unchanged
        (root / "downsampled_nodes.tsv").write_text("task\tnode\trecords\tkept\treason\n" + "".join(sampled))
        catalog["downsampled_nodes"] = len(sampled)
    write_json(root / "outputs.json", catalog)
    return catalog


def make_parser():
    from .run import add_build_arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="action", required=True)
    p = commands.add_parser("prepare", help="freeze source, partition nodes, write config.json and run.sh")
    p.add_argument("--root", required=True, help="new run directory (config, frozen source, partitions, logs, status)")
    p.add_argument("--tensors", help="tensor output directory (default: <root>/tensors)")
    p.add_argument("--gam", required=True)
    p.add_argument("--index", help="default: GAM path + .gai")
    p.add_argument("--nodes", required=True, help="sorted target node list, e.g. discovery output")
    p.add_argument("--node-stats", help="discovery node_stats.json: predict task costs and run expensive tasks first")
    p.add_argument("--supplement-rounds", type=int, default=0,
                   help="after the tasks, up to N rounds of supplement tasks for nodes that left-normalization moved "
                        "target-node indels onto (default 0 = off: targets from `discover --normalized` already "
                        "contain them; use 3 with raw-rule targets)")
    p.add_argument("--supplement-min-records", type=int, default=3,
                   help="supplement nodes need at least this many displaced records")
    p.add_argument("--tasks", type=int, default=512)
    p.add_argument("--processes", type=int, default=32)
    add_build_arguments(p, outputs=False)
    p.set_defaults(gam_cache_mb=8192)
    f = p.add_argument_group("finalize (after every task validated): merge per chromosome, then labels")
    f.add_argument("--merge-shard-size", type=int, default=32768,
                   help="tensors per merged <chrom>_shard_* file (0 = keep the task layout, no merge/labels)")
    f.add_argument("--keep-sources", action="store_true", help="keep the task_* directories after a verified merge")
    f.add_argument("--reference-path", help="tools.graph_prep ref-path-scan directory (GRCh38 coordinates)")
    f.add_argument("--somatic-vcf")
    f.add_argument("--somatic-bed")
    f.add_argument("--germline-vcf")
    f.add_argument("--germline-bed")
    f.add_argument("--reference-fasta", help="GRCh38 FASTA with .fai (labels)")
    f.add_argument("--truth-dir", help="where the truth tables (<set>.graph.tsv) are written")
    r = commands.add_parser("run", help="execute all pending tasks")
    r.add_argument("--root", required=True)
    r.add_argument("--resume", action="store_true", help="skip validated tasks; set aside partial outputs and rerun them")
    z = commands.add_parser("finalize", help="merge + label a completed run (what `run` does at its end)")
    z.add_argument("--root", required=True)
    t = commands.add_parser("task", help="run and validate one task (used by `run`)")
    t.add_argument("--root", required=True)
    t.add_argument("--index", type=int, required=True)
    return parser


def main(argv=None):
    args = make_parser().parse_args(argv)
    if args.action == "prepare":
        prepare(args)
    elif args.action == "task":
        task(Path(args.root).resolve(), args.index)
    elif args.action == "finalize":
        root = Path(args.root).resolve()
        config = read_config(root)
        verify(root, config)
        print(json.dumps(finalize(root, config), indent=2, default=str))
    else:
        run(args.root, args.resume)


if __name__ == "__main__":
    main()

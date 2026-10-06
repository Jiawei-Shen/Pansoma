"""Task queue behaviour, node partitioning, a full prepare -> run -> resume cycle, the package
guard, the frozen closure, standalone finalize, finalize labels with AF floors and static checks of the package."""
import ast
from contextlib import redirect_stderr, redirect_stdout
import csv
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from . import golden
from .fixtures import af_gam, graph_fixture, tiny_gam  # noqa: F401
from ..tensor_postprocessing.chr_index import FIELDS
from .. import native, orchestrate
from ..orchestrate import execute_queue, partition_nodes

PACKAGE_DIR = Path(orchestrate.__file__).resolve().parent
REPO = PACKAGE_DIR.parent
# the frozen run source must import without tests/ or tools/ (never copied by prepare)
RUNTIME_MODULES = ("orchestrate", "build", "run", "native", "tensor_postprocessing.merge_shards",
                   "tensor_postprocessing.truth_labels")


def write_chr_table(path):
    """chr1 = nodes 1-35 (autosome), chrX = nodes 36-99 (non_autosomal)."""
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows([dict(chrom="chr1", first_node=1, last_node=35, nodes=35, dataset="autosome", source="t"),
                          dict(chrom="chrX", first_node=36, last_node=99, nodes=64, dataset="non_autosomal", source="t")])
    return path


def clean_environment(**extra):
    """The environment of a fresh interpreter that must not see the checkout through PYTHONPATH."""
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "SLURM_CPUS_PER_TASK")}
    env.update(extra)
    return env


class QueueTest(unittest.TestCase):
    def test_refill_fresh_processes_and_logs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            script = root / "worker.py"
            script.write_text("import os,sys,time,json\nfrom pathlib import Path\np=Path(sys.argv[1])\nstart=time.time()\n"
                              "time.sleep(float(sys.argv[2]))\nprint('done', sys.argv[1])\n"
                              "p.write_text(json.dumps(dict(pid=os.getpid(),start=start,end=time.time())))\n")
            commands = {i: [sys.executable, str(script), str(root / f"{i}.json"), str(delay)]
                        for i, delay in enumerate([.8, .05, .05, .05])}
            result = execute_queue(commands, root / "queue.json", 2, .01, log_dir=root / "logs")
            rows = [json.loads((root / f"{i}.json").read_text()) for i in range(4)]
            self.assertEqual(result["completed_tasks"], 4)
            self.assertEqual(len({r["pid"] for r in rows}), 4)
            self.assertLess(rows[2]["start"], rows[0]["end"])
            events = sorted([(r["start"], 1) for r in rows] + [(r["end"], -1) for r in rows])
            active = 0
            for _, delta in events:
                active += delta
                self.assertLessEqual(active, 2)
            self.assertIn("done", (root / "logs/task_0003.log").read_text())

    def test_failure_cancels_active_and_skips_pending(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            commands = {5: [sys.executable, "-c", "import time;time.sleep(.1);raise SystemExit(7)"],
                        6: [sys.executable, "-c", "import time;time.sleep(30)"],
                        7: [sys.executable, "-c", "raise SystemExit(0)"]}
            with self.assertRaisesRegex(RuntimeError, "code 7"):
                execute_queue(commands, root / "queue.json", 2, .01)
            status = json.loads((root / "queue.json").read_text())
            self.assertEqual(status["status"], "failed")
            self.assertEqual([status["tasks"][k]["status"] for k in ("5", "6", "7")], ["failed", "cancelled", "pending"])

    def test_given_order_is_the_start_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            commands = {i: [sys.executable, "-c", "import time;time.sleep(.02)"] for i in range(4)}
            ledger = execute_queue(commands, root / "queue.json", 1, .005, order=[2, 0, 3, 1])
            starts = sorted(ledger["tasks"], key=lambda k: ledger["tasks"][k]["started_unix"])
            self.assertEqual(starts, ["2", "0", "3", "1"])
            with self.assertRaisesRegex(ValueError, "exactly once"):
                execute_queue(commands, root / "queue.json", 1, .005, order=[2, 0, 3])

    def test_ledger_is_written_only_when_a_task_starts_or_ends(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            commands = {i: [sys.executable, "-c", "import time;time.sleep(.3)"] for i in range(3)}
            writes = []
            real = orchestrate.write_json
            with patch.object(orchestrate, "write_json", lambda path, value: (writes.append(path), real(path, value))), \
                    redirect_stdout(io.StringIO()):
                ledger = execute_queue(commands, root / "queue.json", 2, .002)
            self.assertEqual(ledger["completed_tasks"], 3)
            # start + at most one write per start/end event + the final one (not one per poll)
            self.assertLessEqual(len(writes), 1 + 2 * len(commands) + 1)
            self.assertEqual(json.loads((root / "queue.json").read_text())["status"], "complete")

    def test_node_costs_scan_discovery_stats(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "node_stats.json"
            stats = {str(n): dict(perfect=n % 7, not_perfect=n * 3, max_read_length=9) for n in (5, 1000, 12, 40)}
            path.write_text(json.dumps(stats, indent=2))  # the layout write_json produces
            self.assertEqual(orchestrate.node_costs(path, np.array([5, 12, 1000])).tolist(), [15, 36, 3000])
            with self.assertRaisesRegex(ValueError, "lacks"):
                orchestrate.node_costs(path, np.array([5, 6]))

    def test_partition_covers_every_node_exactly_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.txt"
            source.write_text("".join(f"{i * 3}\n" for i in range(1, 10001)))
            parts = partition_nodes(source, root / "parts", 512)
            self.assertEqual(len(parts), 512)
            self.assertEqual("".join(Path(p["nodes_file"]).read_text() for p in parts), source.read_text())
            self.assertEqual(sum(p["nodes"] for p in parts), 10000)
            self.assertEqual(parts[0]["first_node"], 3)
            source.write_text("3\n3\n")
            with self.assertRaisesRegex(ValueError, "sorted"):
                partition_nodes(source, root / "bad", 1)


class EndToEndTest(unittest.TestCase):
    def prepare(self, root, gam, graph, nodes, **overrides):
        argv = ["prepare", "--root", str(root), "--gam", str(gam), "--nodes", str(nodes), "--graph-index", str(graph),
                "--haplotypes", "90", "--tasks", "3", "--processes", "2", "--gam-cache-mb", "1", "--batch-nodes", "2",
                "--shard-size", "2",
                "--min-variants", str(overrides.get("min_variants", 1)),
                "--merge-shard-size", str(overrides.get("merge_shard_size", 0)),
                "--snv-min-af", overrides.get("snv_min_af", ".06"), "--indel-min-af", overrides.get("indel_min_af", ".08")]
        if overrides.get("chr_index"):
            argv += ["--chr-index", str(overrides["chr_index"]), "--chromosomes", overrides.get("chromosomes", "all")]
        if overrides.get("keep_sources"):
            argv += ["--keep-sources"]
        if overrides.get("tensors"):
            argv += ["--tensors", str(overrides["tensors"])]
        if overrides.get("node_stats"):
            argv += ["--node-stats", str(overrides["node_stats"])]
        if overrides.get("decoder"):
            argv += ["--decoder", overrides["decoder"]]
        if overrides.get("max_node_reads"):
            argv += ["--max-node-reads", str(overrides["max_node_reads"])]
        argv += overrides.get("extra", [])
        with redirect_stdout(io.StringIO()), patch.dict(os.environ, dict(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")):
            orchestrate.main(argv)
        return json.loads((root / "config.json").read_text())

    def test_split_run_follows_the_cost_order_and_resumes_damaged_tasks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = tiny_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 86) for n in (10, 20, 30, 1000)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n1000\n")
            stats = root / "node_stats.json"
            stats.write_text(json.dumps({str(n): dict(perfect=1, not_perfect=c, max_read_length=6)
                                         for n, c in ((10, 1), (20, 2), (30, 9), (1000, 5), (7, 99))}, indent=2))
            run = root / "run"
            config = self.prepare(run, gam, graph, nodes, node_stats=stats)
            self.assertEqual([p["nodes"] for p in config["parts"]], [2, 1, 1])
            self.assertEqual([p["predicted_cost"] for p in config["parts"]], [3, 9, 5])
            self.assertTrue(config["schedule"]["order"].startswith("predicted cost descending"))
            self.assertEqual(config["package"], orchestrate.PACKAGE)
            self.assertEqual(config["native_decoder"]["available"], native.load()[0] is not None)
            self.assertNotIn("supplement", config)
            self.assert_frozen_closure(run)
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                orchestrate.run(run)
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status["tensors"], status["tensors_by_type"]),
                             ("complete", 4, dict(SNV=4, INDEL=0)))
            ledger = json.loads((run / "queue_status.json").read_text())["tasks"]
            self.assertEqual(sorted(ledger, key=lambda k: ledger[k]["started_unix"]), ["1", "2", "0"])
            catalog = json.loads((run / "outputs.json").read_text())
            self.assertEqual([e["tensors"] for e in catalog["outputs"]["SNV"]], [2, 1, 1])
            self.assertEqual(catalog["tensors_dir"], str(run / "tensors"))
            shards = [np.load(f) for i in range(3) for f in sorted((run / "tensors" / "SNV" / f"task_{i:04d}").glob("shard_*_data.npy"))]
            combined = np.concatenate(shards)
            self.assertEqual(combined.shape, (4, 8, 200, 101))
            self.assertTrue((run / "logs/task_0000.log").exists())
            self.assertTrue((run / "memory.ndjson").exists())
            # A second plain run refuses to overwrite; --resume redoes only the damaged task.
            with self.assertRaisesRegex(ValueError, "resume"):
                orchestrate.run(run)
            (run / "tensors/SNV/task_0001/manifest.json").unlink()  # interrupted before completion
            (run / "tensors/SNV/task_0002/shard_00000_data.npy").unlink()  # complete manifest, damaged shard
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                orchestrate.run(run, resume=True)
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["previously_completed"], status["pending"], status["tensors"]), (1, 2, 4))
            self.assertEqual(len(status["set_aside"]), 6)  # shared, SNV and INDEL of both tasks
            self.assertTrue(list((run / "incomplete").rglob("task_0001")))
            resumed = np.concatenate([np.load(f) for i in range(3)
                                      for f in sorted((run / "tensors" / "SNV" / f"task_{i:04d}").glob("shard_*_data.npy"))])
            np.testing.assert_array_equal(resumed, combined)
            # A root whose config names another package is refused before anything runs.
            original = (run / "config.json").read_text()
            (run / "config.json").write_text(json.dumps(dict(json.loads(original), package="other_package")))
            with self.assertRaisesRegex(ValueError, "prepared by other_package"):
                orchestrate.run(run, resume=True)
            (run / "config.json").write_text(original)
            # Any change to a frozen input is refused.
            with nodes.open("a") as stream:
                stream.write("\n")
            with self.assertRaisesRegex(ValueError, "Input changed"):
                orchestrate.run(run, resume=True)

    def test_split_run_applies_independent_thresholds_and_validates_types(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = af_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n40\n50\n")
            run = root / "run"
            tensors = root / "elsewhere"
            self.prepare(run, gam, graph, nodes, min_variants=3, tensors=tensors)
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                orchestrate.run(run)
            status = json.loads((run / "status.json").read_text())
            self.assertEqual(status["tensors_by_type"], dict(SNV=1, INDEL=2))
            records = [json.loads(l) for i in range(3)
                       for l in (tensors / "INDEL" / f"task_{i:04d}" / "variant_summary.ndjson").read_text().splitlines()]
            self.assertEqual([(m["node_id"], m["event_type"], m["af"]) for m in records], [(30, "DEL", .08), (50, "INS", .08)])
            for i in range(3):
                for kind in ("SNV", "INDEL"):
                    report = json.loads((tensors / kind / f"task_{i:04d}" / "validation_report.json").read_text())
                    self.assertTrue(report["passed"])
                    self.assertEqual(report["dtype"], "int8")
            self.assertFalse(list((tensors / "shared/task_0000").glob("shard_*")))
            self.assertFalse((run / "SNV").exists())


    def test_autosome_selection_then_automatic_merge(self):
        """--chromosomes autosome drops chrX targets before partitioning; `run` ends with a verified merge."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = af_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n40\n50\n")
            table = write_chr_table(root / "chr.tsv")
            plain = root / "plain"
            self.prepare(plain, gam, graph, nodes, min_variants=3, chr_index=table, chromosomes="autosome")
            run = root / "run"
            config = self.prepare(run, gam, graph, nodes, min_variants=3, chr_index=table,
                                  chromosomes="autosome", merge_shard_size=4)
            self.assertEqual(sum(p["nodes"] for p in config["parts"]), 3)
            self.assertEqual((config["chromosome_selection"]["nodes_kept"], config["chromosome_selection"]["removed"]),
                             (3, {"chrX": 2}))
            self.assertIn("--chromosomes autosome", " ".join(orchestrate.build_command(run, config, 0)))
            for folder in (plain, run):
                with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                    orchestrate.run(folder)
            manifest = json.loads((plain / "tensors/SNV/task_0000/manifest.json").read_text())
            self.assertEqual(manifest["chromosome_selection"]["selection"], "autosome")
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status["merged"]), ("finalized", True))
            self.assertFalse(list((run / "tensors").glob("*/task_*")))  # default: sources deleted after verification
            for kind in ("SNV", "INDEL"):
                before = [np.load(f) for i in range(3)
                          for f in sorted((plain / "tensors" / kind / f"task_{i:04d}").glob("shard_*_data.npy"))]
                merged = json.loads((run / "tensors" / kind / "manifest.json").read_text())
                after = [np.load(run / "tensors" / kind / s["file"]) for s in merged["chromosomes"]["chr1"]["shards"]]
                np.testing.assert_array_equal(np.concatenate(after), np.concatenate(before))
                self.assertEqual(set(merged["chromosomes"]), {"chr1"})
            with self.assertRaisesRegex(ValueError, "already merged"):
                orchestrate.run(run, resume=True)
            with redirect_stdout(io.StringIO()):
                self.assertEqual(orchestrate.finalize(run), {})  # nothing left to do
            self.assertFalse((run / "finalize_report.json").exists())

    def test_nodes_over_the_read_cap_are_listed_validated_and_merged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = af_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n40\n50\n")
            stats = root / "node_stats.json"
            stats.write_text(json.dumps({str(n): dict(perfect=47, not_perfect=3, max_read_length=6)
                                         for n in (10, 20, 30, 40, 50)}, indent=2))
            run = root / "run"  # 50 records per node, capped at 45
            config = self.prepare(run, gam, graph, nodes, min_variants=1, chr_index=write_chr_table(root / "chr.tsv"),
                                  chromosomes="autosome", merge_shard_size=4, node_stats=stats, max_node_reads=45,
                                  snv_min_af=".01", indel_min_af=".01")
            self.assertNotIn("downsample", config)
            self.assertEqual(config["builder"]["max_node_reads"], 45)
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                orchestrate.run(run)
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status["merged"]), ("finalized", True))
            task = {n: next(i for i, p in enumerate(config["parts"]) if p["first_node"] <= n <= p["last_node"])
                    for n in (10, 20, 30)}  # 40 and 50 are chrX, left out by --chromosomes autosome
            self.assertEqual((run / "downsampled_nodes.tsv").read_text(), "task\tnode\trecords\tkept\treason\n" +
                             "".join(f"{task[n]}\t{n}\t50\t45\tmax_node_reads\n" for n in (10, 20, 30)))
            self.assertEqual(json.loads((run / "outputs.json").read_text())["downsampled_nodes"], 3)
            merged = [json.loads(l) for f in (run / "tensors").glob("*/chr1*.ndjson") for l in f.read_text().splitlines()]
            self.assertTrue(merged)
            self.assertTrue(all(r["downsampled_from"] == 50 for r in merged))

    def test_standalone_finalize_completes_an_interrupted_finalize(self):
        """`run` fails inside finalize; `python -m <package>.orchestrate finalize` (a fresh process
        from the checkout) merges, and a second call has nothing left to do."""
        from ..tensor_postprocessing import merge_shards
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = af_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n40\n50\n")
            table = write_chr_table(root / "chr.tsv")
            run = root / "run"
            self.prepare(run, gam, graph, nodes, min_variants=3, chr_index=table, chromosomes="autosome",
                         merge_shard_size=4)
            with patch.object(merge_shards, "merge", side_effect=RuntimeError("interrupted merge")), \
                    patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(RuntimeError, "interrupted merge"):
                    orchestrate.run(run)
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status.get("merged")), ("failed", None))
            self.assertNotIn("merge", json.loads((run / "outputs.json").read_text()))
            command = [sys.executable, "-m", f"{orchestrate.PACKAGE}.orchestrate", "finalize", "--root", str(run)]
            first = subprocess.run(command, cwd=REPO, env=clean_environment(), capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual(list(json.loads(first.stdout)), ["merge"])
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status["merged"]), ("finalized", True))
            merged = json.loads((run / "tensors/SNV/manifest.json").read_text())
            self.assertEqual(set(merged["chromosomes"]), {"chr1"})
            second = subprocess.run(command, cwd=REPO, env=clean_environment(), capture_output=True, text=True)
            self.assertEqual((second.returncode, json.loads(second.stdout)), (0, {}), second.stderr)
            self.assertFalse((run / "finalize_report.json").exists())

    def test_run_without_finalize_then_a_separate_finalize(self):
        """`run --no-finalize` stops at `complete` with the task outputs in place; `finalize` (another process,
        as a smaller Slurm job) then merges them like `run` would have."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = af_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n40\n50\n")
            run = root / "run"
            self.prepare(run, gam, graph, nodes, min_variants=3, chr_index=write_chr_table(root / "chr.tsv"),
                         chromosomes="autosome", merge_shard_size=4)
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                orchestrate.main(["run", "--root", str(run), "--no-finalize"])
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status.get("merged")), ("complete", None))
            self.assertNotIn("merge", json.loads((run / "outputs.json").read_text()))
            self.assertEqual(len(list((run / "tensors").glob("SNV/task_*"))), 3)
            command = [sys.executable, "-m", f"{orchestrate.PACKAGE}.orchestrate", "finalize", "--root", str(run)]
            result = subprocess.run(command, cwd=REPO, env=clean_environment(), capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(list(json.loads(result.stdout)), ["merge"])
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status["merged"]), ("finalized", True))
            self.assertFalse(list((run / "tensors").glob("*/task_*")))
            self.assertEqual(set(json.loads((run / "tensors/SNV/manifest.json").read_text())["chromosomes"]), {"chr1"})

    def test_prepare_refuses_merge_or_selection_without_chr_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = tiny_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 86) for n in (10, 20, 30, 1000)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n1000\n")
            with self.assertRaisesRegex(ValueError, "chr-index"):
                self.prepare(root / "run", gam, graph, nodes, merge_shard_size=32768)
            self.assertFalse((root / "run").exists())

    def test_prepare_refuses_label_af_floors_without_labels_or_outside_0_1(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = tiny_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 86) for n in (10, 20, 30, 1000)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n1000\n")
            for extra, message in ((["--label-snv-min-af", "0.07"], "--label-snv-min-af needs the label options"),
                                   (["--label-indel-min-af", "0.1"], "--label-indel-min-af needs the label options"),
                                   (["--label-snv-min-af", "1.5"], r"--label-snv-min-af must be in \[0, 1\]"),
                                   (["--label-indel-min-af=-0.1"], r"--label-indel-min-af must be in \[0, 1\]")):
                with self.subTest(extra=extra):
                    with self.assertRaisesRegex(ValueError, message):
                        self.prepare(root / "run", gam, graph, nodes, extra=extra)
                    self.assertFalse((root / "run").exists())  # refused before anything is created

    def test_finalize_labels_with_the_af_floors_frozen_at_prepare(self):
        """--label-snv-min-af / --label-indel-min-af go into config.postprocess.labels after truth_dir (null when
        unset); finalize labels with them, exactly as `tensor_postprocessing label` with the same floors, and with
        both null only the floored tensors change (O1's world: SNV tensors at AF 0.16-0.8, the germline DEL at 0.2,
        the somatic INS at 0.6)."""
        from ..tensor_postprocessing.__main__ import main as postprocessing
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "inputs").mkdir()
            run = root / "run"
            argv = golden.mini_world(root / "inputs", run)
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                orchestrate.main(["prepare", "--root", str(root / "plain"), *argv[2:]])
                orchestrate.main(["prepare", *argv, "--label-snv-min-af", "0.3", "--label-indel-min-af", "0.3"])
            unset, frozen = (json.loads((r / "config.json").read_text())["postprocess"]["labels"]
                             for r in (root / "plain", run))
            self.assertEqual(list(frozen), [*orchestrate.LABEL_INPUTS, "truth_dir", "snv_min_af", "indel_min_af"])
            self.assertEqual(unset, dict(frozen, snv_min_af=None, indel_min_af=None))
            self.assertEqual((frozen["snv_min_af"], frozen["indel_min_af"]), (0.3, 0.3))
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                orchestrate.run(run)
            self.assertTrue(json.loads((run / "status.json").read_text())["labeled"])
            config = json.loads((run / "config.json").read_text())
            tensors = Path(config["tensors"])

            def labelled():
                """Label .npy/.ndjson bytes; per kind labels.manifest.json without `created`, [(AF, labels.ndjson line)]."""
                files = {str(p.relative_to(tensors)): p.read_bytes() for p in sorted(tensors.glob("*/chr*_labels*"))}
                manifests, lines = {}, {}
                for kind in ("SNV", "INDEL"):
                    manifests[kind] = dict(json.loads((tensors / kind / "labels.manifest.json").read_text()), created=None)
                    lines[kind] = []
                    for f in sorted((tensors / kind).glob("chr*_labels.ndjson")):
                        summary = (tensors / kind / f.name.replace("_labels", "_variant_summary")).read_text().splitlines()
                        lines[kind] += [(json.loads(a)["af"], json.loads(b)) for a, b in zip(summary, f.read_text().splitlines())]
                return files, manifests, lines

            files, manifests, lines = labelled()
            for kind, below in (("SNV", "below_snv_min_af"), ("INDEL", "below_indel_min_af")):
                self.assertEqual((manifests[kind]["snv_min_af"], manifests[kind]["indel_min_af"]), (0.3, 0.3))
                self.assertEqual([line["reason"] == below for _, line in lines[kind]], [af < 0.3 for af, _ in lines[kind]])
                floored = [line for af, line in lines[kind] if af < 0.3]
                self.assertTrue(any(line["somatic"] or line["germline"] for line in floored), kind)  # truth tensors too
            post, labels = config["postprocess"], config["postprocess"]["labels"]
            with redirect_stdout(io.StringIO()):
                postprocessing(["label", "--tensors", str(tensors), "--reference-path", post["reference_path"],
                                "--fasta", labels["reference_fasta"], "--somatic-vcf", labels["somatic_vcf"],
                                "--somatic-bed", labels["somatic_bed"], "--germline-vcf", labels["germline_vcf"],
                                "--germline-bed", labels["germline_bed"], "--truth-dir", labels["truth_dir"],
                                "--snv-min-af", "0.3", "--indel-min-af", "0.3"])
            self.assertEqual(labelled(), (files, manifests, lines))
            # Both floors null (a prepare without them): finalize labels again; only the floored tensors differ.
            labels.update(snv_min_af=None, indel_min_af=None)
            (run / "config.json").write_text(json.dumps(config))
            for kind in ("SNV", "INDEL"):
                (tensors / kind / "labels.manifest.json").unlink()
            with redirect_stdout(io.StringIO()):
                self.assertEqual(list(orchestrate.finalize(run)), ["labels"])
            _, unfloored, again = labelled()
            for kind in ("SNV", "INDEL"):
                self.assertEqual((unfloored[kind]["snv_min_af"], unfloored[kind]["indel_min_af"]), (None, None))
                for (af, line), (_, unset_line) in zip(lines[kind], again[kind]):
                    self.assertEqual(unset_line == line, af >= 0.3, line["candidate_id"])
                truths = [line for af, line in again[kind] if af < 0.3 and (line["somatic"] or line["germline"])]
                self.assertEqual([line["label"] for line in truths], [1 if line["somatic"] else 2 for line in truths])
            # A finalize that failed in the recall scan (labels.manifest.json written, truth_recall.json not yet): a
            # repeated finalize labels again and writes the recall reports, then has nothing left to do.
            (tensors / "truth_recall.json").unlink()
            with redirect_stdout(io.StringIO()):
                self.assertEqual(list(orchestrate.finalize(run)), ["labels"])
                self.assertTrue((tensors / "truth_recall.json").exists())
                self.assertEqual(orchestrate.finalize(run), {})

    def test_prepare_refuses_an_unusable_native_decoder_and_warns_under_auto(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = tiny_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 86) for n in (10, 20, 30, 1000)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n1000\n")
            with patch.object(native, "_LOADED", (None, dict(reason="not built here"))):
                with self.assertRaisesRegex(ValueError, "--decoder native.*not built here"):
                    self.prepare(root / "native", gam, graph, nodes, decoder="native")
                self.assertFalse((root / "native").exists())  # refused before anything is created
                warning = io.StringIO()
                with redirect_stderr(warning):
                    config = self.prepare(root / "auto", gam, graph, nodes)
                self.assertIn("decode in Python", warning.getvalue())
                self.assertIn("not built here", warning.getvalue())
                self.assertEqual(config["native_decoder"], dict(available=False, reason="not built here"))
                self.assertEqual(config["builder"]["decoder"], "auto")
            config = self.prepare(root / "python", gam, graph, nodes, decoder="python")
            self.assertIsNone(config["native_decoder"]["available"])

    def test_roots_prepared_by_another_package_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for config, named in ((dict(tasks=1), r"another package \(no config.package\)"),
                                  (dict(package="other_package", tasks=1), "other_package")):
                (root / "config.json").write_text(json.dumps(config))
                calls = dict(run=lambda: orchestrate.run(root), resume=lambda: orchestrate.run(root, resume=True),
                             task=lambda: orchestrate.task(root, 0), finalize=lambda: orchestrate.finalize(root),
                             cli=lambda: orchestrate.main(["finalize", "--root", str(root)]))
                for name, call in calls.items():
                    with self.subTest(config=named, call=name), self.assertRaisesRegex(ValueError, "prepared by " + named):
                        call()

    def assert_frozen_closure(self, run):
        """<root>/source holds the runtime modules only, and they import from there alone."""
        source = run / "source"
        frozen = source / orchestrate.PACKAGE
        self.assertTrue((frozen / "build.py").exists())
        self.assertTrue((frozen / "fastdecode.cpp").exists())
        for name in ("tests", "tools"):
            self.assertFalse((frozen / name).exists(), name)
        self.assertFalse(list(frozen.rglob(".fastdecode-build-*")) + list(frozen.rglob("__pycache__")))
        code = ("import importlib, json, sys\n"
                f"files = [importlib.import_module(n).__file__ for n in {[f'{orchestrate.PACKAGE}.{m}' for m in RUNTIME_MODULES]!r}]\n"
                f"from {orchestrate.PACKAGE} import native\n"
                "print(json.dumps(dict(files=files, native=native.load()[0] is not None,\n"
                "                      modules=sorted(m for m in sys.modules if m.startswith('indexed_gam_pipeline')))))\n")
        result = subprocess.run([sys.executable, "-c", code], cwd=source, env=clean_environment(PYTHONDONTWRITEBYTECODE="1"),
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        found = json.loads(result.stdout)
        for path in found["files"]:
            self.assertTrue(Path(path).resolve().is_relative_to(frozen.resolve()), path)
        self.assertTrue(all(m == orchestrate.PACKAGE or m.startswith(orchestrate.PACKAGE + ".") for m in found["modules"]))
        self.assertFalse([m for m in found["modules"] if m.split(".")[1:2] in (["tools"], ["tests"])])
        self.assertEqual(found["native"], bool(list(PACKAGE_DIR.glob("_fastdecode*.so"))))


class StaticTest(unittest.TestCase):
    """Package-wide source checks (the runtime/tools/tests boundary, no module-search-path edits,
    no other pipeline package names)."""

    def sources(self, *suffixes):
        return sorted(p for p in PACKAGE_DIR.rglob("*") if p.suffix in suffixes and "__pycache__" not in p.parts)

    def test_runtime_modules_never_import_tools_or_tests(self):
        for path in self.sources(".py"):
            relative = path.relative_to(PACKAGE_DIR)
            if relative.parts[0] in ("tools", "tests"):
                continue
            package = [orchestrate.PACKAGE, *relative.parts[:-1]]
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node, ast.ImportFrom):
                    base = package[:len(package) - node.level + 1] if node.level else []
                    module = base + (node.module.split(".") if node.module else [])
                    targets = [module + [a.name] for a in node.names] if not node.module else [module]
                elif isinstance(node, ast.Import):
                    targets = [a.name.split(".") for a in node.names]
                else:
                    continue
                for target in targets:
                    if target[:1] == [orchestrate.PACKAGE]:
                        self.assertNotIn(target[1:2], (["tools"], ["tests"]), f"{relative}: {'.'.join(target)}")

    def test_no_path_hacks_and_no_other_package_names(self):
        hack, package_name = "sys." + "path", re.compile(r"indexed_gam_pipeline\w+")
        for path in self.sources(".py", ".cpp", ".sh"):
            relative, text = path.relative_to(PACKAGE_DIR), path.read_text()
            self.assertNotIn(hack, text, str(relative))
            self.assertEqual(set(package_name.findall(text)) - {orchestrate.PACKAGE}, set(), str(relative))


if __name__ == "__main__":
    unittest.main()

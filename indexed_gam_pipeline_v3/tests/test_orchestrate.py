"""Task queue behaviour, node partitioning, a full prepare -> run -> resume cycle, the package
guard, the frozen closure, standalone finalize and static checks of the package."""
import ast
from contextlib import redirect_stderr, redirect_stdout
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

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
            # start + at most one write per start/end event + the final one (v2 wrote on every poll)
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
                "--tasks", "3", "--processes", "2", "--gam-cache-mb", "1", "--batch-nodes", "2", "--shard-size", "2",
                "--min-variants", str(overrides.get("min_variants", 1)),
                "--merge-shard-size", str(overrides.get("merge_shard_size", 0)),
                "--snv-min-af", ".06", "--indel-min-af", ".08"]
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
            self.assertEqual(config["supplement"]["max_tasks"], 8)  # 4 x processes, frozen
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

    def test_run_appends_supplement_tasks_for_indels_normalized_off_the_targets(self):
        from .fixtures import mirror, spec_alignment, write_gam
        sequences = {1: "C", 2: "AT", 3: "AT", 4: "G"}
        right = [(1, 0, False, [(1, 1, "")]), (2, 0, False, [(2, 2, "")]), (3, 0, False, [(2, 2, ""), (0, 2, "AT")]),
                 (4, 0, False, [(1, 1, "")])]
        left = [(1, 0, False, [(1, 1, "")]), (2, 0, False, [(0, 2, "AT"), (2, 2, "")]), (3, 0, False, [(2, 2, "")]),
                (4, 0, False, [(1, 1, "")])]
        ref = [(n, 0, False, [(len(s), len(s), "")]) for n, s in sequences.items()]
        rows = ([spec_alignment(right, sequences, name=f"f{i}") for i in range(3)]
                + [spec_alignment(mirror(left, sequences), sequences, name=f"r{i}") for i in range(3)]
                + [spec_alignment(ref, sequences, name=f"ref{i}") for i in range(4)])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam = write_gam(root / "g.gam", rows)
            graph = graph_fixture(root / "graph.sqlite", [(n, s, 5) for n, s in sequences.items()])
            nodes = root / "nodes.txt"
            nodes.write_text("3\n4\n")  # vg's forward placement made node 3 a target; node 2 is none
            run = root / "run"
            argv = ["prepare", "--root", str(run), "--gam", str(gam), "--nodes", str(nodes), "--graph-index", str(graph),
                    "--tasks", "2", "--processes", "1", "--gam-cache-mb", "1", "--batch-nodes", "2", "--shard-size", "2",
                    "--min-variants", "1", "--width", "11", "--max-indel-len", "5", "--merge-shard-size", "0",
                    "--snv-min-af", ".05", "--indel-min-af", ".05"]
            with redirect_stdout(io.StringIO()), patch.dict(os.environ, dict(OMP_NUM_THREADS="1")):
                orchestrate.main(argv + ["--supplement-rounds", "3"])  # raw-rule targets: supplement on
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="1")), redirect_stdout(io.StringIO()):
                orchestrate.run(run)
            config = json.loads((run / "config.json").read_text())
            rounds = config["supplement"]["rounds"]
            self.assertEqual([(r["round"], r["nodes"], r["tasks"]) for r in rounds], [(1, 1, 1), (2, 0, 0)])
            self.assertEqual((config["supplement"]["max_tasks"], config["supplement"]["min_records"]), (4, 3))
            self.assertEqual((config["tasks"], Path(config["parts"][2]["nodes_file"]).read_text()), (3, "2\n"))
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status["tensors"]), ("complete", 1))
            (site,) = [json.loads(l) for l in (run / "tensors/INDEL/task_0002/variant_summary.ndjson").read_text().splitlines()]
            self.assertEqual((site["site_id"], site["alt_count"], site["site_counts"]), ("2:0:INDEL", 6, {"A1": 6, "REF": 4}))
            catalog = json.loads((run / "outputs.json").read_text())
            self.assertEqual([e["tensors"] for e in catalog["outputs"]["INDEL"]], [0, 0, 1])
            # the default (--supplement-rounds 0, for normalized targets) keeps the plain run
            argv[2] = str(root / "plain")
            with redirect_stdout(io.StringIO()), patch.dict(os.environ, dict(OMP_NUM_THREADS="1")):
                orchestrate.main(argv)
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="1")), redirect_stdout(io.StringIO()):
                orchestrate.run(root / "plain")
            self.assertEqual(json.loads((root / "plain/config.json").read_text())["tasks"], 2)

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
                             displaced=lambda: orchestrate.displaced_nodes(root, root / "nodes.txt"),
                             cli=lambda: orchestrate.main(["finalize", "--root", str(root)]))
                for name, call in calls.items():
                    with self.subTest(config=named, call=name), self.assertRaisesRegex(ValueError, "prepared by " + named):
                        call()
            self.assertFalse((root / "nodes.txt").exists())

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
    """Package-wide source checks (the runtime/tools/tests boundary and leftovers of v2)."""

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
        hack, other = "sys." + "path", "indexed_gam_pipeline_v" + "2"
        allowed = {Path("tests/golden.py")}  # records the goldens by running the other package (docs, argument)
        for path in self.sources(".py", ".cpp"):
            relative, text = path.relative_to(PACKAGE_DIR), path.read_text()
            self.assertNotIn(hack, text, str(relative))
            if relative not in allowed:
                self.assertNotIn(other, text, str(relative))


if __name__ == "__main__":
    unittest.main()

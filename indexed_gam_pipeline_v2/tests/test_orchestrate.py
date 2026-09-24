"""Task queue behaviour, node partitioning and a full prepare -> run -> resume cycle."""
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from fixtures import af_gam, graph_fixture, tiny_gam  # noqa: F401  (sets sys.path)
import csv
from indexed_gam_pipeline_v2.tensor_postprocessing.chr_index import FIELDS
from indexed_gam_pipeline_v2 import orchestrate
from indexed_gam_pipeline_v2.orchestrate import execute_queue, partition_nodes


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
                "--merge-shard-size", str(overrides.get("merge_shard_size", 0))]
        if overrides.get("chr_index"):
            argv += ["--chr-index", str(overrides["chr_index"]), "--chromosomes", overrides.get("chromosomes", "all")]
        if overrides.get("keep_sources"):
            argv += ["--keep-sources"]
        if overrides.get("split"):
            argv += ["--snv-min-af", ".06", "--indel-min-af", ".08"]
        if overrides.get("tensors"):
            argv += ["--tensors", str(overrides["tensors"])]
        if overrides.get("node_stats"):
            argv += ["--node-stats", str(overrides["node_stats"])]
        with redirect_stdout(io.StringIO()), patch.dict(os.environ, dict(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")):
            orchestrate.main(argv)
        return json.loads((root / "config.json").read_text())

    def test_single_output_run_matches_one_unsplit_build_and_resumes(self):
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
            self.assertTrue((run / "source" / "indexed_gam_pipeline_v2" / "build.py").exists())
            self.assertFalse((run / "source" / "indexed_gam_pipeline_v2" / "tests").exists())
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                orchestrate.run(run)
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status["tensors"]), ("complete", 4))
            ledger = json.loads((run / "queue_status.json").read_text())["tasks"]
            self.assertEqual(sorted(ledger, key=lambda k: ledger[k]["started_unix"]), ["1", "2", "0"])
            catalog = json.loads((run / "outputs.json").read_text())
            self.assertEqual([e["tensors"] for e in catalog["outputs"]["ALL"]], [2, 1, 1])
            self.assertEqual(catalog["tensors_dir"], str(run / "tensors"))
            shards = [np.load(f) for i in range(3) for f in sorted((run / "tensors" / "ALL" / f"task_{i:04d}").glob("shard_*_data.npy"))]
            combined = np.concatenate(shards)
            self.assertEqual(combined.shape, (4, 8, 200, 101))
            self.assertTrue((run / "logs/task_0000.log").exists())
            self.assertTrue((run / "memory.ndjson").exists())
            # A second plain run refuses to overwrite; --resume redoes only the damaged task.
            with self.assertRaisesRegex(ValueError, "resume"):
                orchestrate.run(run)
            (run / "tensors/ALL/task_0001/manifest.json").unlink()  # interrupted before completion
            (run / "tensors/ALL/task_0002/shard_00000_data.npy").unlink()  # complete manifest, damaged shard
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), redirect_stdout(io.StringIO()):
                orchestrate.run(run, resume=True)
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["previously_completed"], status["pending"], status["tensors"]), (1, 2, 4))
            self.assertEqual(len(status["set_aside"]), 2)
            self.assertTrue(list((run / "incomplete").rglob("task_0001")))
            resumed = np.concatenate([np.load(f) for i in range(3)
                                      for f in sorted((run / "tensors" / "ALL" / f"task_{i:04d}").glob("shard_*_data.npy"))])
            np.testing.assert_array_equal(resumed, combined)
            # Any change to a frozen input is refused.
            with nodes.open("a") as stream:
                stream.write("\n")
            with self.assertRaisesRegex(ValueError, "Input changed"):
                orchestrate.run(run, resume=True)

    def test_run_appends_supplement_tasks_for_indels_normalized_off_the_targets(self):
        from fixtures import spec_alignment, write_gam
        from test_left_align import mirror
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
                    "--min-variants", "1", "--width", "11", "--max-indel-len", "5", "--merge-shard-size", "0"]
            with redirect_stdout(io.StringIO()), patch.dict(os.environ, dict(OMP_NUM_THREADS="1")):
                orchestrate.main(argv)
            with patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="1")), redirect_stdout(io.StringIO()):
                orchestrate.run(run)
            config = json.loads((run / "config.json").read_text())
            rounds = config["supplement"]["rounds"]
            self.assertEqual([(r["round"], r["nodes"], r["tasks"]) for r in rounds], [(1, 1, 1), (2, 0, 0)])
            self.assertEqual((config["tasks"], Path(config["parts"][2]["nodes_file"]).read_text()), (3, "2\n"))
            status = json.loads((run / "status.json").read_text())
            self.assertEqual((status["status"], status["tensors"]), ("complete", 1))
            (site,) = [json.loads(l) for l in (run / "tensors/ALL/task_0002/variant_summary.ndjson").read_text().splitlines()]
            self.assertEqual((site["site_id"], site["alt_count"], site["site_counts"]), ("2:0:INDEL", 6, {"A1": 6, "REF": 4}))
            catalog = json.loads((run / "outputs.json").read_text())
            self.assertEqual([e["tensors"] for e in catalog["outputs"]["ALL"]], [0, 0, 1])
            # --supplement-rounds 0 keeps the plain run
            argv[2] = str(root / "plain")
            with redirect_stdout(io.StringIO()), patch.dict(os.environ, dict(OMP_NUM_THREADS="1")):
                orchestrate.main(argv + ["--supplement-rounds", "0"])
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
            self.prepare(run, gam, graph, nodes, split=True, min_variants=3, tensors=tensors)
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
            table = root / "chr.tsv"
            with table.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows([dict(chrom="chr1", first_node=1, last_node=35, nodes=35, dataset="autosome", source="t"),
                                  dict(chrom="chrX", first_node=36, last_node=99, nodes=64, dataset="non_autosomal", source="t")])
            plain = root / "plain"
            self.prepare(plain, gam, graph, nodes, split=True, min_variants=3, chr_index=table, chromosomes="autosome")
            run = root / "run"
            config = self.prepare(run, gam, graph, nodes, split=True, min_variants=3, chr_index=table,
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


if __name__ == "__main__":
    unittest.main()

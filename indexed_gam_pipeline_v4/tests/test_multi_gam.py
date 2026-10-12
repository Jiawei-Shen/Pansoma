"""Several GAMs read as one (e.g. one sorted GAM per library): the records of a GAM split into two GAMs must give
what the single GAM gives, in the reader (fetch, fetch_capped with the read cap over both), discovery, build and a
whole prepare -> run; and one --gam keeps recording a plain path."""
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from .fixtures import graph_fixture, split_args, write_gam
from .golden import random_chain_world
from .test_discovery import needs_module, synthetic_world, write_grouped_gam
from .test_record_index import gamsort_key
from .. import discovery, gam_record_index, orchestrate, vg_pb2
from ..build import build
from ..gam_reader import IndexedGam, MultiGam, gam_pairs, open_gams
from ..gam_record_index import ENVIRONMENT
from ..run import make_parser
from ..tools.compare_runs import compare

GAM_ARGUMENTS = ("arguments.gam", "arguments.index")


def split_world(directory, records, parts=2, order=False):
    """records -> one GAM of all of them and `parts` GAMs of every parts-th record (each in gamsort order when
    `order`); the single GAM holds the parts one after the other unless `order`, so file order matches."""
    directory = Path(directory)
    raws = [a.SerializeToString(deterministic=True) for a in records]
    if order:
        raws.sort(key=gamsort_key)
    pieces = [raws[k::parts] for k in range(parts)]
    whole = sorted(raws, key=gamsort_key) if order else [r for piece in pieces for r in piece]
    load = lambda rs: [vg_pb2.Alignment.FromString(r) for r in rs]  # noqa: E731
    single = write_gam(directory / "single.gam", load(whole))
    split = [write_gam(directory / f"part{k}.gam", load(piece)) for k, piece in enumerate(pieces)]
    return single, split


def outcome(fetched):
    return [(a.SerializeToString(deterministic=True), sorted(out)) for a, out in fetched]


class ReaderTest(unittest.TestCase):
    def test_fetch_and_capped_fetch_equal_the_single_gam(self):
        with tempfile.TemporaryDirectory() as tmp:
            sequences, records = random_chain_world(20261012, nodes=60, records=600)
            single, split = split_world(tmp, records, parts=3)
            with patch("sys.stderr", io.StringIO()):
                for path in (single, *split):
                    gam_record_index.build(path, processes=1, force=True)
            nodes = sorted(sequences)
            for mode in ("off", "require", "memory"):  # GAI walk, record path from <gam>.gri, record path in memory
                with patch.dict(os.environ, {ENVIRONMENT: mode}), patch("sys.stderr", io.StringIO()):
                    one = IndexedGam(single, cache_bytes=4096)
                    many = open_gams([str(p) for p in split], cache_bytes=4096)
                    self.assertIsInstance(many, MultiGam)
                    rng = random.Random(3)
                    for trial in range(20):
                        wanted = set(rng.sample(nodes, rng.randint(1, 12)))
                        for cap in (0, 1, 4, 9):
                            with self.subTest(mode=mode, trial=trial, cap=cap):
                                a, b = {}, {}
                                self.assertEqual(outcome(many.fetch_capped(wanted, cap, b)),
                                                 outcome(one.fetch_capped(wanted, cap, a)))
                                self.assertEqual(a["sampled"], b["sampled"])
                                self.assertEqual(a["returned_alignments"], b["returned_alignments"])
                        self.assertEqual([x.SerializeToString() for x in many.fetch(wanted)],
                                         [x.SerializeToString() for x in one.fetch(wanted)])
                    stats = many.cache_stats["gams"]
                    self.assertEqual(len(stats), 3)
                    if mode != "off":
                        self.assertTrue(all(s["record_fetches"] > 0 and s["fallback_fetches"] == 0 for s in stats))

    def test_pairs(self):
        self.assertEqual(gam_pairs("a.gam"), [("a.gam", None)])
        self.assertEqual(gam_pairs(["a.gam", "b.gam"], ["a.gai", "b.gai"]), [("a.gam", "a.gai"), ("b.gam", "b.gai")])
        with self.assertRaisesRegex(ValueError, "twice"):
            gam_pairs(["a.gam", "./a.gam"])
        with self.assertRaisesRegex(ValueError, "one per GAM"):
            gam_pairs(["a.gam", "b.gam"], "a.gai")

    def test_one_gam_is_recorded_as_a_path(self):
        parser = make_parser()
        args = parser.parse_args(["build", "--gam", "a.gam", "--index", "a.gai", "--output", "o", "--nodes", "n",
                                  "--graph-index", "g", "--haplotypes", "82", "--snv-min-af", ".1",
                                  "--indel-min-af", ".1", "--snv-output", "s", "--indel-output", "i"])
        self.assertEqual((args.gam, args.index), ("a.gam", "a.gai"))
        args = parser.parse_args(["discover", "--gam", "a.gam", "b.gam", "--output", "o", "--graph-index", "g"])
        self.assertEqual((args.gam, args.index), (["a.gam", "b.gam"], None))


@needs_module
class DiscoveryTest(unittest.TestCase):
    def test_counts_and_outputs_equal_the_single_gam(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sequences, records = synthetic_world(12, 300)
            single = write_grouped_gam(root / "single.gam", records, 5)
            split = [write_grouped_gam(root / f"part{k}.gam", records[k::2], 4) for k in range(2)]
            graph = graph_fixture(root / "graph.sqlite", [(n, s, 1) for n, s in sequences.items()])
            one = discovery.native_counts(single, None, 2, str(graph), 5, 50)
            many = discovery.native_counts([str(p) for p in split], None, 2, str(graph), 5, 50)
            table = lambda c: {n: (p, q, m) for n, p, q, m in zip(*(x.tolist() for x in c[:4]))}  # noqa: E731
            self.assertEqual(table(many), table(one))
            self.assertEqual(many[4]["alignments"], len(records))
            outputs = {}
            for name, gam in (("one", str(single)), ("many", [str(p) for p in split])):
                out = root / name
                args = type("Args", (), dict(gam=gam, output=str(out), min_mapq=5, node_alt=0.05,
                                             graph_index=str(graph), processes=2))()
                with redirect_stdout(io.StringIO()):
                    discovery.discover(args)
                outputs[name] = ((out / "target_nodes.txt").read_text(),
                                 json.loads((out / "node_stats.json").read_text()))
                report = json.loads((out / "discovery_report.json").read_text())
                self.assertEqual(report["gam"], str(single.resolve()) if name == "one"
                                 else [str(p.resolve()) for p in split])
            self.assertEqual(outputs["many"], outputs["one"])


class BuildTest(unittest.TestCase):
    def test_builds_equal_the_single_gam(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sequences, records = random_chain_world(20261013, nodes=80, records=900)
            (root / "inputs").mkdir()
            single, split = split_world(root / "inputs", records, parts=2)
            (root / "sorted").mkdir()
            ordered, _ = split_world(root / "sorted", records, parts=2, order=True)
            rng = random.Random(7)
            graph = graph_fixture(root / "inputs" / "graph.sqlite",
                                  [(n, s, rng.choice([1, 5, 82, 300])) for n, s in sequences.items()])
            nodes = root / "inputs" / "nodes.txt"
            nodes.write_text("".join(f"{n}\n" for n in sorted(sequences)))
            options = dict(nodes=str(nodes), graph_index=str(graph), haplotypes=82, batch_nodes=7, max_node_span=40,
                           max_node_reads=6, rows=20, width=21, max_indel_len=10, shard_size=3, debug_rows=True,
                           min_variants=2)
            with redirect_stdout(io.StringIO()), patch("sys.stderr", io.StringIO()):
                build(split_args(root, "one", gam=str(single), **options))
                build(split_args(root, "ordered", gam=str(ordered), **options))
                manifest = build(split_args(root, "many", gam=[str(p) for p in split], **options))
            self.assertEqual(compare(root / "one", root / "many", masks=GAM_ARGUMENTS), [])
            # A single GAM in gamsort order: the same shards and labels; what records file order only differs:
            # the audit records' record_index (a record's place in its batch's fetch) and the debug-rows
            # selection_audit list (the windows in fetch order).
            differences = compare(root / "ordered", root / "many", masks=GAM_ARGUMENTS)
            self.assertTrue(all(d.endswith(("unsupported_events.ndjson", "variant_summary.ndjson"))
                                for d in differences), differences)

            def rows(run, name, stream):
                out = []
                for line in (root / run / name / stream).read_text().splitlines():
                    row = json.loads(line)
                    row.pop("record_index", None)
                    if "selection_audit" in row:
                        row["selection_audit"].sort(key=lambda w: json.dumps(w, sort_keys=True))
                    out.append(json.dumps(row, sort_keys=True))
                return sorted(out) if stream == "unsupported_events.ndjson" else out
            for name in ("shared", "INDEL", "SNV"):
                for stream in ("unsupported_events.ndjson", "variant_summary.ndjson"):
                    if (root / "many" / name / stream).exists():
                        self.assertEqual(rows("ordered", name, stream), rows("many", name, stream), (name, stream))
            self.assertGreater(manifest["tensors"], 0)
            self.assertEqual(manifest["arguments"]["gam"], [str(p) for p in split])
            self.assertEqual(len(manifest["gam_group_cache"]["gams"]), 2)


class OrchestrateTest(unittest.TestCase):
    def test_prepare_and_run_equal_the_single_gam(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sequences, records = random_chain_world(20261014, nodes=50, records=500)
            (root / "inputs").mkdir()
            single, split = split_world(root / "inputs", records, parts=2)
            graph = graph_fixture(root / "inputs" / "graph.sqlite", [(n, s, 40) for n, s in sequences.items()])
            nodes = root / "inputs" / "nodes.txt"
            nodes.write_text("".join(f"{n}\n" for n in sorted(sequences)))
            configs = {}
            for name, gams in (("one", [single]), ("many", split)):
                argv = ["prepare", "--root", str(root / name), "--gam", *map(str, gams), "--nodes", str(nodes),
                        "--graph-index", str(graph), "--haplotypes", "82", "--tasks", "3", "--processes", "2",
                        "--gam-cache-mb", "1", "--batch-nodes", "5", "--shard-size", "2", "--min-variants", "2",
                        "--max-node-reads", "6", "--merge-shard-size", "0", "--snv-min-af", ".05", "--indel-min-af", ".05"]
                with redirect_stdout(io.StringIO()), patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")):
                    orchestrate.main(argv)
                    orchestrate.run(root / name)
                configs[name] = json.loads((root / name / "config.json").read_text())
            inputs = configs["many"]["inputs"]
            self.assertEqual(orchestrate.config_gams(configs["many"]), [(str(p), str(p) + ".gai") for p in split])
            self.assertNotIn("gam_2", configs["one"]["inputs"])
            self.assertEqual((inputs["gam_2"]["path"], inputs["index_2"]["path"]), (str(split[1]), str(split[1]) + ".gai"))
            self.assertEqual(configs["many"]["record_index"], [None, None])
            self.assertEqual(compare(root / "one", root / "many", masks=GAM_ARGUMENTS,
                                     anchor_a=root / "one", anchor_b=root / "many"), [])
            status = json.loads((root / "many" / "status.json").read_text())
            self.assertEqual(status["status"], "complete")
            self.assertGreater(status["tensors"], 0)


if __name__ == "__main__":
    unittest.main()

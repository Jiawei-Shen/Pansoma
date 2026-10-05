"""Whole `run build`s on a vg-style multi-record GAM: the reader of 0886f70 (gam_reader_oracle, patched into
build) against the GAI walk on BgzfReader, the record path from <gam>.gri and the record path in memory. Auto
batches, a small --max-batch-alignments (split batches), a small --max-node-reads (capped nodes) and --debug-rows;
tools.compare_runs must find no difference."""
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from . import gam_reader_oracle as oracle
from .fixtures import graph_fixture, split_args
from .golden import random_chain_world
from .test_record_index import World, struct_field
from .. import build as build_module, gam_record_index
from ..build import build
from ..gam_record_index import ENVIRONMENT
from ..tools.compare_runs import compare


class WholeBuildTest(unittest.TestCase):
    def test_builds_with_every_reader_are_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sequences, records = random_chain_world(20261005, nodes=90, records=1000)
            rng = random.Random(5)
            raws = [a.SerializeToString(deterministic=True) + (struct_field(rng) if rng.random() < 0.7 else b"")
                    for a in records]
            (root / "inputs").mkdir()
            world = World(root / "inputs", "chain", 9, raws=raws, foreign=0.05)
            graph = graph_fixture(root / "inputs" / "graph.sqlite",
                                  [(n, s, rng.choice([1, 5, 90, 300, 464, 1000])) for n, s in sequences.items()])
            nodes = root / "inputs" / "nodes.txt"
            nodes.write_text("".join(f"{n}\n" for n in sorted(sequences)))
            with patch("sys.stderr", io.StringIO()):
                gam_record_index.build(world.path, processes=2, force=True)
            options = dict(gam=str(world.path), nodes=str(nodes), graph_index=str(graph), haplotypes=464,
                           batch_nodes="auto", max_node_span=12, max_batch_alignments=20, max_node_reads=6,
                           rows=20, width=21, max_indel_len=10, shard_size=3, debug_rows=True, min_variants=2)
            outputs, caches = {}, {}
            for name, mode in (("oracle", "off"), ("walk", "off"), ("file", "require"), ("memory", "memory")):
                args = split_args(root, name, **options)
                reader = oracle.IndexedGam if name == "oracle" else build_module.IndexedGam
                with patch.dict(os.environ, {ENVIRONMENT: mode}), patch.object(build_module, "IndexedGam", reader), \
                        redirect_stdout(io.StringIO()), patch("sys.stderr", io.StringIO()):
                    caches[name] = build(args)["gam_group_cache"]
                outputs[name] = root / name
            for name in ("walk", "file", "memory"):
                self.assertEqual(compare(outputs["oracle"], outputs[name]), [], name)
            for name in ("file", "memory"):
                self.assertGreater(caches[name]["record_fetches"], 0)
                self.assertEqual(caches[name]["fallback_fetches"], 0)
            self.assertNotIn("record_fetches", caches["oracle"])
            timing = [json.loads(line) for line in
                      (root / "oracle/shared/batch_timing.ndjson").read_text().splitlines()]
            self.assertTrue(any(row["batch_plan"].get("split") for row in timing))
            self.assertTrue((root / "oracle/shared/downsampled_nodes.tsv").exists())
            summary = (root / "oracle/SNV/variant_summary.ndjson").read_text().splitlines()
            self.assertGreater(len(summary), 10)


if __name__ == "__main__":
    unittest.main()

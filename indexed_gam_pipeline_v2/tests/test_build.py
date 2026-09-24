"""End-to-end builds on synthetic data: output contract, split mode, limits and the auditor."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from fixtures import af_gam, build_args, graph_fixture, tiny_gam  # noqa: F401  (sets sys.path)
from indexed_gam_pipeline_v2 import build as build_module
from indexed_gam_pipeline_v2.build import build
from indexed_gam_pipeline_v2.candidates import BASES, CHANNELS, decode_alignment, encode_count
from indexed_gam_pipeline_v2.orchestrate import validate_shards
from indexed_gam_pipeline_v2.validate_examples import validate


def quiet_build(args):
    with redirect_stdout(io.StringIO()):
        return build(args)


def files(folder, names=("variant_summary.ndjson", "filtered_candidates.ndjson", "unsupported_events.ndjson")):
    folder = Path(folder)
    return {p.name: p.read_bytes() for p in list(folder.glob("shard_*_data.npy")) + [folder / n for n in names if (folder / n).exists()]}


class BuildTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.gam, _ = tiny_gam(self.root)
        self.nodes = self.root / "nodes.txt"
        self.nodes.write_text("10\n20\n30\n1000\n")
        self.graph = graph_fixture(self.root / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 1000)])
        self.base = dict(gam=str(self.gam), nodes=str(self.nodes), graph_index=str(self.graph))

    def tearDown(self):
        self.tmp.cleanup()

    def args(self, **overrides):
        return build_args(**dict(self.base, **overrides))

    def test_shards_summary_manifest_contract_and_audit(self):
        args = self.args(output=str(self.root / "out"), debug_rows=True)
        manifest = quiet_build(args)
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual((manifest["tensors"], manifest["shards"], manifest["shape"]), (4, 2, [8, 200, 101]))
        summary = [json.loads(s) for s in (self.root / "out/variant_summary.ndjson").read_text().splitlines()]
        self.assertEqual(len(summary), 4)
        for i, meta in enumerate(summary):
            self.assertEqual((meta["shard_index"], meta["index_within_shard"]), (i // 2, i % 2))
            x = np.load(self.root / f"out/shard_{i // 2:05d}_data.npy")[i % 2]
            self.assertEqual((list(x.shape), str(x.dtype)), (manifest["shape"], "int8"))
            self.assertEqual(int(x[6].max()), encode_count(331))
            self.assertTrue(np.all(x[7][x[0] != 0] == 1) and not x[7][x[0] == 0].any())  # forward-only fixture
            self.assertEqual(meta["channels"], CHANNELS)
            self.assertEqual(meta["coverage"], meta["alt_count"] + meta["ref_count"] + meta["other_count"])
            self.assertEqual(meta["parameters"]["min_af"], 0.05)
        self.assertEqual(json.loads((self.root / "out/manifest.json").read_text())["tensors"], 4)
        self.assertEqual(len((self.root / "out/batch_timing.ndjson").read_text().splitlines()), 3)  # [10,20] [30] [1000]
        self.assertEqual(validate_shards(self.root / "out", 2)["tensors"], 4)
        audit = validate(self.root / "out", str(self.gam), None, str(self.graph))
        self.assertTrue(audit["passed"])
        self.assertEqual(audit["examples"], 4)
        # The auditor must catch a corrupted graph-reference channel.
        shard = self.root / "out/shard_00000_data.npy"
        corrupted = np.load(shard)
        corrupted[0, 5, 0, 49] = BASES["C"]
        np.save(shard, corrupted)
        with self.assertRaisesRegex(ValueError, "graph base"):
            validate(self.root / "out", str(self.gam), None, str(self.graph))

    def test_target_node_restriction_does_not_change_outputs(self):
        # The Python decoder, patched (tests/test_native_decoder.py covers the native one).
        args = self.args(batch_nodes=1, max_node_span=10000, debug_rows=True, rows=4, decoder="python")
        unrestricted = lambda *a, **kw: decode_alignment(*a, **dict(kw, target_nodes=None))  # noqa: E731
        args.output = str(self.root / "baseline")
        with patch.object(build_module, "decode_alignment", unrestricted):
            quiet_build(args)
        args.output = str(self.root / "restricted")
        quiet_build(args)
        self.assertEqual(files(self.root / "baseline"), files(self.root / "restricted"))

    def test_max_tensors_and_missing_graph_node(self):
        args = self.args(output=str(self.root / "limited"), max_tensors=3)
        manifest = quiet_build(args)
        self.assertEqual((manifest["tensors"], manifest["shards"]), (3, 2))
        self.assertEqual(validate_shards(self.root / "limited", 2)["last_shard_tensors"], 1)
        graph_fixture(self.root / "partial.sqlite", [(n, "AAAAAA", 1) for n in (10, 20)])
        args = self.args(output=str(self.root / "missing"), graph_index=str(self.root / "partial.sqlite"))
        with self.assertRaisesRegex(ValueError, "missing from GBZ graph index"):
            quiet_build(args)

    def test_argument_validation(self):
        for overrides, message in ((dict(max_indel_len=51), "max-indel-len"),
                                   (dict(width=10, max_indel_len=20), "width"),
                                   (dict(snv_output="x"), "Split output needs"),
                                   (dict(snv_output="x", indel_output="y", snv_min_af=.1, indel_min_af=.1, variant_type="snp"), "variant-type all"),
                                   (dict(snv_output="o", indel_output="o", snv_min_af=.1, indel_min_af=.1), "distinct")):
            args = self.args(**dict(dict(output="o"), **overrides))
            with self.assertRaisesRegex(ValueError, message):
                quiet_build(args)


class SplitOutputTest(unittest.TestCase):
    def test_one_shared_pass_equals_two_independent_builds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = af_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n40\n50\n")
            common = dict(gam=str(gam), nodes=str(nodes), graph_index=str(graph), batch_nodes=5,
                          shard_size=1, min_variants=3, min_af=.05)
            split = build_args(**dict(common, output=str(root / "shared"), snv_output=str(root / "SNV"),
                                      indel_output=str(root / "INDEL"), snv_min_af=.06, indel_min_af=.08,
                                      decoder="python"))  # counts calls of the (patched) Python decoder
            fetch_calls, decode_calls = [], []
            reader_fetch = build_module.IndexedGam.fetch

            def counting_fetch(self, *a, **kw):
                fetch_calls.append(1)
                return reader_fetch(self, *a, **kw)

            def counting_decode(*a, **kw):
                decode_calls.append(1)
                return decode_alignment(*a, **kw)

            with patch.object(build_module.IndexedGam, "fetch", counting_fetch), \
                 patch.object(build_module, "decode_alignment", counting_decode):
                manifest = quiet_build(split)
            self.assertEqual((len(fetch_calls), len(decode_calls)), (1, 250))
            self.assertEqual(manifest["tensors_by_type"], dict(SNV=1, INDEL=2))
            self.assertEqual(manifest["output_layout"], "split")
            self.assertFalse(list((root / "shared").glob("shard_*")))
            snv = [json.loads(l) for l in (root / "SNV/variant_summary.ndjson").read_text().splitlines()]
            indel = [json.loads(l) for l in (root / "INDEL/variant_summary.ndjson").read_text().splitlines()]
            self.assertEqual([(m["node_id"], m["event_type"], m["af"]) for m in snv], [(10, "SNP", .06)])
            self.assertEqual([(m["node_id"], m["event_type"], m["af"]) for m in indel], [(30, "DEL", .08), (50, "INS", .08)])
            self.assertEqual(snv[0]["parameters"]["variant_type"], "snp")
            for kind, variant_type, af in (("SNV", "snp", .06), ("INDEL", "indel", .08)):
                single = build_args(**dict(common, output=str(root / f"only_{kind}"), variant_type=variant_type, min_af=af))
                quiet_build(single)
                # Shards and summaries must match byte for byte. (Audit streams differ by design: a
                # --variant-type snp build logs INDEL candidates as rejected, split mode routes them.)
                self.assertEqual(files(root / kind, names=("variant_summary.ndjson",)),
                                 files(root / f"only_{kind}", names=("variant_summary.ndjson",)), kind)
                self.assertEqual(validate_shards(root / kind, 1)["tensors"], len(snv) if kind == "SNV" else len(indel))
            # DEL 3/50 = 0.06 fails the 0.08 INDEL threshold but would pass the SNV one.
            filtered = [json.loads(l) for l in (root / "INDEL/filtered_candidates.ndjson").read_text().splitlines()]
            self.assertIn((20, ["min_af"]), [(m["node_id"], m["reasons"]) for m in filtered])
            limited = build_args(**dict(common, output=str(root / "lim"), snv_output=str(root / "lim_SNV"),
                                        indel_output=str(root / "lim_INDEL"), snv_min_af=.06, indel_min_af=.08, max_tensors=1))
            self.assertEqual(quiet_build(limited)["tensors"], 1)
            empty = build_args(**dict(common, output=str(root / "empty"), snv_output=str(root / "empty_SNV"),
                                      indel_output=str(root / "empty_INDEL"), snv_min_af=1, indel_min_af=1))
            self.assertEqual(quiet_build(empty)["tensors_by_type"], dict(SNV=0, INDEL=0))


if __name__ == "__main__":
    unittest.main()

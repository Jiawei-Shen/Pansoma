"""The native (C++) decoder: identical to candidates.decode_alignment; safe selection and fallback.

Equivalence and build tests run when the module is built (python -m indexed_gam_pipeline_v2.native
compile); the selection / fallback tests that need no module always run.
"""
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from fixtures import af_gam, build_args, graph_fixture, spec_alignment  # noqa: F401  (sets sys.path)
from indexed_gam_pipeline_v2 import native
from indexed_gam_pipeline_v2.build import build
from indexed_gam_pipeline_v2.candidates import decode_alignment

MODULE, INFO = native.load()
needs_module = unittest.skipIf(MODULE is None, f"native decoder unavailable: {INFO.get('reason')}")
VOLATILE = ("decoder", "arguments", "timing", "graph_index_performance")


def outputs(folder):
    """{relative path: bytes} of data files, {relative path: manifest without volatile keys} of manifests."""
    folder = Path(folder)
    result = {}
    for p in sorted(folder.rglob("*")):
        if p.name == "batch_timing.ndjson" or not p.is_file():
            continue
        if p.name == "manifest.json":
            manifest = json.loads(p.read_text().replace(str(folder), "<OUT>"))
            result[str(p.relative_to(folder))] = {k: v for k, v in manifest.items() if k not in VOLATILE}
        else:
            result[str(p.relative_to(folder))] = p.read_bytes()
    return result


@needs_module
class EquivalenceTest(unittest.TestCase):
    def test_random_records_are_decoded_identically(self):
        """Scattered and repeat-chain records, both orientations, every edit kind, N bases, missing
        or malformed qualities, max_indel 2/3/50, with and without target nodes and left-normalization:
        every column (with value types), visit, observation, move, unsupported event and error type."""
        rng = random.Random(7)
        for index in range(30000):
            case = native.synthetic_case(rng, index)
            expected = native._outcome(decode_alignment, case)
            actual = native._outcome(lambda *a, **kw: native.decode_native(MODULE, *a, **kw), case)
            if expected != actual:
                self.fail(f"synthetic record {index} differs: sequences {case[1]}, alignment {case[0]}")

    def test_column_array_behaves_like_the_column_list(self):
        sequences = {1: "ACGTACGTAC" * 30, 2: "GGGCCC"}  # several 128-column blocks
        a = spec_alignment([(1, 0, False, [(150, 150, ""), (0, 2, "TT"), (150, 150, "")]),
                            (2, 0, True, [(6, 6, "")])], sequences)
        reference = decode_alignment(a, sequences)[0].columns
        columns = native.decode_native(MODULE, a, sequences)[0].columns
        self.assertIsInstance(columns, native.ColumnArray)
        self.assertEqual(len(columns), len(reference))
        for key in (0, 127, 128, 200, -1, -len(reference)):
            self.assertEqual(columns[key], reference[key])
        for key in (slice(None), slice(0, 0), slice(120, 260), slice(-10, None), slice(5, 300, 7), slice(400, 500),
                    slice(300, 100)):
            self.assertEqual(columns[key], reference[key])
        self.assertEqual(list(columns), reference)
        self.assertTrue(columns == reference and reference == columns)
        for key in (len(reference), -len(reference) - 1):
            with self.assertRaises(IndexError):
                columns[key]

    def test_builds_are_identical(self):
        """Split-mode builds with --debug-rows: every shard, summary and audit stream byte for byte."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam, _ = af_gam(root)
            graph = graph_fixture(root / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
            nodes = root / "nodes.txt"
            nodes.write_text("10\n20\n30\n40\n50\n")
            results = {}
            for decoder in ("python", "native"):
                out = root / decoder
                args = build_args(gam=str(gam), nodes=str(nodes), graph_index=str(graph), batch_nodes=2,
                                  shard_size=1, min_variants=3, min_af=.05, debug_rows=True, output=str(out / "shared"),
                                  snv_output=str(out / "SNV"), indel_output=str(out / "INDEL"), snv_min_af=.06,
                                  indel_min_af=.08, decoder=decoder)
                with redirect_stdout(io.StringIO()):
                    manifest = build(args)
                self.assertEqual(manifest["decoder"]["used"], decoder)
                results[decoder] = outputs(out)
            self.assertEqual(results["python"], results["native"])
            self.assertEqual(manifest["decoder"]["native_record_fallbacks"], 0)
            self.assertEqual(manifest["tensors_by_type"], dict(SNV=1, INDEL=2))


class SelectionTest(unittest.TestCase):
    def setUp(self):
        environment = patch.dict(os.environ)
        environment.start()
        self.addCleanup(environment.stop)
        os.environ.pop("PANSOMA_DECODER", None)

    def test_explicit_python_environment_override_and_bad_choice(self):
        decode, info = native.select_decoder("python")
        self.assertIs(decode, decode_alignment)
        self.assertEqual(info, dict(requested="python", used="python"))
        os.environ["PANSOMA_DECODER"] = "python"
        self.assertEqual(native.select_decoder("auto")[1], dict(requested="auto", environment="PANSOMA_DECODER=python",
                                                               used="python"))
        with self.assertRaisesRegex(ValueError, "--decoder"):
            native.select_decoder("fast")

    def test_unavailable_module_falls_back_or_fails_loudly(self):
        with patch.object(native, "_LOADED", (None, dict(reason="not built"))):
            decode, info = native.select_decoder("auto")
            self.assertIs(decode, decode_alignment)
            self.assertEqual(info, dict(requested="auto", used="python", reason="not built"))
            with self.assertRaisesRegex(ValueError, "not built"):
                native.select_decoder("native")

    def test_import_failure_is_reported(self):
        with patch.object(native.importlib, "import_module", side_effect=ImportError("no such module")):
            module, info = native._load()
        self.assertIsNone(module)
        self.assertIn("not built", info["reason"])

    @needs_module
    def test_stale_build_or_failed_self_test_is_not_used(self):
        with patch.object(native, "sha256_file", return_value="0" * 64):
            module, info = native._load()
        self.assertIsNone(module)
        self.assertIn("stale", info["reason"])
        with patch.object(native, "self_test", return_value="differs on synthetic record 3"):
            module, info = native._load()
        self.assertIsNone(module)
        self.assertIn("record 3", info["reason"])

    def test_a_record_the_native_decoder_fails_on_is_decoded_in_python(self):
        class Broken:
            @staticmethod
            def decode(*args):
                raise RuntimeError("native failure")

        decoder = native.NativeDecoder(Broken)
        sequences = {1: "ACGTAC"}
        a = spec_alignment([(1, 0, False, [(2, 2, ""), (1, 1, "T"), (3, 3, "")])], sequences)
        self.assertEqual(native.snapshot(decoder(a, sequences)), native.snapshot(decode_alignment(a, sequences)))
        self.assertEqual(decoder.fallbacks, 1)
        a.quality = b"\x1e"  # malformed: the Python decoder's own error is raised
        with self.assertRaisesRegex(ValueError, "quality length"):
            decoder(a, sequences)


if __name__ == "__main__":
    unittest.main()

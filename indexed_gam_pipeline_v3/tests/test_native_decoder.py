"""The native (C++) decoder: identical to candidates.decode_alignment; safe selection and fallback.

Equivalence and build tests run when the module is built (python -m indexed_gam_pipeline_v3.native
compile); the selection / fallback tests that need no module always run. A module that is built
but does not load is a failure, not a skip.
"""
from contextlib import redirect_stdout
import io
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from .fixtures import af_gam, build_args, graph_fixture, spec_alignment  # noqa: F401
from .. import native
from ..build import build
from ..candidates import decode_alignment
from ..common import sha256_file
from ..tools.compare_runs import compare, fingerprint

MODULE, INFO = native.load()
needs_module = unittest.skipIf(MODULE is None, f"native decoder unavailable: {INFO.get('reason')}")


class BuiltModuleTest(unittest.TestCase):
    def test_a_compiled_module_next_to_native_py_must_load(self):
        """Catches a stale build (compiled from another fastdecode.cpp), a module built for another
        interpreter and a module of another package loaded into this process."""
        built = sorted(native.HERE.glob(native.MODULE + "*.so"))
        if not built:
            self.skipTest(f"no compiled {native.MODULE} next to native.py (python -m {native.__package__}.native compile)")
        self.assertIsNotNone(MODULE, f"{built[0].name} is present but not used: {INFO.get('reason')}")
        self.assertEqual(Path(MODULE.__file__).resolve(), built[0].resolve())
        self.assertEqual(MODULE.__name__, f"{native.__package__}.{native.MODULE}")
        self.assertEqual(INFO, dict(module=built[0].name, source_sha256=sha256_file(native.SOURCE),
                                    self_test_records=native.SELF_TEST_CASES))
        self.assertFalse(list(native.HERE.glob(native.MODULE + ".build.json")))  # compile writes no build.json


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

    def test_column_array_keeps_only_recent_blocks(self):
        sequences = {1: "ACGTACGTAC" * 300}  # 3000 columns = 24 blocks
        a = spec_alignment([(1, 0, False, [(3000, 3000, "")])], sequences)
        reference = decode_alignment(a, sequences)[0].columns
        columns = native.decode_native(MODULE, a, sequences)[0].columns
        self.assertEqual(list(columns), reference)  # a full walk builds every block once
        self.assertEqual(len(columns.blocks), native.CACHED_BLOCKS)
        self.assertEqual(list(columns.blocks), list(range(24 - native.CACHED_BLOCKS, 24)))
        self.assertEqual(columns[0:300], reference[0:300])  # evicted blocks come back equal
        self.assertEqual(list(columns.blocks)[-3:], [0, 1, 2])
        recent = columns.blocks[0]
        columns[5]  # a hit refreshes the block, no rebuild
        self.assertEqual(list(columns.blocks)[-1], 0)
        self.assertIs(columns.blocks[0], recent)
        self.assertEqual(len(columns.blocks), native.CACHED_BLOCKS)

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
                results[decoder] = out
            # every file, byte for byte or (manifests, batch timing) without the decoder and timing keys
            self.assertEqual(compare(results["python"], results["native"]), [])
            self.assertLessEqual({"shared/manifest.json", "shared/batch_timing.ndjson", "SNV/shard_00000_data.npy",
                                  "INDEL/shard_00001_data.npy", "INDEL/variant_summary.ndjson"},
                                 set(fingerprint(results["native"])))
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

"""This package against the goldens recorded from v2 (tests/golden.py, golden_hashes.json), and the
normalization of tools.compare_runs that the goldens rely on."""
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np

from . import golden
from .fixtures import af_gam, graph_fixture, split_args
from ..build import build
from ..tools.compare_runs import compare

PACKAGE = __package__.split(".")[0]
NATIVE, REASON = golden.native_module(PACKAGE)


@unittest.skipIf(os.environ.get("PANSOMA_DECODER"), "PANSOMA_DECODER is set; the golden check runs both decoders")
class GoldenTest(unittest.TestCase):
    """Every golden case rebuilt by this package (in subprocesses) reproduces v2's recorded fingerprints."""

    @classmethod
    def setUpClass(cls):
        cls.recorded = golden.load()
        decoders = ["python"] + (["native"] if NATIVE else [])
        if not NATIVE:
            print(f"\ngolden native builds skipped: {REASON}")
        cls.work = tempfile.mkdtemp(prefix="golden-test-")
        cls.results = golden.run_jobs(PACKAGE, golden.jobs(decoders=decoders), cls.work)
        cls.problems = golden.check_results(cls.recorded, cls.results)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def assert_matches(self, labels):
        self.assertTrue(labels)
        found = [f"{label}: {line}" for label in labels for line in self.problems.get(label, [])]
        if found:
            drift = golden.environment_drift(self.recorded)
            self.fail("golden mismatch:\n" + "\n".join(found)
                      + ("\nenvironment differs from the recording:\n" + "\n".join(drift) if drift else ""))

    def test_build_cases_with_the_python_decoder(self):
        self.assert_matches([f"{case}-python" for case in golden.BUILD_CASES])

    @unittest.skipUnless(NATIVE, REASON)
    def test_build_cases_with_the_native_decoder(self):
        self.assert_matches([f"{case}-native" for case in golden.BUILD_CASES])

    def test_orchestrated_run(self):
        """O1 under --decoder auto: native when built, and Python under PANSOMA_DECODER=python."""
        self.assert_matches([label for label in self.results if label.startswith("O1-")])

    def test_random_world_is_independent_of_the_hash_seed(self):
        self.assert_matches([label for label in self.results if "-seed" in label])

    def test_output_constants(self):
        self.assertEqual(golden.compare_constants(self.recorded["constants"], golden.constants(PACKAGE)), [])


class CompareRunsTest(unittest.TestCase):
    """A one-byte change of a shard, a summary line or a labels.npy is a difference; timing, created and
    decoder values are not."""

    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="compare-runs-"))
        gam, _ = af_gam(self.root)
        graph = graph_fixture(self.root / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
        (self.root / "nodes.txt").write_text("10\n20\n30\n40\n50\n")
        with redirect_stdout(io.StringIO()):
            build(split_args(self.root, "a", gam=str(gam), graph_index=str(graph), nodes=str(self.root / "nodes.txt"),
                             min_variants=3, snv_min_af=0.06, indel_min_af=0.08, shard_size=1))
        np.save(self.root / "a" / "SNV" / "chr1_shard_00000_labels.npy", np.array([0, 1, -1, 2], dtype=np.int8))
        shutil.copytree(self.root / "a", self.root / "b")
        self.a, self.b = self.root / "a", self.root / "b"

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def differences(self):
        return compare(self.a, self.b, anchor_a=self.root, anchor_b=self.root)

    def flip_last_byte(self, path):
        data = bytearray(path.read_bytes())
        data[-1] ^= 1
        path.write_bytes(bytes(data))

    def edit_json(self, path, change):
        value = json.loads(path.read_text())
        change(value)
        path.write_text(json.dumps(value, indent=2) + "\n")

    def test_identical_copies_have_no_differences(self):
        self.assertEqual(self.differences(), [])

    def test_one_byte_in_a_shard_is_a_difference(self):
        self.flip_last_byte(self.b / "INDEL" / "shard_00001_data.npy")
        self.assertEqual(self.differences(), ["differs: INDEL/shard_00001_data.npy"])

    def test_one_summary_line_is_a_difference(self):
        path = self.b / "INDEL" / "variant_summary.ndjson"
        lines = path.read_text().splitlines(keepends=True)
        record = json.loads(lines[-1])
        record["alt_count"] += 1
        path.write_text("".join(lines[:-1]) + json.dumps(record) + "\n")
        self.assertEqual(self.differences(), ["differs: INDEL/variant_summary.ndjson"])

    def test_one_byte_in_a_labels_file_is_a_difference(self):
        self.flip_last_byte(self.b / "SNV" / "chr1_shard_00000_labels.npy")
        self.assertEqual(self.differences(), ["differs: SNV/chr1_shard_00000_labels.npy"])

    def test_timing_created_and_decoder_values_are_not_differences(self):
        def volatile(m):
            m["timing"] = {k: v + 1 for k, v in m["timing"].items()}
            m["decoder"] = dict(used="elsewhere")
            m["arguments"]["decoder"] = "native" if m["arguments"]["decoder"] != "native" else "python"
            m["graph_index_performance"] = dict(queries=-1)
            m["created"] = "2000-01-01T00:00:00"
        for kind in ("shared", "SNV", "INDEL"):
            self.edit_json(self.b / kind / "manifest.json", volatile)
        path = self.b / "shared" / "batch_timing.ndjson"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        for row in rows:
            row["elapsed_seconds"] += 5
            row["cumulative_stage_seconds"] = {k: v * 2 for k, v in row["cumulative_stage_seconds"].items()}
        path.write_text("".join(json.dumps(row) + "\n" for row in rows))
        for ignored in ("memory.ndjson", "logs/task_0000.log", "source/x.py", "queue_status.json", "run.sh"):
            (self.b / ignored).parent.mkdir(parents=True, exist_ok=True)
            (self.b / ignored).write_text("anything")
        self.assertEqual(self.differences(), [])

    def test_manifest_values_key_order_and_extra_files_are_differences(self):
        self.edit_json(self.b / "SNV" / "manifest.json", lambda m: m["parameters"].update(min_af=0.07))
        path = self.b / "INDEL" / "manifest.json"
        value = json.loads(path.read_text())
        path.write_text(json.dumps(dict(reversed(list(value.items())))))
        (self.b / "shared" / "extra.txt").write_text("x")
        batch = self.b / "shared" / "batch_timing.ndjson"
        batch.write_text(batch.read_text().replace('"alignments": ', '"alignments": 1', 1))
        self.assertEqual(self.differences(), [
            "only in B: shared/extra.txt", "differs: INDEL/manifest.json (key order)",
            "differs: SNV/manifest.json (keys: parameters)", "differs: shared/batch_timing.ndjson"])


if __name__ == "__main__":
    unittest.main()

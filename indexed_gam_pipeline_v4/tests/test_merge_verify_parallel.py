"""tensor_postprocessing.merge_shards verification in the pool: the spot-check lines drawn in the parent are the ones
rng.sample over each summary's lines takes, the spot count is the serial one, and a damaged merged shard, summary or
source tensor stops the merge with the error of the serial verification (kept here, verbatim, as the oracle)."""
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest import mock

import numpy as np

from ..tensor_postprocessing import merge_shards
from ..tensor_postprocessing.merge_shards import merge, spot_picks, verify_shard, verify_summary
from .test_tensor_postprocessing_merge import fake_run, write_chr_index

SEED = 20260923


def serial_spot_check(directories, results, source_paths, spots, seed):
    """merge_shards' serial spot check, verbatim (spot_check_sources)."""
    rng = random.Random(seed)
    checked = 0
    for group, result in sorted(results.items()):
        directory = directories[result["dataset"]]
        with (directory / result["summary"]).open() as stream:
            lines = stream.readlines()
        for line in rng.sample(lines, min(spots, len(lines))):
            r = json.loads(line)
            out = np.load(directory / r["shard_file"], mmap_mode="r")[r["index_within_shard"]]
            folder = Path(source_paths[r["source_task"]])
            source = np.load(folder / f"shard_{r['source_shard_index']:05d}_data.npy",
                             mmap_mode="r")[r["source_index_within_shard"]]
            if not np.array_equal(out, source):
                raise ValueError(f"Spot check failed: {group} {r['shard_file']}[{r['index_within_shard']}]")
            checked += 1
    return checked


def serial_verify(written, work, outputs, sources_of, spots, seed):
    """merge's serial step 2 (summaries and totals per kind, shards, spot checks); spot_checked."""
    jobs = []
    for kind, (reference, results) in written.items():
        for group, result in results.items():
            verify_summary(work[kind][result["dataset"]], group, result)
            jobs += [(work[kind][result["dataset"]] / s["file"], s["tensors"], reference["shape"], s["sha256"])
                     for s in result["shards"]]
        total = sum(r["tensors"] for r in results.values())
        if total != outputs["tensors_by_type"][kind] or total != sum(s["tensors"] for s in sources_of[kind]):
            raise ValueError(f"{kind}: merged {total} tensors, outputs.json lists {outputs['tensors_by_type'][kind]}")
    for job in jobs:
        verify_shard(*job)
    return {kind: serial_spot_check(work[kind], results, {s["task"]: s["path"] for s in sources_of[kind]}, spots, seed)
            for kind, (reference, results) in written.items()}


def flip(path, offset):
    data = bytearray(Path(path).read_bytes())
    data[offset] ^= 1
    Path(path).write_bytes(bytes(data))


class SpotPickTest(unittest.TestCase):
    def test_picks_are_the_lines_rng_sample_takes(self):
        """Group sizes on both sides of random.sample's set/list switch (21 + 4 ** ceil(log4(3k))), groups given out
        of name order; one generator per kind."""
        sizes = (1, 2, 5, 49, 50, 51, 199, 200, 201, 1044, 1045, 1046, 1047, 4000, 100000)
        for seed in (SEED, 1, 7):
            for spots in (0, 1, 3, 50, 200):
                names = [f"chr{c}" for c in random.Random(seed).sample(range(1, 23), len(sizes))]
                results = {name: dict(tensors=n) for name, n in zip(names, sizes)}
                rng = random.Random(seed)
                expected = {}
                for group, result in sorted(results.items()):
                    lines = [f"{i}\n" for i in range(result["tensors"])]
                    expected[group] = [int(line) for line in rng.sample(lines, min(spots, len(lines)))]
                picks = spot_picks(results, spots, seed)
                self.assertEqual(list(picks.items()), list(expected.items()), (seed, spots))


class VerifyTest(unittest.TestCase):
    def merge_with(self, damage=None, spots=50, workers=2):
        """merge() of test_tensor_postprocessing_merge.fake_run with `damage(written, work, tensors)` applied after the
        copy: (spot_checked or the ValueError message, the same from serial_verify on the damaged files)."""
        original, oracle = merge_shards.write_kinds, {}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tensors, _ = fake_run(root)

            def write_kinds(sources_of, index, shard_size, work, reference_path, pool, workers):
                written = original(sources_of, index, shard_size, work, reference_path, pool, workers)
                if damage:
                    written = damage(written, work, tensors) or written
                outputs = json.loads((root / "outputs.json").read_text())
                try:
                    oracle["result"] = serial_verify(written, work, outputs, sources_of, spots, SEED)
                except (ValueError, OSError) as error:
                    oracle["result"] = f"{type(error).__name__}: {error}"
                return written

            with mock.patch.object(merge_shards, "write_kinds", write_kinds):
                try:
                    found = merge(root, write_chr_index(root / "chr.tsv"), shard_size=4, workers=workers, spots=spots,
                                  seed=SEED)["spot_checked"]
                except (ValueError, OSError) as error:
                    found = f"{type(error).__name__}: {error}"
        return found, oracle["result"]

    def test_spot_count_as_serial(self):
        for spots, workers in ((50, 2), (3, 1), (1, 4), (0, 2)):
            found, expected = self.merge_with(spots=spots, workers=workers)
            self.assertEqual(found, expected)
            self.assertEqual(list(found), ["SNV", "INDEL"])
        self.assertEqual(self.merge_with(spots=50)[0], {"SNV": 13, "INDEL": 10})  # every tensor of fake_run

    def test_damage_raises_the_serial_error(self):
        def merged_tensor(written, work, tensors):  # caught by the shard hash and by the spot check: the hash first
            flip(work["SNV"]["autosome"] / "chr1_shard_00001_data.npy", -1)

        def merged_truncated(written, work, tensors):
            path = work["INDEL"]["autosome"] / "chr2_shard_00000_data.npy"
            path.write_bytes(path.read_bytes()[:-7])

        def summary_record(written, work, tensors):
            path = work["SNV"]["autosome"] / "chr2_variant_summary.ndjson"
            lines = path.read_text().splitlines(keepends=True)
            lines[1] = json.dumps(dict(json.loads(lines[1]), af=0.5)) + "\n"
            path.write_text("".join(lines))

        def summary_order(written, work, tensors):
            path = work["INDEL"]["autosome"] / "chr2_variant_summary.ndjson"
            lines = path.read_text().splitlines(keepends=True)
            path.write_text("".join([lines[1], lines[0]] + lines[2:]))

        def source_tensor(written, work, tensors):  # only the spot check sees it
            array = np.load(tensors / "SNV" / "task_0001" / "shard_00001_data.npy", mmap_mode="r+")
            array[0, 0, 0, 0] ^= 1  # node 150: chr2's second tensor
            array.flush()
            del array

        def source_missing(written, work, tensors):
            (tensors / "INDEL" / "task_0001" / "shard_00002_data.npy").unlink()

        def totals(written, work, tensors):
            reference, results = written["INDEL"]
            return dict(written, INDEL=(reference, {g: r for g, r in results.items() if g != "chr1"}))

        for damage in (merged_tensor, merged_truncated, summary_record, summary_order, source_tensor, source_missing,
                       totals):
            with self.subTest(damage=damage.__name__):
                found, expected = self.merge_with(damage)
                self.assertIsInstance(found, str)
                self.assertEqual(found, expected)
        self.assertIn("Shard bytes differ", self.merge_with(merged_tensor)[0])
        self.assertIn("Spot check failed: chr2 chr2_shard_00000_data.npy[1]", self.merge_with(source_tensor)[0])


if __name__ == "__main__":
    unittest.main()

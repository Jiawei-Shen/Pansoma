import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "prepare_training_shards.py"


def write_shard(directory: Path, labels: list[int], offset: int = 0) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    data = np.zeros((len(labels), 5, 2, 3), dtype=np.float32)
    data[:, 0, 0, 0] = np.arange(offset, offset + len(labels))
    np.save(directory / "shard_00000_data.npy", data)
    np.save(directory / "shard_00000_labels.npy", np.asarray(labels, dtype=np.int8))


class PrepareTrainingShardsTest(unittest.TestCase):
    def test_random_split_filters_non_linear_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "tensors_chr1"
            output = root / "dataset"
            write_shard(source, [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1])

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input-dir",
                    str(source),
                    "--output-dir",
                    str(output),
                    "--val-fraction",
                    "0.2",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            train_labels = np.load(next((output / "train").glob("*_labels.npy")))
            val_labels = np.load(next((output / "val").glob("*_labels.npy")))
            self.assertEqual(len(train_labels), 8)
            self.assertEqual(len(val_labels), 2)
            self.assertEqual(set(np.unique(train_labels)), {0, 1})
            self.assertEqual(set(np.unique(val_labels)), {0, 1})
            summary = json.loads((output / "split_summary.json").read_text())
            self.assertEqual(summary["totals"]["ignored"], 2)

    def test_explicit_directories_remain_separate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            train_source = root / "tensors_chr2"
            val_source = root / "tensors_chr1"
            output = root / "dataset"
            write_shard(train_source, [0, 1, -1], offset=10)
            write_shard(val_source, [1, 0, -1], offset=20)

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--train-input-dir",
                    str(train_source),
                    "--val-input-dir",
                    str(val_source),
                    "--output-dir",
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            train_data = np.load(next((output / "train").glob("*_data.npy")))
            val_data = np.load(next((output / "val").glob("*_data.npy")))
            self.assertEqual(set(train_data[:, 0, 0, 0]), {10.0, 11.0})
            self.assertEqual(set(val_data[:, 0, 0, 0]), {20.0, 21.0})


if __name__ == "__main__":
    unittest.main()

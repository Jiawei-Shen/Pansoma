"""postprocess.merge_shards: per-chromosome shards, byte identity, provenance, bookkeeping, deletion."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ..tensor_postprocessing.chr_index import FIELDS, ChrIndex
from ..tensor_postprocessing.merge_shards import (LAYOUT, merge, verify_shard,
                                                              verify_summary)

SHAPE = [8, 4, 5]
BLOCKS = [dict(chrom="chr1", first_node=1, last_node=100, nodes=100, dataset="autosome", source="t"),
          dict(chrom="chr2", first_node=101, last_node=200, nodes=100, dataset="autosome", source="t"),
          dict(chrom="chrX", first_node=201, last_node=300, nodes=100, dataset="non_autosomal", source="t")]
# task -> nodes of its tensors (node-sorted, as the builder writes them); task 1 crosses chr1 -> chr2 -> chrX.
TASKS = {"SNV": {0: [3, 5, 5, 9, 60], 1: [99, 100, 101, 150, 201, 202, 250], 2: [260]},
         "INDEL": {0: [7], 1: [120, 130, 140, 141, 142, 143, 144, 145, 146], 2: []}}


def write_chr_index(path):
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(BLOCKS)
    return path


def fake_run(root, source_shard=3):
    """A finished split run: tensors/<kind>/task_NNNN with shards of `source_shard` tensors."""
    rng = np.random.default_rng(1)
    tensors, outputs = root / "tensors", {}
    expected = {}
    for kind, tasks in TASKS.items():
        outputs[kind] = []
        for task, nodes in tasks.items():
            folder = tensors / kind / f"task_{task:04d}"
            folder.mkdir(parents=True)
            data = rng.integers(-1, 127, size=(len(nodes), *SHAPE), dtype=np.int8)
            shards = [data[i:i + source_shard] for i in range(0, len(nodes), source_shard)]
            with (folder / "variant_summary.ndjson").open("w") as summary:
                for si, shard in enumerate(shards):
                    np.save(folder / f"shard_{si:05d}_data.npy", shard)
                    for i in range(len(shard)):
                        n = nodes[si * source_shard + i]
                        record = dict(candidate_id=f"{n}:{i}:SNP:A>C", node_id=n, start=i, af=0.1 * (i + 1),
                                      alleles=[dict(alt="C")], shard_index=si, index_within_shard=i)
                        summary.write(json.dumps(record) + "\n")
                        expected.setdefault(kind, []).append((n, shard[i].copy(), record))
            (folder / "filtered_candidates.ndjson").write_text(f'{{"task": {task}, "kind": "{kind}"}}\n')
            (folder / "unsupported_events.ndjson").write_text("")
            manifest = dict(status="complete", schema_version=5, tensor_format_version="indexed-gam-candidate-v5",
                            tensor_storage_version="int8-count-log2x14-v1", shape=SHAPE, dtype="int8",
                            channels=list("abcdefgh"), parameters=dict(min_af=0.06), sample_unit="site-v1",
                            tensors=len(nodes), shards=len(shards), graph_index=dict(path="g"))
            (folder / "manifest.json").write_text(json.dumps(manifest))
            outputs[kind].append(dict(task=task, path=str(folder), tensors=len(nodes), shards=len(shards)))
    for task in (0, 1, 2):
        shared = tensors / "shared" / f"task_{task:04d}"
        shared.mkdir(parents=True)
        (shared / "batch_timing.ndjson").write_text(json.dumps(dict(batch=1, elapsed_seconds=task)) + "\n")
    catalog = dict(root=str(root), tensors_dir=str(tensors), outputs=outputs,
                   tensors_by_type={k: sum(len(v) for v in t.values()) for k, t in TASKS.items()})
    catalog["tensors"] = sum(catalog["tensors_by_type"].values())
    (root / "outputs.json").write_text(json.dumps(catalog))
    (root / "status.json").write_text(json.dumps(dict(status="complete")))
    return tensors, expected


def group_of(node):
    return "chr1" if node <= 100 else "chr2" if node <= 200 else "chrX"


class MergeTest(unittest.TestCase):
    def run_merge(self, root, **kwargs):
        tensors, expected = fake_run(root)
        report = merge(root, write_chr_index(root / "chr.tsv"), shard_size=4, workers=2, spots=50, **kwargs)
        return tensors, expected, report

    def test_per_chromosome_shards_are_byte_identical_and_traceable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tensors, expected, report = self.run_merge(root, keep_sources=True)
            for kind, items in expected.items():
                for dataset, directory in (("autosome", tensors / kind), ("non_autosomal", tensors / "non_autosomal" / kind)):
                    manifest = json.loads((directory / "manifest.json").read_text())
                    self.assertEqual((manifest["layout"], manifest["shard_size"], manifest["dataset"]), (LAYOUT, 4, dataset))
                    for group, info in manifest["chromosomes"].items():
                        want = [e for e in items if group_of(e[0]) == group]
                        self.assertEqual(info["tensors"], len(want))
                        sizes = [s["tensors"] for s in info["shards"]]
                        self.assertEqual(sizes, [4] * (len(want) // 4) + ([len(want) % 4] if len(want) % 4 else []))
                        got = np.concatenate([np.load(directory / s["file"]) for s in info["shards"]])
                        np.testing.assert_array_equal(got, np.stack([e[1] for e in want]))  # source order kept
                        records = [json.loads(l) for l in (directory / info["summary"]).read_text().splitlines()]
                        for i, (r, (_, _, source)) in enumerate(zip(records, want)):
                            self.assertEqual((r["chrom"], r["shard_index"], r["index_within_shard"]), (group, i // 4, i % 4))
                            self.assertEqual(r["shard_file"], f"{group}_shard_{i // 4:05d}_data.npy")
                            self.assertEqual({k: v for k, v in r.items() if k not in (
                                "chrom", "shard_file", "shard_index", "index_within_shard", "source_task",
                                "source_shard_index", "source_index_within_shard")},
                                {k: v for k, v in source.items() if k not in ("shard_index", "index_within_shard")})
                            self.assertEqual((r["source_shard_index"], r["source_index_within_shard"]),
                                             (source["shard_index"], source["index_within_shard"]))
                            src = Path(tensors / kind / f"task_{r['source_task']:04d}")
                            np.testing.assert_array_equal(
                                np.load(directory / r["shard_file"])[r["index_within_shard"]],
                                np.load(src / f"shard_{r['source_shard_index']:05d}_data.npy")[r["source_index_within_shard"]])
                self.assertEqual({g for g in json.loads((tensors / "non_autosomal" / kind / "manifest.json").read_text())["chromosomes"]},
                                 {group_of(e[0]) for e in items} & {"chrX"})
                audit = (tensors / kind / "filtered_candidates.ndjson").read_text().splitlines()
                self.assertEqual([json.loads(l)["task"] for l in audit], [0, 1, 2])
                self.assertFalse((tensors / kind / ".merging").exists())
                self.assertTrue((tensors / kind / "task_0001" / "shard_00000_data.npy").exists())  # kept
            outputs = json.loads((root / "outputs.json").read_text())
            self.assertEqual(outputs["merge"]["kinds"]["SNV"]["autosome"]["chromosomes"], {"chr1": 7, "chr2": 2})
            self.assertEqual(outputs["merge"]["kinds"]["INDEL"]["autosome"]["chromosomes"], {"chr1": 1, "chr2": 9})
            self.assertTrue(outputs["merge"]["sources_kept"])
            self.assertTrue((root / "outputs.pre_merge.json").exists())
            status = json.loads((root / "status.json").read_text())
            self.assertEqual((status["status"], status["merged"]), ("complete", True))
            timing = [json.loads(l) for l in (root / "batch_timing.ndjson").read_text().splitlines()]
            self.assertEqual([t["task"] for t in timing], [0, 1, 2])
            with self.assertRaisesRegex(ValueError, "already merged"):
                merge(root, root / "chr.tsv", shard_size=4)

    def test_sources_deleted_only_when_asked_and_finalized(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tensors, _, _ = self.run_merge(root, keep_sources=False)
            self.assertFalse(list(tensors.glob("*/task_*")))
            self.assertTrue((tensors / "SNV" / "chr1_shard_00000_data.npy").exists())
            self.assertEqual(json.loads((root / "status.json").read_text())["status"], "finalized")

    def test_verification_detects_changed_bytes_and_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tensors, _, _ = self.run_merge(root, keep_sources=True)
            manifest = json.loads((tensors / "SNV" / "manifest.json").read_text())
            info = manifest["chromosomes"]["chr1"]
            shard = info["shards"][0]
            path = tensors / "SNV" / shard["file"]
            verify_shard(path, shard["tensors"], SHAPE, shard["sha256"])
            x = np.load(path)
            x[0, 0, 0, 0] ^= 1
            np.save(path, x)
            with self.assertRaisesRegex(ValueError, "differ from their source"):
                verify_shard(path, shard["tensors"], SHAPE, shard["sha256"])
            summary = tensors / "SNV" / info["summary"]
            verify_summary(tensors / "SNV", "chr1", dict(info, dataset="autosome"))
            lines = summary.read_text().splitlines()
            lines[1] = json.dumps(dict(json.loads(lines[1]), af=0.99))
            summary.write_text("\n".join(lines) + "\n")
            with self.assertRaisesRegex(ValueError, "differ from the source"):
                verify_summary(tensors / "SNV", "chr1", dict(info, dataset="autosome"))

    def test_nodes_outside_every_block_fail_before_publishing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            TASKS["SNV"][2].append(999)
            try:
                tensors, _ = fake_run(root)
                with self.assertRaisesRegex(ValueError, "outside every chromosome block"):
                    merge(root, write_chr_index(root / "chr.tsv"), shard_size=4, workers=1)
            finally:
                TASKS["SNV"][2].remove(999)
            self.assertFalse(list((tensors / "SNV").glob("chr*")))
            self.assertNotIn("merge", json.loads((root / "outputs.json").read_text()))

    def test_select_nodes(self):
        from ..tensor_postprocessing.chr_index import select_nodes
        with tempfile.TemporaryDirectory() as tmp:
            table = write_chr_index(Path(tmp) / "chr.tsv")
            nodes = [1, 50, 101, 150, 250, 999]
            kept, report = select_nodes(nodes, "all")
            self.assertEqual((kept.tolist(), report["nodes_kept"]), (nodes, 6))
            kept, report = select_nodes(nodes, "autosome", table)
            self.assertEqual((kept.tolist(), report["removed"]), ([1, 50, 101, 150], {"chrX": 1, "no_block": 1}))
            kept, report = select_nodes(nodes, "chr2,chrX", table)
            self.assertEqual((kept.tolist(), report["chromosomes"]), ([101, 150, 250], ["chr2", "chrX"]))
            with self.assertRaisesRegex(ValueError, "Unknown chromosome"):
                select_nodes(nodes, "chr3", table)
            with self.assertRaisesRegex(ValueError, "needs --chr-index"):
                select_nodes(nodes, "autosome")

    def test_chr_index_lookup(self):
        with tempfile.TemporaryDirectory() as tmp:
            index = ChrIndex(write_chr_index(Path(tmp) / "chr.tsv"))
            self.assertEqual(index.names_of([0, 1, 100, 101, 300, 301]), [None, "chr1", "chr1", "chr2", "chrX", None])


if __name__ == "__main__":
    unittest.main()

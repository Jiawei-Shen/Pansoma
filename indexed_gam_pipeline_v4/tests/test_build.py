"""End-to-end builds on synthetic data: the shared/SNV/INDEL output contract, limits and the auditor."""
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from .fixtures import af_gam, build_args, graph_fixture, simple_alignment, split_args, tiny_gam, write_gam
from ..gam_reader import record_key
from .. import build as build_module
from ..build import ARGUMENTS, PARAMETERS, build
from ..candidates import (BASES, CHANNELS, ROW_SELECTION, TENSOR_FORMAT, TENSOR_STORAGE, WINDOW_ENCODING,
                          decode_alignment, encode_count)
from ..orchestrate import validate_shards
from ..run import make_parser
from ..tools.validate_examples import validate


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

    def args(self, name, **overrides):
        """Outputs under root/name/{shared,SNV,INDEL}, AF 0.05 for both types."""
        return split_args(self.root, name, **dict(self.base, **overrides))

    def test_shards_summary_manifest_contract_and_audit(self):
        manifest = quiet_build(self.args("out", debug_rows=True))  # tiny_gam: four SNV sites, no indel
        out = self.root / "out"
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual((manifest["tensors"], manifest["shards"], manifest["shape"]), (4, 2, [8, 200, 101]))
        self.assertEqual((manifest["tensors_by_type"], manifest["output_layout"]), (dict(SNV=4, INDEL=0), "split"))
        self.assertEqual(manifest["variant_outputs"], {k: str(out / k) for k in ("SNV", "INDEL")})
        self.assertEqual(sorted(p.name for p in (out / "shared").iterdir()),
                         ["batch_timing.ndjson", "filtered_candidates.ndjson", "manifest.json",
                          "target_nodes.txt", "unsupported_events.ndjson"])  # no shards, no summary
        summary = [json.loads(s) for s in (out / "SNV/variant_summary.ndjson").read_text().splitlines()]
        self.assertEqual(len(summary), 4)
        for i, meta in enumerate(summary):
            self.assertEqual((meta["shard_index"], meta["index_within_shard"]), (i // 2, i % 2))
            x = np.load(out / f"SNV/shard_{i // 2:05d}_data.npy")[i % 2]
            self.assertEqual((list(x.shape), str(x.dtype)), (manifest["shape"], "int8"))
            self.assertEqual(int(x[6].max()), encode_count(331, 90))  # fixture default H 90
            self.assertEqual(meta["parameters"]["haplotypes"], 90)
            self.assertTrue(np.all(x[7][x[0] != 0] == 1) and not x[7][x[0] == 0].any())  # forward-only fixture
            self.assertEqual(meta["channels"], CHANNELS)
            self.assertEqual(meta["coverage"], meta["alt_count"] + meta["ref_count"] + meta["other_count"])
            self.assertEqual((meta["parameters"]["min_af"], meta["parameters"]["variant_type"]), (0.05, "snp"))
            self.assertEqual((meta["sample_unit"], meta["site_id"]), ("site", f"{meta['node_id']}:{meta['start']}:SNV"))
            self.assertEqual((meta["tensor_format"], meta["tensor_storage"], meta["row_selection"],
                              meta["window_encoding"]), (TENSOR_FORMAT, TENSOR_STORAGE, ROW_SELECTION, WINDOW_ENCODING))
            self.assertEqual([k for k in meta if "version" in k], [])
        typed = json.loads((out / "SNV/manifest.json").read_text())
        self.assertEqual((typed["status"], typed["tensors"], typed["shared_output"]), ("complete", 4, str(out / "shared")))
        # The format names are plain; gai_version is the GAI's own format number (vg's).
        self.assertEqual((TENSOR_FORMAT, TENSOR_STORAGE, ROW_SELECTION, WINDOW_ENCODING),
                         ("indexed-gam-candidate", "int8-count-haplotypes100", "site-allele-blocks-uniform-similarity",
                          "site-layout-columns"))
        self.assertEqual([k for k in typed if "version" in k], ["gai_version"])
        # parameters/arguments: the shared manifest's key tables; a typed one adds its kind (variant_type)
        # and records its own AF threshold and output.
        self.assertEqual((list(manifest["parameters"]), list(manifest["arguments"])), (list(PARAMETERS), list(ARGUMENTS)))
        for kind, variant_type in (("SNV", "snp"), ("INDEL", "indel")):
            m = json.loads((out / kind / "manifest.json").read_text())
            self.assertEqual(list(m["parameters"].items()),
                             list(dict(manifest["parameters"], min_af=0.05, variant_type=variant_type).items()))
            self.assertEqual(list(m["arguments"].items()),
                             list(dict(manifest["arguments"], output=str(out / kind), min_af=0.05).items()))
        self.assertEqual((json.loads((out / "INDEL/manifest.json").read_text())["tensors"],
                          (out / "INDEL/variant_summary.ndjson").read_text()), (0, ""))
        self.assertEqual(len((out / "shared/batch_timing.ndjson").read_text().splitlines()), 3)  # [10,20] [30] [1000]
        self.assertEqual(validate_shards(out / "SNV", 2)["tensors"], 4)
        audit = validate(out / "SNV", str(self.gam), None, str(self.graph))
        self.assertTrue(audit["passed"])
        self.assertEqual(audit["examples"], 4)
        # The auditor must catch a corrupted graph-reference channel.
        shard = out / "SNV/shard_00000_data.npy"
        corrupted = np.load(shard)
        corrupted[0, 5, 0, 49] = BASES["C"]
        np.save(shard, corrupted)
        with self.assertRaisesRegex(ValueError, "graph base"):
            validate(out / "SNV", str(self.gam), None, str(self.graph))

    def test_target_node_restriction_does_not_change_outputs(self):
        # The Python decoder, patched (tests/test_native_decoder.py covers the native one).
        options = dict(batch_nodes=1, max_node_span=10000, debug_rows=True, rows=4, decoder="python")
        unrestricted = lambda *a, **kw: decode_alignment(*a, **dict(kw, target_nodes=None))  # noqa: E731
        with patch.object(build_module, "decode_alignment", unrestricted):
            quiet_build(self.args("baseline", **options))
        quiet_build(self.args("restricted", **options))
        for kind in ("shared", "SNV", "INDEL"):
            self.assertEqual(files(self.root / "baseline" / kind), files(self.root / "restricted" / kind), kind)

    def test_haplotypes_set_channel_6_and_only_channel_6(self):
        """--haplotypes reaches every channel-6 cell (331 paths: 108 at H 90, 100 at H 331, 45 at H 464) and the manifest;
        nothing else in the outputs depends on it."""
        built = {}
        for h, code in ((90, 108), (331, 100), (464, 45)):
            manifest = quiet_build(self.args(f"h{h}", haplotypes=h))
            self.assertEqual(manifest["parameters"]["haplotypes"], h)
            self.assertIn(f"H = {h} haplotype paths", manifest["encodings"]["node_distinct_gbwt_path_count"])
            x = np.concatenate([np.load(f) for f in sorted((self.root / f"h{h}" / "SNV").glob("shard_*_data.npy"))])
            self.assertEqual(sorted(np.unique(x[:, 6]).tolist()), [0, code])
            built[h] = x
        for h in (331, 464):
            np.testing.assert_array_equal(np.delete(built[h], 6, axis=1), np.delete(built[90], 6, axis=1))
            np.testing.assert_array_equal(built[h][:, 6] != 0, built[90][:, 6] != 0)
        with self.assertRaisesRegex(ValueError, "--haplotypes"):
            quiet_build(self.args("h0", haplotypes=0))

    def test_missing_graph_node(self):
        graph_fixture(self.root / "partial.sqlite", [(n, "AAAAAA", 1) for n in (10, 20)])
        args = self.args("missing", graph_index=str(self.root / "partial.sqlite"))
        with self.assertRaisesRegex(ValueError, "missing from GBZ graph index"):
            quiet_build(args)

    def test_argument_validation(self):
        for overrides, message in ((dict(max_indel_len=51), "max-indel-len"),
                                   (dict(width=10, max_indel_len=20), "width"),
                                   (dict(gam_cache_mb=0), "gam-cache-mb must be positive"),
                                   (dict(snv_output=None), "Split output needs"),
                                   (dict(indel_min_af=None), "Split output needs"),
                                   (dict(indel_output=str(self.root / "o/SNV")), "distinct"),
                                   (dict(snv_output=str(self.root / "o/shared/SNV")), "non-nested")):
            with self.assertRaisesRegex(ValueError, message):
                quiet_build(self.args("o", **overrides))
        self.assertFalse((self.root / "o").exists())  # refused before any output directory exists
        # The CLI refuses a 0 MiB cache and a missing split option before any work.
        argv = ["build", "--gam", "g", "--output", "o", "--nodes", "n", "--graph-index", "x", "--snv-output", "s",
                "--indel-output", "i", "--haplotypes", "90", "--snv-min-af", ".1", "--indel-min-af", ".1"]
        self.assertEqual(make_parser().parse_args(argv).gam_cache_mb, 1024)
        for bad in (argv + ["--gam-cache-mb", "0"], argv[:-2], argv[:9] + argv[11:], argv[:13] + argv[15:],
                    argv + ["--haplotypes", "0"]):
            with self.assertRaises(SystemExit), redirect_stderr(io.StringIO()):
                make_parser().parse_args(bad)

    def test_every_stream_is_closed_before_status_complete(self):
        states, save = [], build_module.OutputDir.save

        def recording_save(output, **extra):
            save(output, **extra)
            if output.manifest["status"] == "complete":
                states.append((output.path.name, {k: s.closed for k, s in output.streams.items()}))

        with patch.object(build_module.OutputDir, "save", recording_save):
            quiet_build(self.args("closed"))
        self.assertEqual([name for name, _ in states], ["SNV", "INDEL", "shared"])
        for name, closed in states:
            self.assertTrue(all(closed.values()), (name, closed))
        self.assertEqual(set(states[0][1]), {"filtered", "unsupported", "summary"})


class SplitOutputTest(unittest.TestCase):
    def test_one_shared_pass_routes_types_with_independent_thresholds(self):
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
            reader_fetch = build_module.IndexedGam.fetch_capped

            def counting_fetch(self, *a, **kw):
                fetch_calls.append(1)
                return reader_fetch(self, *a, **kw)

            def counting_decode(*a, **kw):
                decode_calls.append(1)
                return decode_alignment(*a, **kw)

            with patch.object(build_module.IndexedGam, "fetch_capped", counting_fetch), \
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
            for kind in ("SNV", "INDEL"):
                self.assertEqual(validate_shards(root / kind, 1)["tensors"], len(snv) if kind == "SNV" else len(indel))
            # The INDEL threshold never touches an SNV byte (shards, summary and audit streams).
            again = build_args(**dict(common, output=str(root / "again"), snv_output=str(root / "again_SNV"),
                                      indel_output=str(root / "again_INDEL"), snv_min_af=.06, indel_min_af=.5))
            self.assertEqual(quiet_build(again)["tensors_by_type"], dict(SNV=1, INDEL=0))
            self.assertEqual(files(root / "SNV"), files(root / "again_SNV"))
            self.assertNotEqual(files(root / "INDEL"), files(root / "again_INDEL"))
            # DEL 3/50 = 0.06 fails the 0.08 INDEL threshold but would pass the SNV one.
            filtered = [json.loads(l) for l in (root / "INDEL/filtered_candidates.ndjson").read_text().splitlines()]
            self.assertIn((20, ["min_af"]), [(m["node_id"], m["reasons"]) for m in filtered])
            empty = build_args(**dict(common, output=str(root / "empty"), snv_output=str(root / "empty_SNV"),
                                      indel_output=str(root / "empty_INDEL"), snv_min_af=1, indel_min_af=1))
            self.assertEqual(quiet_build(empty)["tensors_by_type"], dict(SNV=0, INDEL=0))


if __name__ == "__main__":
    unittest.main()


class AdaptiveBatchTest(unittest.TestCase):
    """--batch-nodes auto changes only how targets are grouped: tensors and summaries stay identical."""

    def build_af(self, root, name, **options):
        inputs = root / (name + "_inputs")
        inputs.mkdir()
        gam, _ = af_gam(inputs)
        graph = graph_fixture(inputs / "graph.sqlite", [(n, "AAAAAA", 331) for n in (10, 20, 30, 40, 50)])
        nodes = inputs / "nodes.txt"
        nodes.write_text("10\n20\n30\n40\n50\n")
        args = build_args(gam=str(gam), nodes=str(nodes), graph_index=str(graph), shard_size=1, min_variants=3,
                          output=str(root / name / "shared"), snv_output=str(root / name / "SNV"),
                          indel_output=str(root / name / "INDEL"), snv_min_af=.06, indel_min_af=.08, **options)
        quiet_build(args)
        tensors = {k: files(root / name / k, names=("variant_summary.ndjson",)) for k in ("SNV", "INDEL")}
        timing = [json.loads(l) for l in (root / name / "shared/batch_timing.ndjson").read_text().splitlines()]
        return tensors, timing

    def test_auto_batches_keep_tensors_and_record_the_plan(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixed, _ = self.build_af(root, "fixed", batch_nodes=1)
            auto, timing = self.build_af(root, "auto", batch_nodes="auto")
            self.assertEqual(fixed, auto)
            self.assertEqual([t["target_nodes"] for t in timing], [5])
            self.assertEqual(timing[0]["batch_plan"]["size"], 512)
            self.assertEqual(sorted(timing[0]["batch_plan"]["bytes_per_node"]), ["1024", "2048", "512"])

    def test_auto_splits_a_batch_over_the_alignment_limit_instead_of_failing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixed, _ = self.build_af(root, "fixed", batch_nodes=1)
            auto, timing = self.build_af(root, "auto", batch_nodes="auto", max_batch_alignments=60)  # 50 per node
            self.assertEqual(fixed, auto)
            self.assertEqual([t["target_nodes"] for t in timing], [1, 1, 1, 1, 1])
            self.assertTrue(all(t["batch_plan"]["split"] for t in timing))
            with self.assertRaisesRegex(ValueError, "Batch alignment limit"):
                self.build_af(root, "strict", batch_nodes=5, max_batch_alignments=60)

    def test_adaptive_batches_prefer_large_batches_only_when_they_save_reads(self):
        class Reader:  # every node its own run of 100 blocks: no sharing, larger batches do not help
            def __init__(self, shared):
                self.shared = shared

            def ranges(self, batch):
                if self.shared:  # all nodes in one run: one large batch reads it once
                    return [(0, 100 << 16)]
                return [(n * 1000 << 16, (n * 1000 + 100) << 16) for n in batch]

        nodes = list(range(1, 5001))
        sizes = lambda shared: [len(b) for b, _ in build_module.adaptive_batches(nodes, Reader(shared), 10**9)]  # noqa: E731
        self.assertEqual(set(sizes(False)[:-1]), {512})
        self.assertEqual(set(sizes(True)[:-1]), {2048})
        self.assertEqual(sum(sizes(True)), 5000)


class CappedFetchTest(unittest.TestCase):
    """The read cap applied while reading the GAM: a node over --max-node-reads is built from its records with
    the smallest record digest only, exactly as if the GAM held nothing else for it, whatever the batching."""

    def build(self, root, name, rows, nodes, **options):
        inputs = root / (name + "_inputs")
        inputs.mkdir()
        gam = write_gam(inputs / "reads.gam", rows)
        graph = graph_fixture(inputs / "graph.sqlite", [(n, "AAAAAA", 331) for n in nodes])
        (inputs / "nodes.txt").write_text("".join(f"{n}\n" for n in nodes))
        options = dict(dict(batch_nodes=2, min_variants=1, snv_min_af=.01, indel_min_af=.01), **options)
        quiet_build(build_args(gam=str(gam), nodes=str(inputs / "nodes.txt"), graph_index=str(graph), shard_size=1,
                               output=str(root / name / "shared"), snv_output=str(root / name / "SNV"),
                               indel_output=str(root / name / "INDEL"), **options))
        return root / name

    @staticmethod
    def sites(folder, drop=("shard_index", "index_within_shard")):
        """{site_id: (summary without `drop` and the recorded max_node_reads, tensor bytes)} of both typed outputs."""
        found = {}
        for kind in ("SNV", "INDEL"):
            for line in (folder / kind / "variant_summary.ndjson").read_text().splitlines():
                meta = json.loads(line)
                meta["parameters"].pop("max_node_reads")
                shard = np.load(folder / kind / f"shard_{meta['shard_index']:05d}_data.npy")
                found[meta["site_id"]] = ({k: v for k, v in meta.items() if k not in drop},
                                          shard[meta["index_within_shard"]].tobytes())
        return found

    @staticmethod
    def af_rows():
        with tempfile.TemporaryDirectory() as tmp:
            return af_gam(Path(tmp))[1]  # five nodes x 50 records

    def test_nodes_within_the_cap_are_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, rows, nodes = Path(tmp), self.af_rows(), (10, 20, 30, 40, 50)
            capped = self.build(root, "capped", rows, nodes, max_node_reads=50)
            uncapped = self.build(root, "uncapped", rows, nodes, max_node_reads=0)
            self.assertEqual(self.sites(capped), self.sites(uncapped))
            for f in ("SNV/filtered_candidates.ndjson", "INDEL/filtered_candidates.ndjson"):
                self.assertEqual((capped / f).read_text(), (uncapped / f).read_text())
            self.assertFalse((capped / "shared/downsampled_nodes.tsv").exists())
            self.assertNotIn("downsampled_nodes", json.loads((capped / "shared/manifest.json").read_text()))

    def test_a_node_over_the_cap_is_built_from_its_smallest_digest_records_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, rows, nodes = Path(tmp), self.af_rows(), (10, 20, 30, 40, 50)
            kept = []
            for node in nodes:  # each record visits one node: the GAM a 40-record cap leaves for it
                mine = [a for a in rows if a.path.mapping[0].position.node_id == node]
                keys = sorted(record_key(a) for a in mine)[:40]
                kept += [a for a in mine if record_key(a) in set(keys)]
            capped = self.build(root, "capped", rows, nodes, max_node_reads=40)
            alone = self.build(root, "alone", kept, nodes, max_node_reads=0)
            sites = self.sites(capped, drop=("shard_index", "index_within_shard", "downsampled_from"))
            self.assertTrue(sites)
            self.assertEqual(sites, self.sites(alone))
            for f in ("SNV/filtered_candidates.ndjson", "INDEL/filtered_candidates.ndjson"):
                self.assertEqual((capped / f).read_text(), (alone / f).read_text())
            self.assertTrue(all(meta["downsampled_from"] == 50 for meta, _ in self.sites(capped).values()))
            self.assertEqual((capped / "shared/downsampled_nodes.tsv").read_text(), "node\trecords\tkept\treason\n" +
                             "".join(f"{n}\t50\t40\tmax_node_reads\n" for n in nodes))
            self.assertEqual(json.loads((capped / "shared/manifest.json").read_text())["downsampled_nodes"], 5)

    def test_capped_records_do_not_depend_on_the_batching(self):
        rows = [simple_alignment(nodes, mutation=i % 7 < 2, mapq=5 if i % 11 == 0 else 60, name=f"w{i}")
                for i, nodes in enumerate([(11, 12, 13), (12,), (11, 12), (12, 13, 14), (10, 11)][i % 5] for i in range(90))]
        with tempfile.TemporaryDirectory() as tmp:
            root, nodes = Path(tmp), (10, 11, 12, 13, 14)
            one = self.build(root, "one", rows, nodes, max_node_reads=12, batch_nodes=1)
            for name, batch in (("two", 2), ("five", 5), ("auto", "auto")):
                other = self.build(root, name, rows, nodes, max_node_reads=12, batch_nodes=batch)
                self.assertEqual(self.sites(other), self.sites(one), name)
                self.assertEqual((other / "shared/downsampled_nodes.tsv").read_text(),
                                 (one / "shared/downsampled_nodes.tsv").read_text())

    def test_a_collapsed_repeat_decodes_one_capped_set(self):
        """Nodes sharing all their reads keep the same capped records: a batch decodes the cap, not the depth."""
        rows = [simple_alignment((11, 12, 13, 14), mutation=i % 5 == 0, name=f"c{i}") for i in range(300)]
        with tempfile.TemporaryDirectory() as tmp:
            out = self.build(Path(tmp), "repeat", rows, (11, 12, 13, 14), max_node_reads=20, batch_nodes=4)
            timing = [json.loads(l) for l in (out / "shared/batch_timing.ndjson").read_text().splitlines()]
            self.assertEqual([(t["target_nodes"], t["alignments"], t["gam_query"]["sampled_nodes"]) for t in timing],
                             [(4, 20, 4)])

    def test_the_batch_limit_counts_capped_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, rows, nodes = Path(tmp), self.af_rows(), (10, 20, 30, 40, 50)
            single = self.build(root, "single", rows, nodes, max_node_reads=40, batch_nodes=1)
            auto = self.build(root, "auto", rows, nodes, max_node_reads=40, batch_nodes="auto", max_batch_alignments=45)
            self.assertEqual(self.sites(auto), self.sites(single))  # 5 x 40 records split down to single nodes
            with self.assertRaisesRegex(ValueError, "Batch alignment limit"):  # a fixed multi-node batch fails
                self.build(root, "fixed", rows, nodes, max_node_reads=40, batch_nodes=2, max_batch_alignments=45)
            with self.assertRaisesRegex(ValueError, "Batch alignment limit"):  # no cap: one node over the limit
                self.build(root, "uncapped", rows, nodes, max_node_reads=0, batch_nodes=1, max_batch_alignments=45)

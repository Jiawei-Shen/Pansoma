"""Site units, the AF upper-bound prefilter and the per-node read cap."""
import contextlib
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from .fixtures import (SEQ, graph_fixture, mixed_af_rows, overlap, site_rows, spec_alignment,
                       split_args, write_gam)
from ..build import ALLELE_FIELDS, ARGUMENTS, KIND, PARAMETERS, build, capped_reads, candidate_units, recorded
from ..candidates import BASES, Candidate, decode_alignment, exact_coverage
from ..orchestrate import BUILDER_OPTIONS, PACKAGE, build_command, validate_shards
from ..orchestrate import make_parser as prepare_parser
from ..run import make_parser
from ..tools.validate_examples import validate


def run(root, name, make_rows, **overrides):
    """A build of `make_rows` over nodes 10 and 20 into root/name/{shared,SNV,INDEL}."""
    gam = root / f"{make_rows.__name__}.gam"
    if not gam.exists():
        write_gam(gam, make_rows())
        graph_fixture(root / "graph.sqlite", [(n, s, 7) for n, s in SEQ.items()])
    (root / "nodes.txt").write_text("10\n20\n")
    options = dict(gam=str(gam), nodes=str(root / "nodes.txt"), graph_index=str(root / "graph.sqlite"),
                   batch_nodes=512, shard_size=2048)
    options.update(overrides)
    with redirect_stdout(io.StringIO()):
        build(split_args(root, name, **options))
    return root / name


def records(folder):
    metas = [json.loads(line) for line in (folder / "variant_summary.ndjson").read_text().splitlines()]
    shards = sorted(folder.glob("shard_*_data.npy"))
    return metas, (np.concatenate([np.load(p) for p in shards]) if shards else np.zeros((0, 8, 200, 101)))


def filtered(folder):
    return [json.loads(line) for line in (folder / "filtered_candidates.ndjson").read_text().splitlines()]


def without(meta, *keys):
    return {k: v for k, v in meta.items() if k not in keys}


class ExactCoverageTest(unittest.TestCase):
    def test_matches_overlap_for_every_position_event_strand_and_repeat(self):
        seq = {1: "ACGT", 2: "TT"}
        specs = [
            ([(1, 0, False, [(1, 1, ""), (1, 1, "T"), (2, 2, "")]), (2, 0, False, [(2, 2, "")]),
              (1, 0, False, [(1, 1, ""), (1, 1, "T"), (2, 2, "")])], "repeat"),
            ([(1, 0, True, [(2, 2, ""), (1, 1, "A"), (1, 1, "")])], "reverse"),
            ([(1, 0, False, [(1, 1, ""), (2, 0, ""), (1, 1, "")])], "del"),
            ([(1, 0, False, [(2, 2, ""), (0, 2, "TT"), (2, 2, "")])], "ins"),
            ([(1, 1, False, [(2, 2, "")])], "partial"),
            ([(1, 3, False, [(1, 1, "")]), (2, 0, False, [(2, 2, "")])], "edge"),
        ]
        reads = [decode_alignment(spec_alignment(s, seq, name=n), seq)[0] for s, n in specs]
        by_node = {n: [r for r in reads if any(v.node == n for v in r.visits)] for n in seq}
        probes = set()
        for node, s in seq.items():
            for pos in range(len(s) + 1):
                probes.add(Candidate(node, pos, "", "G", "INS"))
                for length in (1, 2, 3):
                    if pos + length <= len(s):
                        probes.add(Candidate(node, pos, s[pos:pos + length], "", "DEL"))
                if pos < len(s):
                    probes.add(Candidate(node, pos, s[pos], "C" if s[pos] != "C" else "A", "SNP"))
        coverage = exact_coverage(probes, reads + reads[:2])  # duplicate records count separately
        for c in probes:
            expected = sum(overlap(r, c, 10) is not None for r in by_node.get(c.node, []) + reads[:2]
                           if any(v.node == c.node for v in r.visits))
            self.assertEqual(coverage[c], expected, c)
        self.assertEqual(exact_coverage([Candidate(99, 0, "A", "C", "SNP")], reads), {Candidate(99, 0, "A", "C", "SNP"): 0})


class ReadCapTest(unittest.TestCase):
    def test_cap_is_deterministic_order_free_and_optional(self):
        seq = {1: "ACGT"}
        reads = [decode_alignment(spec_alignment([(1, 0, False, [(1, 1, ""), (1, 1, "T"), (2, 2, "")])], seq,
                                                 name=f"r{i}"), seq)[0] for i in range(6)]
        smallest = sorted(r.digest for r in reads)[:3]
        self.assertEqual([r.digest for r in capped_reads(reads, 3)], smallest)
        self.assertEqual([r.digest for r in capped_reads(list(reversed(reads)), 3)], smallest)
        self.assertIs(capped_reads(reads, 0), reads)
        self.assertIs(capped_reads(reads, 6), reads)  # not deeper than the cap: untouched

    def test_capped_build_passes_cap_aware_audit(self):
        seq = {1: "ACGTACGTAC"}
        rows = [spec_alignment([(1, 0, False, [(4, 4, ""), (1, 1, "C" if i % 2 else ""), (5, 5, "")])], seq,
                               name=f"deep{i}") for i in range(12)]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gam = write_gam(root / "deep.gam", rows)
            graph = graph_fixture(root / "graph.sqlite", [(1, seq[1], 3)])
            (root / "nodes.txt").write_text("1\n")
            common = dict(gam=str(gam), nodes=str(root / "nodes.txt"), graph_index=str(graph), debug_rows=True, rows=4)
            for cap, coverage in ((5, 5), (0, 12)):
                out = root / f"cap{cap}" / "SNV"
                with redirect_stdout(io.StringIO()):
                    build(split_args(root, f"cap{cap}", **common, max_node_reads=cap))
                (meta,), _ = records(out)
                self.assertEqual(meta["coverage"], coverage)
                self.assertEqual(meta["parameters"]["max_node_reads"], cap)
                report = validate(out, str(gam), None, str(graph))
                self.assertTrue(report["passed"])
                self.assertEqual(report["capped_nodes"], int(cap > 0))
            # The auditor re-applies the recorded cap: a manifest that hides it must fail.
            manifest = json.loads((root / "cap5/SNV/manifest.json").read_text())
            manifest["parameters"]["max_node_reads"] = 0
            (root / "cap5/SNV/manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "independent (representative|site) coverage"):
                validate(root / "cap5/SNV", str(gam), None, str(graph))


class EarlyAfFilterTest(unittest.TestCase):
    def test_without_cap_outputs_are_unchanged_and_rejections_equal(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            common = dict(max_node_reads=0, snv_min_af=.3, indel_min_af=.3, min_variants=1, batch_nodes=2,
                          shard_size=2, rows=4, debug_rows=True)
            base = run(root, "base", mixed_af_rows, early_af_filter=False, **common)
            fast = run(root, "fast", mixed_af_rows, early_af_filter=True, **common)
            manifest = json.loads((fast / "shared/manifest.json").read_text())
            self.assertEqual(manifest["early_af_rejected"], 4)  # G>T, DEL, INS on node 10; G>C on node 20
            self.assertEqual(manifest["early_rejected"], 0)
            strip = lambda m: dict(m, parameters=without(m["parameters"], "early_af_filter"))  # noqa: E731
            for kind in ("SNV", "INDEL"):
                shards = sorted(p.name for p in (base / kind).glob("shard_*_data.npy"))
                self.assertEqual(shards, sorted(p.name for p in (fast / kind).glob("shard_*_data.npy")))
                for name in shards:
                    self.assertEqual((base / kind / name).read_bytes(), (fast / kind / name).read_bytes())
                self.assertEqual([strip(m) for m in records(base / kind)[0]], [strip(m) for m in records(fast / kind)[0]])
            # A>C 5/8 on node 10 and C>A 2/5 on node 20 pass 0.3; every indel is rejected early.
            self.assertEqual([m["site_id"] for m in records(fast / "SNV")[0]], ["10:4:SNV", "20:3:SNV"])
            self.assertEqual(records(fast / "INDEL")[0], [])
            ids = lambda f: sorted(m["candidate_id"] for m in filtered(f / "shared"))  # noqa: E731
            self.assertEqual(ids(base), ids(fast))
            early = [m for m in filtered(fast / "shared") if m.get("support_not_evaluated")]
            self.assertEqual(len(early), 4)
            for m in early:
                self.assertEqual(m["reasons"], ["min_af"])
                self.assertLess(m["af_upper_bound"], .3)
                self.assertEqual(m["af_upper_bound"], m["alt_support_upper_bound"] / m["coverage"])
                self.assertIn("site_id", m)


class SiteTest(unittest.TestCase):
    def test_site_tensor_holds_every_passing_allele(self):
        # (the bytes of this build are frozen as golden G3.)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            site = run(root, "site", site_rows, snv_min_af=.2, indel_min_af=.1, min_variants=2)
            observed = {o.candidate.metadata()["candidate_id"]: o.candidate.kind for r in site_rows()
                        for o in decode_alignment(r, SEQ)[0].observations}  # every allele on the target nodes
            for kind, expected in (("SNV", {"10:4:SNV": ["A>C", "A>G"], "20:3:SNV": ["C>A"]}),
                                   ("INDEL", {"10:7:INDEL": ["DEL", "INS"]})):
                smeta, stensor = records(site / kind)
                got = {m["site_id"]: [a["ref"] + ">" + a["alt"] if kind == "SNV" else a["event_type"]
                                      for a in m["alleles"]] for m in smeta}
                self.assertEqual(got, expected)
                for i, m in enumerate(smeta):
                    self.assertEqual((m["sample_unit"], m["allele_count"]), ("site-v2", len(m["alleles"])))
                    self.assertEqual(list(m["allele_labels"]), [a["label"] for a in m["alleles"]])
                    for a in m["alleles"]:
                        self.assertEqual(list(a), list(ALLELE_FIELDS) + ["label"])
                        self.assertEqual(a["coverage"], a["alt_count"] + a["ref_count"] + a["other_count"])
                        self.assertEqual(a["af"], a["alt_count"] / a["coverage"])
                    # top-level identity and counts are A1's (the representative)
                    for k in ("candidate_id", "coverage", "alt_count", "ref_count", "other_count", "af"):
                        self.assertEqual(m[k], m["alleles"][0][k])
                    self.assertEqual(m["site_counts"]["A1"], m["alleles"][0]["alt_count"])
                    if len(m["alleles"]) > 1:  # every passing allele has its own row block
                        self.assertGreater(m["site_counts"]["A2"], 0)
                        self.assertGreaterEqual(m["site_coverage"], m["coverage"])
                        for g in m["row_groups"]:  # each block spells its own allele in channel 2
                            label = g["allele"]
                            if label.startswith("A"):
                                carried = m["alleles"][int(label[1:]) - 1]
                                rows = stensor[i][:, g["start_row"]:g["end_row"]]
                                lo, hi = m["candidate_columns"]
                                cells = rows[0, :, lo:hi] != 0
                                if carried["event_type"] == "SNP":
                                    self.assertTrue(np.all(rows[2, :, lo][cells[:, 0]] == BASES[carried["alt"]]))
                                self.assertTrue(np.all(rows[2, :, lo:hi][cells] == rows[0, :, lo:hi][cells]))
                    self.assertEqual(m["second_allele_af"], m["alleles"][1]["af"] if len(m["alleles"]) > 1 else 0.0)
                # Every allele appears exactly once: in some site's alleles[] or in the filtered log.
                site_ids = [a["candidate_id"] for m in smeta for a in m["alleles"]] + \
                           [m["candidate_id"] for m in filtered(site / kind)]
                self.assertEqual(len(site_ids), len(set(site_ids)))
                self.assertEqual(sorted(site_ids), sorted(c for c, k in observed.items() if KIND[k] == kind))
                self.assertTrue(all("site_id" in m for m in filtered(site / kind)))
                self.assertEqual(validate_shards(site / kind, 2048)["sites"], len(expected))
            snv = {m["site_id"]: m for m in records(site / "SNV")[0]}
            self.assertEqual([(a["alt"], a["alt_count"]) for a in snv["10:4:SNV"]["alleles"]], [("C", 4), ("G", 3)])
            self.assertEqual(json.loads((site / "SNV/manifest.json").read_text())["sample_unit"], "site-v2")
            # A>T (1 read) fails min_variants in the prefilter and carries its site id.
            (a_t,) = [m for m in filtered(site / "SNV") if m["candidate_id"] == "10:4:SNP:A>T"]
            self.assertEqual((a_t["reasons"], a_t["site_id"]), (["min_variants"], "10:4:SNV"))

    def test_multiallelic_debug_build_passes_the_independent_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = run(root, "debug", site_rows, snv_min_af=.2, indel_min_af=.1, min_variants=2, debug_rows=True, rows=7)
            gam = root / "site_rows.gam"
            for kind in ("SNV", "INDEL"):
                report = validate(out / kind, str(gam), None, str(root / "graph.sqlite"))
                self.assertTrue(report["passed"])
                self.assertTrue(any(r["alleles"] > 1 for r in report["candidates"]), kind)
            # A tampered site-allele cell is caught.
            (meta,) = [m for m in records(out / "INDEL")[0] if m["allele_count"] > 1]
            shard = out / "INDEL" / f"shard_{meta['shard_index']:05d}_data.npy"
            x = np.load(shard)
            row = next(r for r in range(meta["selected_alignments"])
                       if x[meta["index_within_shard"], 0, r, meta["anchor_column"]])  # a cell with evidence
            x[meta["index_within_shard"], 2, row, meta["anchor_column"]] += 1
            np.save(shard, x)
            with self.assertRaisesRegex(ValueError, "site allele channel"):
                validate(out / "INDEL", str(gam), None, str(root / "graph.sqlite"))

    def test_mixed_indel_site_layout_keeps_flanks_aligned(self):
        seq = {1: "GCATGACTGA"}
        ins = [(4, 4, ""), (0, 3, "CCC"), (6, 6, "")]
        short = [(4, 4, ""), (0, 1, "C"), (6, 6, "")]
        long = [(4, 4, ""), (0, 5, "CCCCC"), (6, 6, "")]
        dele = [(4, 4, ""), (2, 0, ""), (4, 4, "")]
        ref = [(10, 10, "")]
        reads = [decode_alignment(spec_alignment([(1, 0, False, e)], seq, name=n), seq)[0]
                 for n, e in (("ins", ins), ("short", short), ("long", long), ("del", dele), ("ref", ref))]
        from ..candidates import NodeReads, SiteLayout, make_site_tensor
        a1, a2 = Candidate(1, 4, "", "CCC", "INS"), Candidate(1, 4, "GA", "", "DEL")
        self.assertEqual(SiteLayout.of([a1, a2]), SiteLayout(1, 4, 3, 2, "GA"))
        node = NodeReads(1, reads, 15)
        x, m = make_site_tensor([a1, a2], [node.eligible(a1, 10), node.eligible(a2, 10)], {1: 5}, rows=6, width=15,
                                debug=True)
        rows = {r["read_name"]: (i, r) for i, r in enumerate(m["rows"])}
        self.assertEqual({n: r["site_allele"] for n, (_, r) in rows.items()},
                         dict(ins="A1", short="OTHER", long="OTHER", ref="REF", **{"del": "A2"}))
        lo, hi = m["candidate_columns"]
        self.assertEqual((lo, hi), (7, 12))  # 3 insertion slots + 2 spanned bases
        for name, (ri, r) in rows.items():
            channel = "".join(".ACGTN-"[v] for v in x[2, ri, lo:hi])
            expected = {"A1": "CCCGA", "A2": "-----", "REF": "---GA", "OTHER": "....."}[r["site_allele"]]
            self.assertEqual(channel, expected, name)
            # the right flank (graph bases after the span) sits at the same columns in every row
            self.assertEqual("".join(".ACGTN-"[v] for v in x[5, ri, hi:hi + 3]), "CTG", name)
        self.assertEqual(rows["long"][1]["cropped_inserted_bases"], 2)
        self.assertEqual(x[0, rows["short"][0], lo:lo + 3].tolist(), [BASES["C"], 6, 6])

    def test_units_group_ins_and_del_but_not_snv(self):
        c = [Candidate(1, 5, "", "G", "INS"), Candidate(1, 5, "A", "", "DEL"), Candidate(1, 5, "A", "C", "SNP"),
             Candidate(1, 5, "A", "T", "SNP"), Candidate(1, 6, "C", "G", "SNP")]
        self.assertEqual(candidate_units(sorted(c)),  # sites (1,5,INDEL) < (1,5,SNV) < (1,6,SNV);
                         [(c[0], c[1]), (c[2], c[3]), (c[4],)])  # alleles keep candidate order

    def test_validation_rejects_duplicate_sites_and_bad_allele_lists(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = run(Path(tmp), "site", site_rows, snv_min_af=.2, indel_min_af=.1, min_variants=2) / "SNV"
            original = (folder / "variant_summary.ndjson").read_text()
            first, second = (json.loads(line) for line in original.splitlines())
            for tamper, message in ((dict(second, site_id=first["site_id"], node_id=first["node_id"], start=first["start"],
                                          alleles=first["alleles"], allele_count=first["allele_count"],
                                          allele_labels=first["allele_labels"],
                                          candidate_id=first["candidate_id"]), "Duplicate site"),
                                    (dict(second, allele_count=2), "allele list mismatch")):
                (folder / "variant_summary.ndjson").write_text(json.dumps(first) + "\n" + json.dumps(tamper) + "\n")
                with self.assertRaisesRegex(ValueError, message):
                    validate_shards(folder, 2048)


class OrchestratorOptionsTest(unittest.TestCase):
    # a production-shaped config.builder; min_allele_bq is the int 10 (argparse default of a float option)
    E2E_BUILDER = dict(rows=200, width=101, gam_cache_mb=6144, batch_nodes=512, max_node_span=10000,
                       max_batch_alignments=20000, shard_size=2048, min_mapq=10, min_af=0.05, min_variants=3,
                       min_allele_bq=10, max_indel_len=50, chromosomes="autosome", max_node_reads=800,
                       early_af_filter=True, decoder="auto")
    MANIFEST_ARGUMENTS = ["command", "gam", "output", "nodes", "index", "graph_index", "snv_min_af", "indel_min_af",
                    "snv_output", "indel_output", "debug_rows", "max_tensors", "variant_type", "rows", "width",
                    "gam_cache_mb", "batch_nodes", "max_node_span", "max_batch_alignments", "shard_size", "min_mapq",
                    "min_af", "min_variants", "min_allele_bq", "max_indel_len", "candidate_unit", "max_node_reads",
                    "chromosomes", "chr_index", "early_af_filter", "decoder"]

    def config(self, **builder):
        return dict(python="python", tensors="/t", variant_outputs=dict(SNV=0.06, INDEL=0.08),
                    builder=dict(self.E2E_BUILDER, **builder), parts=[dict(nodes_file="/n")],
                    inputs=dict(gam=dict(path="/g"), index=dict(path="/g.gai"), graph_index=dict(path="/x"),
                                chr_index=dict(path="/c.tsv")))

    def test_frozen_command_carries_cap_prefilter_and_decoder_explicitly(self):
        joined = " ".join(build_command(Path("/root"), self.config(), 0))
        self.assertIn("--chromosomes autosome --max-node-reads 800 --early-af-filter --decoder auto", joined)
        self.assertNotIn("--variant-type", joined)
        self.assertNotIn("--candidate-unit", joined)
        joined = " ".join(build_command(Path("/root"), self.config(max_node_reads=0, early_af_filter=False,
                                                                  decoder="python"), 0))
        self.assertIn("--max-node-reads 0 --no-early-af-filter --decoder python", joined)

    def test_every_builder_option_is_a_dest_of_the_run_parser(self):
        args = make_parser().parse_args(build_command(Path("/root"), self.config(), 0)[3:])
        self.assertEqual({k: getattr(args, k) for k in BUILDER_OPTIONS}, dict(self.E2E_BUILDER, min_allele_bq=10.0))
        self.assertNotIn("candidate_unit", BUILDER_OPTIONS)

    def test_every_builder_option_is_a_dest_of_the_prepare_parser(self):
        """prepare freezes {k: getattr(args, k) for k in BUILDER_OPTIONS}, with the run parser's
        defaults except --gam-cache-mb 8192 (min_allele_bq stays the int 10 when not given)."""
        base = ["prepare", "--root", "r", "--gam", "g", "--nodes", "n", "--graph-index", "x"]
        args = prepare_parser().parse_args(base + ["--snv-min-af", ".06", "--indel-min-af", ".08"])
        run_args = make_parser().parse_args(["build", "--gam", "g", "--output", "o", "--nodes", "n", "--graph-index", "x",
                                             "--snv-output", "s", "--indel-output", "i", "--snv-min-af", ".06",
                                             "--indel-min-af", ".08"])
        frozen = {k: getattr(args, k) for k in BUILDER_OPTIONS}
        self.assertEqual(frozen, dict({k: getattr(run_args, k) for k in BUILDER_OPTIONS}, gam_cache_mb=8192))
        self.assertEqual(frozen["batch_nodes"], "auto")
        self.assertEqual((type(frozen["min_allele_bq"]), frozen["decoder"]), (int, "auto"))
        # the typed thresholds are required, so config.variant_outputs is always {SNV, INDEL}
        for missing in (["--snv-min-af", ".06"], ["--indel-min-af", ".08"]):
            with redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaises(SystemExit):
                prepare_parser().parse_args(base + missing)

    def test_parameters_and_arguments_bytes_are_the_manifest_format(self):
        """config.builder -> build_command (str round trip) -> run parser -> recorded(): the manifest's exact bytes
        (the key tables every HG008 dataset was written with)."""
        command = build_command(Path("/root"), self.config(), 0)
        self.assertEqual(command[1:4], ["-m", f"{PACKAGE}.run", "build"])
        args = make_parser().parse_args(command[3:])
        shared = recorded(args, PARAMETERS)
        self.assertEqual(json.dumps(shared), '{"min_mapq": 10, "min_af": 0.05, "min_variants": 3, "min_allele_bq": 10.0, '
                         '"max_indel_len": 50, "variant_type": "all", "rows": 200, "width": 101, "max_node_reads": 800, '
                         '"candidate_unit": "site", "early_af_filter": true}')
        # build.build derives the typed parameters with dict(parameters, variant_type=..., min_af=...)
        self.assertEqual(json.dumps(dict(shared, variant_type="snp", min_af=args.snv_min_af)),
                         '{"min_mapq": 10, "min_af": 0.06, "min_variants": 3, "min_allele_bq": 10.0, "max_indel_len": 50, '
                         '"variant_type": "snp", "rows": 200, "width": 101, "max_node_reads": 800, '
                         '"candidate_unit": "site", "early_af_filter": true}')
        self.assertEqual(json.dumps(dict(shared, variant_type="indel", min_af=args.indel_min_af)),
                         '{"min_mapq": 10, "min_af": 0.08, "min_variants": 3, "min_allele_bq": 10.0, "max_indel_len": 50, '
                         '"variant_type": "indel", "rows": 200, "width": 101, "max_node_reads": 800, '
                         '"candidate_unit": "site", "early_af_filter": true}')
        arguments = recorded(args, ARGUMENTS)
        self.assertEqual(list(ARGUMENTS), self.MANIFEST_ARGUMENTS)
        self.assertEqual(json.dumps(arguments), '{"command": "build", "gam": "/g", "output": "/t/shared/task_0000", '
                         '"nodes": "/n", "index": "/g.gai", "graph_index": "/x", "snv_min_af": 0.06, "indel_min_af": 0.08, '
                         '"snv_output": "/t/SNV/task_0000", "indel_output": "/t/INDEL/task_0000", "debug_rows": false, '
                         '"max_tensors": null, "variant_type": "all", "rows": 200, "width": 101, "gam_cache_mb": 6144, '
                         '"batch_nodes": 512, "max_node_span": 10000, "max_batch_alignments": 20000, "shard_size": 2048, '
                         '"min_mapq": 10, "min_af": 0.05, "min_variants": 3, "min_allele_bq": 10.0, "max_indel_len": 50, '
                         '"candidate_unit": "site", "max_node_reads": 800, "chromosomes": "autosome", '
                         '"chr_index": "/c.tsv", "early_af_filter": true, "decoder": "auto"}')
        # A standalone build without --min-allele-bq records the int default.
        args = make_parser().parse_args(["build", "--gam", "g", "--output", "o", "--nodes", "n", "--graph-index", "x",
                                         "--snv-output", "s", "--indel-output", "i", "--snv-min-af", ".1",
                                         "--indel-min-af", ".1"])
        self.assertIn('"min_allele_bq": 10,', json.dumps(recorded(args, PARAMETERS)))
        self.assertIn('"min_allele_bq": 10,', json.dumps(recorded(args, ARGUMENTS)))


if __name__ == "__main__":
    unittest.main()

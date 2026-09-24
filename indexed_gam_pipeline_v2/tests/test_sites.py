"""Site units, the AF upper-bound prefilter and the per-node read cap."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from fixtures import build_args, graph_fixture, spec_alignment, write_gam  # noqa: F401  (sets sys.path)
from indexed_gam_pipeline_v2.build import build, capped_reads, candidate_units
from indexed_gam_pipeline_v2.candidates import BASES, Candidate, decode_alignment, exact_coverage, overlap
from indexed_gam_pipeline_v2.orchestrate import build_command, validate_shards
from indexed_gam_pipeline_v2.validate_examples import validate

SEQ = {10: "ACGTACGTAC", 20: "TTGCAAGGCT"}


def site_rows():
    """Node 10: position 4 A>C x4, A>G x3, A>T x1; 1-bp DEL at 7 x3 and INS TT at 7 x2 (11 reads;
    neither indel can left-shift: the base before position 7 is G). Node 20: position 3 C>A x3 (5 reads)."""
    snv = ["C"] * 4 + ["G"] * 3 + ["T"] + [None] * 3
    indel = {0: "D", 4: "D", 8: "D", 1: "I", 9: "I"}
    rows = []
    for i, base in enumerate(snv):
        edits = [(4, 4, ""), (1, 1, base or "")]
        edits += {"D": [(2, 2, ""), (1, 0, ""), (2, 2, "")],
                  "I": [(2, 2, ""), (0, 2, "TT"), (3, 3, "")]}.get(indel.get(i), [(5, 5, "")])
        rows.append(spec_alignment([(10, 0, False, edits)], SEQ, name=f"n10_{i}"))
    for i in range(5):
        edits = [(3, 3, ""), (1, 1, "A"), (6, 6, "")] if i < 3 else [(10, 10, "")]
        rows.append(spec_alignment([(20, 0, False, edits)], SEQ, name=f"n20_{i}"))
    return rows


def mixed_af_rows():
    """Node 10 (8 reads): A>C at 4 in 5, G>T at 2 in 1, 1-bp DEL at 6 in 1, INS GG at 8 in 1.
    Node 20 (5 reads): C>A at 3 in 2, G>C at 7 in 1, plus one all-match reverse-strand record."""
    rows = []
    for i in range(8):
        e = [(4, 4, ""), (1, 1, "C" if i < 4 else ""), (5, 5, "")]
        if i == 4:
            e = [(2, 2, ""), (1, 1, "T"), (7, 7, "")]
        if i == 5:
            e = [(4, 4, ""), (1, 1, "C"), (1, 1, ""), (1, 0, ""), (3, 3, "")]
        if i == 6:
            e = [(8, 8, ""), (0, 2, "GG"), (2, 2, "")]
        rows.append(spec_alignment([(10, 0, False, e)], SEQ, name=f"a{i}"))
    for i in range(5):
        e = [(3, 3, ""), (1, 1, "A"), (6, 6, "")] if i < 2 else [(10, 10, "")]
        if i == 2:
            e = [(7, 7, ""), (1, 1, "C"), (2, 2, "")]
        rows.append(spec_alignment([(20, 0, i == 4, e)], SEQ, name=f"b{i}"))
    return rows


def run(root, name, make_rows, split=True, **overrides):
    gam = root / f"{make_rows.__name__}.gam"
    if not gam.exists():
        write_gam(gam, make_rows())
        graph_fixture(root / "graph.sqlite", [(n, s, 7) for n, s in SEQ.items()])
    (root / "nodes.txt").write_text("10\n20\n")
    options = dict(gam=str(gam), nodes=str(root / "nodes.txt"), graph_index=str(root / "graph.sqlite"),
                   output=str(root / name / "shared"), batch_nodes=512, shard_size=2048)
    if split:
        options.update(snv_output=str(root / name / "SNV"), indel_output=str(root / name / "INDEL"))
    options.update(overrides)
    with redirect_stdout(io.StringIO()):
        build(build_args(**options))
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
                out = root / f"cap{cap}"
                with redirect_stdout(io.StringIO()):
                    build(build_args(**common, output=str(out), max_node_reads=cap))
                (meta,), _ = records(out)
                self.assertEqual(meta["coverage"], coverage)
                self.assertEqual(meta["parameters"]["max_node_reads"], cap)
                report = validate(out, str(gam), None, str(graph))
                self.assertTrue(report["passed"])
                self.assertEqual(report["capped_nodes"], int(cap > 0))
            # The auditor re-applies the recorded cap: a manifest that hides it must fail.
            manifest = json.loads((root / "cap5/manifest.json").read_text())
            manifest["parameters"]["max_node_reads"] = 0
            (root / "cap5/manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "independent (representative|site) coverage"):
                validate(root / "cap5", str(gam), None, str(graph))


class EarlyAfFilterTest(unittest.TestCase):
    def test_without_cap_outputs_are_unchanged_and_rejections_equal(self):
        for unit in ("site", "allele"):
            with self.subTest(unit=unit), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                common = dict(split=False, candidate_unit=unit, max_node_reads=0, min_af=.3, min_variants=1,
                              batch_nodes=2, shard_size=2, rows=4, debug_rows=True)
                base = run(root, "base", mixed_af_rows, early_af_filter=False, **common)
                fast = run(root, "fast", mixed_af_rows, early_af_filter=True, **common)
                manifest = json.loads((fast / "shared/manifest.json").read_text())
                self.assertEqual(manifest["early_af_rejected"], 4)  # G>T, DEL, INS on node 10; G>C on node 20
                self.assertEqual(manifest["early_rejected"], 0)
                shards = sorted(p.name for p in (base / "shared").glob("shard_*_data.npy"))
                self.assertEqual(shards, sorted(p.name for p in (fast / "shared").glob("shard_*_data.npy")))
                for name in shards:
                    self.assertEqual((base / "shared" / name).read_bytes(), (fast / "shared" / name).read_bytes())
                strip = lambda m: dict(m, parameters=without(m["parameters"], "early_af_filter"))  # noqa: E731
                self.assertEqual([strip(m) for m in records(base / "shared")[0]],
                                 [strip(m) for m in records(fast / "shared")[0]])
                ids = lambda f: sorted(m["candidate_id"] for m in filtered(f / "shared"))  # noqa: E731
                self.assertEqual(ids(base), ids(fast))
                early = [m for m in filtered(fast / "shared") if m.get("support_not_evaluated")]
                self.assertEqual(len(early), 4)
                for m in early:
                    self.assertEqual(m["reasons"], ["min_af"])
                    self.assertLess(m["af_upper_bound"], .3)
                    self.assertEqual(m["af_upper_bound"], m["alt_support_upper_bound"] / m["coverage"])
                    self.assertEqual("site_id" in m, unit == "site")


class SiteTest(unittest.TestCase):
    def test_site_tensor_holds_every_passing_allele_one_allele_sites_equal_allele_units(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            common = dict(snv_min_af=.2, indel_min_af=.1, min_variants=2)
            site = run(root, "site", site_rows, **common)
            allele = run(root, "allele", site_rows, candidate_unit="allele", **common)
            for kind, expected in (("SNV", {"10:4:SNV": ["A>C", "A>G"], "20:3:SNV": ["C>A"]}),
                                   ("INDEL", {"10:7:INDEL": ["DEL", "INS"]})):
                smeta, stensor = records(site / kind)
                ameta, atensor = records(allele / kind)
                got = {m["site_id"]: [a["ref"] + ">" + a["alt"] if kind == "SNV" else a["event_type"]
                                      for a in m["alleles"]] for m in smeta}
                self.assertEqual(got, expected)
                self.assertEqual(len(ameta), sum(map(len, expected.values())))
                by_id = {m["candidate_id"]: (atensor[i], m) for i, m in enumerate(ameta)}
                for i, m in enumerate(smeta):
                    self.assertEqual((m["sample_unit"], m["allele_count"]), ("site-v2", len(m["alleles"])))
                    self.assertEqual(m["candidate_id"], m["alleles"][0]["candidate_id"])
                    tensor, allele_meta = by_id[m["candidate_id"]]
                    site_only = ("sample_unit", "site_id", "alleles", "allele_count", "second_allele_af",
                                 "index_within_shard", "parameters", "shard_index")
                    if len(m["alleles"]) == 1:  # a one-allele site is exactly that allele's tensor
                        np.testing.assert_array_equal(stensor[i], tensor)
                        self.assertEqual(without(m, *site_only), without(allele_meta, *site_only))
                    else:  # every passing allele is in the site tensor; top-level counts are A1's
                        for k in ("candidate_id", "coverage", "alt_count", "ref_count", "other_count", "af"):
                            self.assertEqual(m[k], allele_meta[k])
                        self.assertEqual(list(m["allele_labels"]), [a["label"] for a in m["alleles"]])
                        self.assertEqual(m["site_counts"]["A1"], m["alleles"][0]["alt_count"])
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
                    for a in m["alleles"]:
                        self.assertEqual(without(a, "label"), {k: by_id[a["candidate_id"]][1][k] for k in a if k != "label"})
                # Every allele appears exactly once: in some site's alleles[] or in the filtered log.
                site_ids = [a["candidate_id"] for m in smeta for a in m["alleles"]] + \
                           [m["candidate_id"] for m in filtered(site / kind)]
                allele_ids = [m["candidate_id"] for m in ameta] + [m["candidate_id"] for m in filtered(allele / kind)]
                self.assertEqual(len(site_ids), len(set(site_ids)))
                self.assertEqual(sorted(site_ids), sorted(allele_ids))
                self.assertTrue(all("site_id" in m for m in filtered(site / kind)))
                self.assertEqual(validate_shards(site / kind, 2048)["sites"], len(expected))
                self.assertIsNone(validate_shards(allele / kind, 2048)["sites"])
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
        from indexed_gam_pipeline_v2.candidates import NodeReads, SiteLayout, make_site_tensor
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
        self.assertEqual(candidate_units(sorted(c), True),  # sites (1,5,INDEL) < (1,5,SNV) < (1,6,SNV);
                         [(c[0], c[1]), (c[2], c[3]), (c[4],)])  # alleles keep candidate order
        self.assertEqual(candidate_units(sorted(c), False), [(x,) for x in sorted(c)])

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
    def test_frozen_command_carries_unit_cap_and_prefilter_explicitly(self):
        builder = dict(rows=200, width=101, gam_cache_mb=8192, batch_nodes=512, max_node_span=10000,
                       max_batch_alignments=20000, shard_size=2048, min_mapq=10, min_af=.05, min_variants=3,
                       min_allele_bq=10, max_indel_len=50, chromosomes="all", candidate_unit="site", max_node_reads=800,
                       early_af_filter=True)
        config = dict(python="python", tensors="/t", variant_outputs=dict(SNV=.06, INDEL=.08), builder=builder,
                      parts=[dict(nodes_file="/n")], inputs=dict(gam=dict(path="/g"), index=dict(path="/g.gai"),
                                                               graph_index=dict(path="/x")))
        command = build_command(Path("/root"), config, 0)
        joined = " ".join(command)
        self.assertIn("--candidate-unit site --max-node-reads 800 --early-af-filter", joined)
        builder.update(candidate_unit="allele", max_node_reads=0, early_af_filter=False)
        joined = " ".join(build_command(Path("/root"), config, 0))
        self.assertIn("--candidate-unit allele --max-node-reads 0 --no-early-af-filter", joined)


if __name__ == "__main__":
    unittest.main()

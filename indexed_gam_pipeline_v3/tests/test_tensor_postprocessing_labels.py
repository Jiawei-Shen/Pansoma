"""tensor_postprocessing: GFA reference path, chromosome blocks, truth -> node keys, labels, recall."""
import contextlib
import io
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

import numpy as np
import pysam

from .fixtures import CHR1, CHR2, CHR3, graph_files, graph_fixture, write_vcf
from ..tensor_postprocessing.chr_index import AUTOSOMES, ChrIndex
from ..tensor_postprocessing.reference_path import ReferencePath, rc
from ..tensor_postprocessing.truth_labels import (Bed, LABELS, Locator, TruthSet,
                                                                        label_directory, placements, recall,
                                                                        split_alleles)
from ..tools import graph_prep
from ..tools.graph_prep import scan


class ReferencePathTest(unittest.TestCase):
    def test_coordinates_orientation_and_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            gfa, _ = graph_files(tmp)
            meta = scan(gfa, Path(tmp) / "rp")
            path = ReferencePath(Path(tmp) / "rp")
            self.assertEqual((meta["walks"], meta["ambiguous_reference_nodes"]), (4, 1))  # node 8 visited twice
            self.assertEqual(path.contigs, ["chr1", "chr2", "chr3"])
            self.assertEqual(int(path.chrom[7]), -1)
            locator = Locator(path)
            sequences = dict(line.split("\t")[1:3] for line in Path(gfa).read_text().splitlines() if line[0] == "S")
            for x in range(len(CHR1)):  # every base, every SNV: node key -> linear -> same base
                for alt in "ACGT":
                    if alt == CHR1[x]:
                        continue
                    (key,) = locator.keys("chr1", x, CHR1[x], alt, "SNP")
                    node, start, kind, alleles = key.split(":")
                    ref, a = alleles.split(">")
                    self.assertEqual(sequences[node][int(start)], ref)
                    lin = path.linear(int(node), int(start), ref, a, kind)
                    self.assertEqual((lin["chrom"], lin["pos0"], lin["ref"], lin["alt"]), ("chr1", x, CHR1[x], alt))
            self.assertEqual(locator.keys("chr2", 8, CHR2[8], "C", "SNP"), ["9:1:SNP:A>C"])  # node 9 is unique ...
            self.assertEqual(locator.keys("chr2", 2, CHR2[2], "C", "SNP"), [])  # ... node 8 (visited twice) is not
            self.assertIsNone(path.linear(8, 0, "G", "C", "SNP"))
            self.assertEqual(locator.keys("chr1", 3, "TT", "", "DEL"), ["1:3:DEL:TT>"])
            self.assertEqual(locator.keys("chr1", 4, "TT", "", "DEL"), [])  # would span nodes 1 and 2
            # Reverse node 2 (linear 5..11): linear DEL of T at 7 is forward offset 4 of rc("TTTGCAA") = "TTGCAAA".
            self.assertEqual(locator.keys("chr1", 7, "T", "", "DEL"), ["2:4:DEL:A>"])
            # Insertion at the node1|node2 junction: a key on both nodes.
            self.assertEqual(locator.keys("chr1", 5, "", "C", "INS"), ["1:5:INS:>C", "2:7:INS:>G"])
            for node, start, ref, alt, kind in ((2, 4, "A", "", "DEL"), (2, 7, "", "G", "INS"), (1, 5, "", "C", "INS")):
                lin = path.linear(node, start, ref, alt, kind)
                self.assertEqual(Locator(path).keys(lin["chrom"], lin["pos0"], lin["ref"], lin["alt"], kind).count(
                    f"{node}:{start}:{kind}:{ref}>{alt}"), 1)

    def test_chr_index_blocks_and_walk_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            gfa, _ = graph_files(tmp)
            scan(gfa, Path(tmp) / "rp")
            components = Path(tmp) / "components"
            for chrom, nodes in (("chr1", range(1, 8)), ("chr2", range(8, 10))):
                (components / chrom).mkdir(parents=True)
                (components / chrom / f"{chrom}.component.nodes.raw.txt").write_text("\n".join(map(str, nodes)) + "\n")
            blocks, meta = graph_prep.build(components, Path(tmp) / "rp", Path(tmp) / "idx.chr_node_ranges",
                                            autosomes=("chr1", "chr2"))
            self.assertEqual([(b["chrom"], b["first_node"], b["last_node"]) for b in blocks],
                             [("chr1", 1, 7), ("chr2", 8, 9), ("unplaced", 10, 11)])  # chr3 is not in this test's autosomes
            self.assertEqual(meta["walk_check"]["total"], 4)
            self.assertEqual(meta["version"], "chr-node-ranges-v1")
            index = ChrIndex(Path(tmp) / "idx.chr_node_ranges.tsv")
            self.assertEqual(index.meta, json.loads(json.dumps(meta)))
            self.assertEqual(index.names_of([7, 8, 10, 12]), ["chr1", "chr2", "unplaced", None])
            (components / "chr1" / "chr1.component.nodes.raw.txt").write_text("1\n2\n4\n")
            with self.assertRaisesRegex(ValueError, "not one contiguous"):
                graph_prep.build(components, Path(tmp) / "rp", Path(tmp) / "bad", autosomes=("chr1", "chr2"))


class GraphPrepTest(unittest.TestCase):
    """tools.graph_prep: the one-time graph preparation moved out of the runtime package (formats unchanged)."""

    def run_cli(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            graph_prep.main([str(a) for a in argv])
        return json.loads(out.getvalue())

    def test_cli_scan_check_and_chr_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gfa, fasta = graph_files(tmp)
            scanned = self.run_cli("ref-path-scan", "--gfa", gfa, "--output", tmp / "rp")
            self.assertNotIn("contigs", scanned)
            self.assertEqual((scanned["version"], scanned["max_node"], scanned["walks"]), ("gfa-reference-path-v1", 11, 4))
            self.assertEqual(ReferencePath(tmp / "rp").contigs, ["chr1", "chr2", "chr3"])
            segments = [line.split("\t")[1:3] for line in gfa.read_text().splitlines() if line[0] == "S"]
            db = graph_fixture(tmp / "g.sqlite", [(int(n), s, 1) for n, s in segments])
            with sqlite3.connect(db) as connection:  # chr-index records the index's GBZ fingerprint
                (value,) = connection.execute("SELECT value FROM graph_metadata").fetchone()
                connection.execute("UPDATE graph_metadata SET value=?",
                                   (json.dumps(dict(json.loads(value), source=dict(path="g.gbz"))),))
            report = self.run_cli("ref-path-check", "--path", tmp / "rp", "--graph-index", db, "--fasta", fasta,
                                  "--samples", 50)
            self.assertTrue(report["passed"])
            self.assertEqual((report["length_mismatches"], report["fasta_sequence_mismatches"]), (0, 0))
            self.assertEqual(ReferencePath(tmp / "rp").meta["checks"]["graph_index_and_fasta"], report)
            # the CLI uses all 22 autosomes: chr1-3 are the graph's, chr4-22 get one node each above it
            spans = dict(chr1=range(1, 8), chr2=range(8, 10), chr3=range(10, 12))
            for k, chrom in enumerate(AUTOSOMES):
                (tmp / "components" / chrom).mkdir(parents=True)
                nodes = spans.get(chrom, range(100 + k, 101 + k))
                (tmp / "components" / chrom / f"{chrom}.component.nodes.raw.txt").write_text(
                    "\n".join(map(str, nodes)) + "\n")
            result = self.run_cli("chr-index", "--components-dir", tmp / "components", "--reference-path", tmp / "rp",
                                  "--output", tmp / "idx", "--graph-index", db)
            self.assertEqual([(b["chrom"], b["first_node"], b["last_node"]) for b in result["blocks"][:3]],
                             [("chr1", 1, 7), ("chr2", 8, 9), ("chr3", 10, 11)])
            self.assertEqual(result["walk_check"]["total"], 4)
            meta = json.loads((tmp / "idx.json").read_text())
            self.assertEqual((meta["version"], meta["graph_index"]["nodes"], meta["graph_index"]["gbz"]),
                             ("chr-node-ranges-v1", 11, dict(path="g.gbz")))
            self.assertEqual(ChrIndex(tmp / "idx.tsv").names_of([1, 9, 11, 118, 121]), ["chr1", "chr2", "chr3", "chr19", "chr22"])

    def test_graph_prep_left_the_runtime_package(self):
        from ..tensor_postprocessing import __main__ as cli, chr_index, reference_path
        gone = {reference_path: ("SEPARATORS", "AWK", "parse_walk", "scan", "check"),
                chr_index: ("VERSION", "NAMED_GROUPS", "UNPLACED", "group_of", "component_block", "build", "SELECTIONS"),
                ReferencePath: ("unique",), ChrIndex: ("is_autosome",)}
        for owner, names in gone.items():
            for name in names:
                self.assertFalse(hasattr(owner, name), f"{owner.__name__}.{name}")
        for command in ("ref-path-scan", "ref-path-check", "chr-index"):
            with self.subTest(command=command), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as raised:
                    cli.main([command])
                self.assertEqual(raised.exception.code, 2)


class CrossNodeDeletionTest(unittest.TestCase):
    """Truth keys of deletions over node junctions equal the v6 builder's candidate_id (Candidate.path)."""

    def decoded(self, gfa, specs):
        from .fixtures import spec_alignment
        from ..candidates import decode_alignment
        sequences = {int(l.split("\t")[1]): l.split("\t")[2] for l in Path(gfa).read_text().splitlines() if l[0] == "S"}
        read, _ = decode_alignment(spec_alignment(specs, sequences), sequences)
        (obs,) = [o for o in read.observations if o.candidate.kind == "DEL"]
        return obs.candidate.metadata()

    def test_forward_and_reverse_node_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            gfa, _ = graph_files(tmp)
            scan(gfa, Path(tmp) / "rp")
            path = ReferencePath(Path(tmp) / "rp")
            locator = Locator(path)
            cases = [
                # chr1 nodes 5 (>, linear 44..60) and 6 (>, 60..80): delete linear 58..62.
                ("chr1", CHR1, 58, 4, [(5, 10, False, [(4, 4, ""), (2, 0, "")]), (6, 0, False, [(2, 0, ""), (6, 6, "")])]),
                # chr3 nodes 10 (<, linear 0..12) and 11 (<, 12..30), read along the walk: delete linear 10..14.
                ("chr3", CHR3, 10, 4, [(10, 6, True, [(4, 4, ""), (2, 0, "")]), (11, 0, True, [(2, 0, ""), (6, 6, "")])]),
            ]
            for chrom, genome, pos, m, specs in cases:
                with self.subTest(chrom=chrom):
                    meta = self.decoded(gfa, specs)
                    keys = {k for p in placements(genome, pos, genome[pos:pos + m], "", "DEL")
                            for k in locator.keys(chrom, *p, "DEL")}
                    self.assertIn(meta["candidate_id"], keys)
                    lin = path.linear(meta["node_id"], meta["start"], meta["ref"], meta["alt"], "DEL", meta["path"])
                    self.assertIn((lin["pos0"], lin["ref"]),
                                  [(p[0], p[1]) for p in placements(genome, pos, genome[pos:pos + m], "", "DEL")])
                    self.assertEqual(lin["chrom"], chrom)
            # The reverse pair is keyed on its forward-first node (11, the linear right one), REF reverse-complemented.
            self.assertEqual(locator.keys("chr3", 10, CHR3[10:14], "", "DEL"), [f"11:16:DEL:{rc(CHR3[10:14])}>@10"])
            self.assertEqual(locator.keys("chr1", 58, CHR1[58:62], "", "DEL"), [f"5:14:DEL:{CHR1[58:62]}>@6"])
            # Mixed orientation (node 1 >, node 2 <): no key, and linear() refuses such a path.
            self.assertEqual(locator.keys("chr1", 4, "TT", "", "DEL"), [])
            self.assertIsNone(path.linear(1, 4, "TT", "", "DEL", [[2, 6, 7]]))
            self.assertIsNone(path.linear(5, 14, "ACGT", "", "DEL", [[3, 0, 2]]))  # not the neighbouring node


class TruthLabelTest(unittest.TestCase):
    def test_split_and_placements(self):
        self.assertEqual(list(split_alleles(10, "GGCC", ("G", "GGCA", "TGCC"))),
                         [(10, "GCC", "", "DEL", 1), (12, "C", "A", "SNP", 2), (9, "G", "T", "SNP", 3)])
        self.assertEqual(list(split_alleles(5, "C", ("CAA",))), [(5, "", "AA", "INS", 1)])
        genome = "GATATATC"
        self.assertEqual(placements(genome, 5, "", "AT", "INS"),
                         [(1, "", "AT"), (2, "", "TA"), (3, "", "AT"), (4, "", "TA"), (5, "", "AT"), (6, "", "TA"),
                          (7, "", "AT")])
        self.assertEqual(placements("GTTTC", 3, "T", "", "DEL"), [(1, "T", ""), (2, "T", ""), (3, "T", "")])

    def test_labels_for_every_rule_and_recall(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gfa, fasta = graph_files(tmp)
            scan(gfa, tmp / "rp")
            path = ReferencePath(tmp / "rp")
            locator = Locator(path)
            # somatic: 1-bp DEL in the T run (VCF anchored at the G before it), SNV at 30.
            # somatic: 1-bp DEL in the T run, SNV at 30, and a filtered SNV at 60 (not PASS -> ignore).
            s_alt = "A" if CHR1[60] != "A" else "C"
            somatic_vcf = write_vcf(tmp / "somatic.vcf", [(3, "GT", "G", "PASS", "0|1"),
                                                          (31, CHR1[30], "A" if CHR1[30] != "A" else "C", "PASS", "0|1"),
                                                          (61, CHR1[60], s_alt, "LowQual", "0|1")])
            g_alt = "A" if CHR1[45] != "A" else "C"
            f_alt = "A" if CHR1[50] != "A" else "C"
            o_alt = "A" if CHR1[72] != "A" else "C"
            germline_vcf = write_vcf(tmp / "germline.vcf", [(46, CHR1[45], g_alt, "PASS", "1|1"),
                                                            (51, CHR1[50], f_alt, "GAP1", "1|0"),
                                                            (73, CHR1[72], o_alt, "PASS", "0|1")])  # PASS, outside germline BED
            (tmp / "somatic.bed").write_text("chr1\t0\t80\n")
            (tmp / "germline.bed").write_text("chr1\t0\t70\n")
            somatic = TruthSet("somatic", somatic_vcf, tmp / "somatic.bed", fasta, locator, chromosomes=("chr1",))
            germline = TruthSet("germline", germline_vcf, tmp / "germline.bed", fasta, locator, chromosomes=("chr1",))
            # The DEL has five placements across nodes 1 and 2; all four single-node ones are keys.
            (deletion, snv, _) = somatic.alleles
            self.assertEqual(deletion["placements"], 5)
            self.assertEqual(deletion["keys"], sorted(["1:3:DEL:T>", "1:4:DEL:T>", "2:4:DEL:A>", "2:5:DEL:A>", "2:6:DEL:A>"]))
            confident = somatic.bed.intersect(germline.bed)

            def candidate(x, alt, kind="SNP", ref=None):
                ref = CHR1[x] if ref is None else ref
                (key,) = locator.keys("chr1", x, ref, alt, kind)
                node, start, kind, alleles = key.split(":")
                r, a = alleles.split(">")
                return dict(candidate_id=key, node_id=int(node), start=int(start), ref=r, alt=a, event_type=kind)

            def other(x):
                return "G" if CHR1[x] not in "G" else "T"

            vcf_del = candidate(5, "", "DEL", "T")        # vg placed it on reverse node 2: still the truth allele
            snv_alt = [a for a in "ACGT" if a not in (CHR1[30], snv["alt"])][0]
            cases = [
                ("somatic_del", vcf_del, [], 1, "representative_allele_is_somatic_truth"),
                ("germline", candidate(45, g_alt), [], 2, "representative_allele_is_germline_truth"),
                ("germline_filtered", candidate(50, f_alt), [], -1, "germline_truth_filtered"),
                ("germline_outside_bed", candidate(72, o_alt), [], 2, "representative_allele_is_germline_truth"),
                ("somatic_filtered", candidate(60, s_alt), [], -1, "somatic_truth_filtered"),
                ("other_allele", candidate(30, snv_alt), [candidate(30, snv["alt"])], -1,
                 "truth_matches_non_representative_allele"),
                ("near", candidate(33, other(33)), [], -1, "near_truth_allele_mismatch"),
                ("non", candidate(18, other(18)), [], 0, "confident_no_truth_allele"),
                ("outside", candidate(75, other(75)), [], -1, "outside_confident_region"),
                ("off_path", dict(candidate_id="7:0:SNP:G>A", node_id=7, start=0, ref="G", alt="A", event_type="SNP"),
                 [], -1, "not_on_unique_grch38_node"),
            ]
            merged = tmp / "SNV"
            merged.mkdir()
            with (merged / "chr1_variant_summary.ndjson").open("w") as out:
                for i, (_, rep, extra, _, _) in enumerate(cases):
                    alleles = [rep] + extra
                    out.write(json.dumps(dict(rep, site_id=f"s{i}", alleles=alleles, shard_file="chr1_shard_00000_data.npy",
                                              shard_index=0, index_within_shard=i)) + "\n")
            (merged / "manifest.json").write_text(json.dumps(dict(layout="chromosome-shards-v1", chromosomes={
                "chr1": dict(summary="chr1_variant_summary.ndjson",
                             shards=[dict(file="chr1_shard_00000_data.npy", tensors=len(cases))])})))
            report, matched = label_directory(merged, somatic, germline, path, confident, {})
            labels = np.load(merged / "chr1_shard_00000_labels.npy")
            self.assertEqual(labels.dtype, np.int8)
            self.assertEqual(labels.tolist(), [c[3] for c in cases])
            lines = [json.loads(l) for l in (merged / "chr1_labels.ndjson").read_text().splitlines()]
            for (name, _, _, value, reason), line in zip(cases, lines):
                self.assertEqual((line["label"], line["reason"]), (value, reason), name)
            self.assertEqual(lines[0]["somatic"][0]["vcf_pos"], 3)
            self.assertEqual(report["totals"], {"somatic": 1, "germline": 2, "ignore": 6, "non": 1})
            self.assertEqual(LABELS, {"ignore": -1, "non": 0, "somatic": 1, "germline": 2})
            summary = recall(somatic, matched["somatic"], {}, tmp)
            self.assertEqual(summary["status"], {"tensor_representative": 2, "tensor_non_representative_allele": 1})
            self.assertTrue((tmp / "somatic.recall.tsv").exists())

    def test_insertion_near_span_covers_the_whole_repeat(self):
        """An insertion in a long homopolymer is 'near' anywhere in the run, not only at its leftmost placement
        (regression: residual edits at the run's right end were labelled non)."""
        class NoKeys:
            def keys(self, *args):
                return []
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "g.fa").write_text(">chr1\nGC" + "A" * 30 + "GTCAGTCAGTCAGTCAGT\n")
            pysam.faidx(str(tmp / "g.fa"))
            (tmp / "b.bed").write_text("chr1\t0\t50\n")
            vcf = write_vcf(tmp / "t.vcf", [(2, "C", "CA", "PASS", "0|1")])
            truth = TruthSet("somatic", vcf, tmp / "b.bed", tmp / "g.fa", NoKeys(), chromosomes=("chr1",))
            self.assertEqual(truth.alleles[0]["placements"], 31)  # boundaries 2..32
            self.assertTrue(truth.near("chr1", 31, 33))            # right end of the run, 29 bp from the left one
            self.assertFalse(truth.near("chr1", 44, 45))           # 12 bp past the run

    def test_bed_intersection_and_containment(self):
        a = Bed({"chr1": [(0, 10), (8, 20), (30, 40)]})
        b = Bed({"chr1": [(5, 35)], "chr2": [(0, 5)]})
        both = a.intersect(b)
        self.assertEqual([(int(s), int(e)) for s, e in zip(*both.data["chr1"])], [(5, 20), (30, 35)])
        self.assertTrue(both.contains("chr1", 5, 20) and not both.contains("chr1", 19, 21) and not both.contains("chr2", 0, 1))


if __name__ == "__main__":
    unittest.main()

"""tensor_postprocessing: GFA reference path, chromosome blocks, truth -> node keys, labels, recall."""
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest

import numpy as np
import pysam

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from indexed_gam_pipeline_v2.tensor_postprocessing import chr_index  # noqa: E402
from indexed_gam_pipeline_v2.tensor_postprocessing.reference_path import ReferencePath, rc, scan  # noqa: E402
from indexed_gam_pipeline_v2.tensor_postprocessing.truth_labels import (Bed, LABELS, Locator, TruthSet,  # noqa: E402
                                                                        label_directory, placements, recall,
                                                                        split_alleles)

rng = random.Random(5)
# chr1: 80 bp with a T run at 3..7 crossing the node1 | node2 junction (node 2 is reverse-oriented).
CHR1 = "ACGTTTTTGCAACACACGTAGGCT" + "".join(rng.choice("ACGT") for _ in range(56))
PIECES = [(1, 5, False), (2, 7, True), (3, 12, False), (4, 20, True), (5, 16, False), (6, 20, False)]
CHR2 = "GATTACAGATTACA"


def graph_files(directory):
    """GFA (chr1 walk over nodes 1-6, chr2 walk 8>9>8, HG1 walk through off-reference node 7) and FASTA."""
    directory = Path(directory)
    segments, walk, cursor = [], "", 0
    for node, length, reverse in PIECES:
        piece = CHR1[cursor:cursor + length]
        segments.append((node, rc(piece) if reverse else piece))
        walk += ("<" if reverse else ">") + str(node)
        cursor += length
    assert cursor == len(CHR1)
    segments += [(7, "GG"), (8, CHR2[:7]), (9, CHR2[7:])]
    lines = ["H\tVN:Z:1.1"] + [f"S\t{n}\t{s}" for n, s in segments]
    lines += [f"W\tGRCh38\t0\tchr1\t0\t{len(CHR1)}\t{walk}",
              f"W\tGRCh38\t0\tchr2\t0\t{len(CHR2) + 7}\t>8>9>8",
              "W\tHG1\t1\tctg1\t0\t9\t>1>7>3"]
    (directory / "g.gfa").write_text("\n".join(lines) + "\n")
    (directory / "g.fa").write_text(f">chr1\n{CHR1}\n>chr2\n{CHR2}{CHR2[:7]}\n")
    pysam.faidx(str(directory / "g.fa"))
    return directory / "g.gfa", directory / "g.fa"


def write_vcf(path, records):
    header = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=80>", "##contig=<ID=chr2,length=21>",
              '##FILTER=<ID=GAP1,Description="gap">', '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
    body = [f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t.\tGT\t{gt}" for pos, ref, alt, filt, gt in records]
    Path(path).write_text("\n".join(header + body) + "\n")
    return path


class ReferencePathTest(unittest.TestCase):
    def test_coordinates_orientation_and_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            gfa, _ = graph_files(tmp)
            meta = scan(gfa, Path(tmp) / "rp")
            path = ReferencePath(Path(tmp) / "rp")
            self.assertEqual((meta["walks"], meta["ambiguous_reference_nodes"]), (3, 1))  # node 8 visited twice
            self.assertEqual(path.contigs, ["chr1", "chr2"])
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
            blocks, meta = chr_index.build(components, Path(tmp) / "rp", Path(tmp) / "idx.chr_node_ranges",
                                           autosomes=("chr1", "chr2"))
            self.assertEqual([(b["chrom"], b["first_node"], b["last_node"]) for b in blocks], [("chr1", 1, 7), ("chr2", 8, 9)])
            self.assertEqual(meta["walk_check"]["total"], 3)
            index = chr_index.ChrIndex(Path(tmp) / "idx.chr_node_ranges.tsv")
            self.assertEqual(index.names_of([7, 8, 10]), ["chr1", "chr2", None])
            (components / "chr1" / "chr1.component.nodes.raw.txt").write_text("1\n2\n4\n")
            with self.assertRaisesRegex(ValueError, "not one contiguous"):
                chr_index.build(components, Path(tmp) / "rp", Path(tmp) / "bad", autosomes=("chr1", "chr2"))


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
            somatic_vcf = write_vcf(tmp / "somatic.vcf", [(3, "GT", "G", "PASS", "0|1"),
                                                          (31, CHR1[30], "A" if CHR1[30] != "A" else "C", "PASS", "0|1")])
            g_alt = "A" if CHR1[45] != "A" else "C"
            f_alt = "A" if CHR1[50] != "A" else "C"
            germline_vcf = write_vcf(tmp / "germline.vcf", [(46, CHR1[45], g_alt, "PASS", "1|1"),
                                                            (51, CHR1[50], f_alt, "GAP1", "1|0")])
            (tmp / "somatic.bed").write_text("chr1\t0\t80\n")
            (tmp / "germline.bed").write_text("chr1\t0\t70\n")
            somatic = TruthSet("somatic", somatic_vcf, tmp / "somatic.bed", fasta, locator, chromosomes=("chr1",))
            germline = TruthSet("germline", germline_vcf, tmp / "germline.bed", fasta, locator, chromosomes=("chr1",))
            # The DEL has five placements across nodes 1 and 2; all four single-node ones are keys.
            (deletion, snv) = somatic.alleles
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
                ("germline_filtered", candidate(50, f_alt), [], -1, "germline_truth_filtered_or_outside_bed"),
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
            self.assertEqual(report["totals"], {"somatic": 1, "germline": 1, "ignore": 5, "non": 1})
            self.assertEqual(LABELS, {"ignore": -1, "non": 0, "somatic": 1, "germline": 2})
            summary = recall(somatic, matched["somatic"], {}, tmp)
            self.assertEqual(summary["status"], {"tensor_representative": 1, "tensor_non_representative_allele": 1})
            self.assertTrue((tmp / "somatic.recall.tsv").exists())

    def test_bed_intersection_and_containment(self):
        a = Bed({"chr1": [(0, 10), (8, 20), (30, 40)]})
        b = Bed({"chr1": [(5, 35)], "chr2": [(0, 5)]})
        both = a.intersect(b)
        self.assertEqual([(int(s), int(e)) for s, e in zip(*both.data["chr1"])], [(5, 20), (30, 35)])
        self.assertTrue(both.contains("chr1", 5, 20) and not both.contains("chr1", 19, 21) and not both.contains("chr2", 0, 1))


if __name__ == "__main__":
    unittest.main()

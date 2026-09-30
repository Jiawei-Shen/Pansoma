"""graph_vcf -> linear_vcf -> vcfeval on a small graph: a GRCh38 walk over nodes 1 2 4 5 6 7, off-reference node 3
(the other allele of node 2), and one tensor per case: SNPs, insertions and deletions in and at the ends of
nodes, a deletion over two nodes, the same insertion on both sides of a node junction, homopolymer shifts,
off-reference and non-contiguous events. The rtg vcfeval test runs when $RTG or rtg on PATH exists."""
import gzip
import json
import os
import shutil
import sqlite3
import subprocess
import tempfile
import unittest
from pathlib import Path

import pysam
import torch

from .. import graph_vcf, linear_vcf, vcfeval
from ..vcfio import import_pipeline, parse_info, read_vcf

import_pipeline()
from indexed_gam_pipeline_v4.graph_index import METRIC, SCHEMA  # noqa: E402
from indexed_gam_pipeline_v4.tools.graph_prep import scan  # noqa: E402

NODES = {1: "ACGTACGTAA", 2: "C", 3: "T", 4: "GGGTTTACGA", 5: "T", 6: "A", 7: "CCATGCATGCAAAAACGT"}
WALK = (1, 2, 4, 5, 6, 7)
CHR1 = "".join(NODES[n] for n in WALK)  # node starts on chr1: 1 at 0, 2 at 10, 4 at 11, 5 at 21, 6 at 22, 7 at 23
RP = linear_vcf.graph_path

# (candidate_id, node, start, ref, alt, kind, path, af, p_somatic)
SNV = [("4:6:SNP:A>C", 4, 6, "A", "C", "SNP", [], 0.40, 0.9),         # chr1 18 A>C, PASS
       ("1:2:SNP:G>T", 1, 2, "G", "T", "SNP", [], 0.30, 0.001),       # below --min-score: no record
       ("1:4:SNP:A>G", 1, 4, "A", "G", "SNP", [], 0.05, 0.95),        # LowAF
       ("7:0:SNP:C>A", 7, 0, "C", "A", "SNP", [], 0.30, 0.3),         # LowQual (threshold 0.5), chr1 24 C>A
       ("3:0:SNP:T>G", 3, 0, "T", "G", "SNP", [], 0.50, 0.8)]         # off-reference node: unplaced
INDEL = [("7:5:INS:>G", 7, 5, "", "G", "INS", [], 0.3, 0.9),          # left-aligned 1 bp: chr1 27 T>TG
         ("4:0:INS:>A", 4, 0, "", "A", "INS", [], 0.3, 0.8),          # node offset 0: chr1 11 C>CA ...
         ("2:1:INS:>A", 2, 1, "", "A", "INS", [], 0.2, 0.7),          # ... the same insertion at node 2's end
         ("7:13:DEL:A>", 7, 13, "A", "", "DEL", [], 0.3, 0.9),        # homopolymer, 3 bp left: chr1 33 CA>C
         ("4:0:DEL:GG>", 4, 0, "GG", "", "DEL", [], 0.3, 0.6),        # offset 0: graph padding after it
         ("5:0:DEL:TA>@6", 5, 0, "TA", "", "DEL", [[6, 0, 1]], 0.3, 0.9),  # two nodes, ends at a node end
         ("5:0:DEL:TC>@7", 5, 0, "TC", "", "DEL", [[7, 0, 1]], 0.3, 0.9)]  # skips node 6: no GRCh38 position


def columns(path):
    return [r for r in read_vcf(path)[2]]


def write_checkpoint(path, at_recall=None, labels=(dict(snv_min_af=0.07, indel_min_af=None),)):
    """The parts of a pansoma_net_v2 checkpoint graph_vcf reads: validation thresholds at fixed recalls and the
    training sets' labels."""
    at_recall = {"0.9": dict(threshold=0.5), "0.8": dict(threshold=0.85)} if at_recall is None else at_recall
    torch.save(dict(val=dict(somatic_at_recall=at_recall), data=[dict(directory="set", labels=dict(x)) for x in labels]),
               path)
    return path


class Fixture:
    def __init__(self, root, fasta_bases=CHR1):
        self.root = Path(root)
        gfa = self.root / "g.gfa"
        gfa.write_text("\n".join(["H\tVN:Z:1.1", *(f"S\t{n}\t{s}" for n, s in NODES.items()),
                                  f"W\tGRCh38\t0\tchr1\t0\t{len(CHR1)}\t" + "".join(f">{n}" for n in WALK),
                                  "W\tHG1\t1\tctg1\t0\t21\t>1>3>4"]) + "\n")
        self.reference_path = self.root / "rp"
        scan(gfa, self.reference_path)
        self.fasta = self.root / "g.fa"
        self.fasta.write_text(f">chr1\n{fasta_bases}\n")
        pysam.faidx(str(self.fasta))
        self.graph_index = self.root / "graph.sqlite"
        with sqlite3.connect(self.graph_index) as db:
            db.execute("CREATE TABLE graph_metadata(value TEXT)")
            db.execute("INSERT INTO graph_metadata VALUES(?)",
                       (json.dumps(dict(schema=SCHEMA, metric=METRIC, status="complete", nodes=len(NODES))),))
            db.execute("CREATE TABLE nodes(node_id INTEGER PRIMARY KEY, seq TEXT, distinct_path_count INTEGER)")
            db.executemany("INSERT INTO nodes VALUES(?,?,?)", [(n, s, 1) for n, s in NODES.items()])
        self.tensors = self.root / "sample" / "tensors"
        self.checkpoint = write_checkpoint(self.root / "best.pth")  # threshold 0.5 at recall 0.9, SNV floor 0.07
        self.predictions = {}
        for kind, cases in (("SNV", SNV), ("INDEL", INDEL)):
            (self.tensors / kind).mkdir(parents=True)
            with open(self.tensors / kind / "chr1_variant_summary.ndjson", "w") as f:
                for cid, node, start, ref, alt, event, path, af, _ in cases:
                    record = dict(candidate_id=cid, node_id=node, start=start, ref=ref, alt=alt,
                                  event_type=event, path=path, coverage=20, alt_count=round(20 * af),
                                  ref_count=20 - round(20 * af), af=af, site_id=f"{node}:{start}:{kind}",
                                  allele_count=1, chrom="chr1")
                    if cid == "7:13:DEL:A>":
                        record["downsampled_from"] = 5000  # a node over the read cap
                    f.write(json.dumps(record) + "\n")
            run = self.root / f"run_{kind}"
            run.mkdir()
            self.predictions[kind] = run / f"sample.tensors.{kind}.predictions.ndjson.gz"
            (run / f"sample.tensors.{kind}.metrics.json").write_text(json.dumps(  # predict's: best-F1 0.5
                dict(threshold=0.5, checkpoint=str(self.checkpoint), labels=dict(snv_min_af=0.06))))
            with gzip.open(self.predictions[kind], "wt") as f:
                for case in cases:
                    p = case[-1]
                    f.write(json.dumps(dict(chrom="chr1", candidate_id=case[0], p_non=round(1 - p, 5), p_somatic=p,
                                            p_germline=0.0, pred="somatic" if p >= 0.5 else "non")) + "\n")

    def graph(self, kind, **options):
        out = self.root / f"{kind}.graph.vcf.gz"
        argv = ["--predictions", str(self.predictions[kind]), "--tensors", str(self.tensors), "--kind", kind,
                "--graph-index", str(self.graph_index), "--output", str(out)]
        for key, value in options.items():
            argv += [f"--{key.replace('_', '-')}", str(value)]
        return graph_vcf.build(graph_vcf.parse_args(argv)), out

    def linear(self, kind, **options):
        _, graph = self.graph(kind, **options)
        out = self.root / f"{kind}.linear.vcf.gz"
        stats = linear_vcf.convert(linear_vcf.parse_args(
            ["--graph-vcf", str(graph), "--reference-path", str(self.reference_path), "--fasta", str(self.fasta),
             "--output", str(out)]))
        return stats, out


class GraphVcfTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.f = Fixture(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_snv_records_filters_and_score(self):
        stats, out = self.f.graph("SNV")
        records = {r[2]: r for r in columns(out)}
        self.assertEqual(sorted(records), ["1:4:SNP:A>G", "3:0:SNP:T>G", "4:6:SNP:A>C", "7:0:SNP:C>A"])
        r = records["4:6:SNP:A>C"]
        self.assertEqual(r[:5] + [r[6]], ["4", "7", "4:6:SNP:A>C", "A", "C", "PASS"])
        self.assertEqual(r[5], "10.000")  # -10 log10(1 - 0.9)
        self.assertEqual(r[8:], ["GT:DP:AD:AF", "0/1:20:12,8:0.4"])
        info = parse_info(r[7])
        self.assertEqual((info["START"], info["KIND"], info["EVREF"], info["EVALT"], info["P_SOMATIC"]),
                         ("6", "SNP", "A", "C", "0.90000"))
        self.assertEqual(records["1:4:SNP:A>G"][6], "LowAF")
        self.assertEqual(records["7:0:SNP:C>A"][6], "LowQual")
        self.assertEqual((stats["threshold"], stats["min_af"], stats["records"]), (0.5, 0.07, 4))
        self.assertEqual((stats["min_af_source"], stats["call"]),
                         ("checkpoint training labels snv_min_af", "p_somatic >= 0.5 (validation recall 0.9)"))
        self.assertEqual(stats["filters"], {"PASS": 2, "LowQual": 1, "LowAF": 1})
        header = read_vcf(out)[0]
        self.assertIn("##contig=<ID=4,length=10>", header)
        self.assertIn("##pansoma_chromosomes=chr1", header)

    def test_indel_padding_in_node_coordinates(self):
        _, out = self.f.graph("INDEL")
        got = {r[2]: (r[0], r[1], r[3], r[4], "NOPAD" in parse_info(r[7])) for r in columns(out)}
        self.assertEqual(got, {
            "7:5:INS:>G": ("7", "5", "G", "GG", False),         # the node base before the insertion
            "4:0:INS:>A": ("4", "1", "G", "AG", False),         # offset 0: the base after it, POS 1
            "2:1:INS:>A": ("2", "1", "C", "CA", False),
            "7:13:DEL:A>": ("7", "13", "AA", "A", False),
            "4:0:DEL:GG>": ("4", "1", "GGG", "G", False),
            "5:0:DEL:TA>@6": ("5", "1", "TAN", "N", True),      # ends at node 6's end: no base known after it
            "5:0:DEL:TC>@7": ("5", "1", "TCC", "C", False)})    # REF runs on over node 7
        self.assertEqual(parse_info(next(r for r in columns(out) if r[2] == "5:0:DEL:TA>@6")[7])["PATH"], "6:0:1")
        self.assertEqual(parse_info(next(r for r in columns(out) if r[2] == "7:13:DEL:A>")[7])["NODE_RECORDS"], "5000")

    def test_threshold_min_score_and_min_af_options(self):
        stats, out = self.f.graph("SNV", threshold=0.95, min_score=0.0005, min_af=0.0)
        filters = {r[2]: r[6] for r in columns(out)}
        self.assertEqual(filters, {"1:2:SNP:G>T": "LowQual", "1:4:SNP:A>G": "PASS", "4:6:SNP:A>C": "LowQual",
                                   "7:0:SNP:C>A": "LowQual", "3:0:SNP:T>G": "LowQual"})
        self.assertEqual(stats["pred_differs_from_threshold"], 2)  # predict's pred used 0.5

    def test_threshold_predict_is_predicts_call_set(self):
        """With --threshold predict a record is PASS when predict called it (pred, made before rounding), whatever
        its rounded p_somatic; by default, when p_somatic reaches the validation threshold."""
        lines = [json.loads(x) for x in gzip.open(self.f.predictions["SNV"], "rt")]
        lines[3].update(p_somatic=0.5, pred="non")      # 7:0 rounded up to the threshold, below it before
        lines[0].update(p_somatic=0.49999, pred="somatic")
        with gzip.open(self.f.predictions["SNV"], "wt") as f:
            f.write("".join(json.dumps(x) + "\n" for x in lines))
        stats, out = self.f.graph("SNV", threshold="predict")
        self.assertEqual({r[2]: r[6] for r in columns(out)}["4:6:SNP:A>C"], "PASS")
        self.assertEqual({r[2]: r[6] for r in columns(out)}["7:0:SNP:C>A"], "LowQual")
        self.assertEqual((stats["pred_differs_from_threshold"], stats["threshold"]), (2, 0.5))
        _, out = self.f.graph("SNV")
        self.assertEqual({r[2]: r[6] for r in columns(out)}["4:6:SNP:A>C"], "LowQual")
        self.assertEqual({r[2]: r[6] for r in columns(out)}["7:0:SNP:C>A"], "PASS")

    def test_offref_rescue_is_called_by_offref_call(self):
        """A record with predict's off-reference rescue fields is called by offref_call, not the threshold, and takes
        QUAL and OFFREF_P_* from the rescored probabilities; the other records keep the threshold."""
        lines = [json.loads(x) for x in gzip.open(self.f.predictions["SNV"], "rt")]
        lines[4].update(offref_call=False, p_offref_non=0.7, p_offref_somatic=0.2, p_offref_germline=0.1)  # 3:0, p 0.8
        lines[1].update(offref_call=True, p_offref_non=0.3, p_offref_somatic=0.6, p_offref_germline=0.1)   # 1:2, p 0.001
        with gzip.open(self.f.predictions["SNV"], "wt") as f:
            f.write("".join(json.dumps(x) + "\n" for x in lines))
        stats, out = self.f.graph("SNV")
        records = {r[2]: r for r in columns(out)}
        self.assertEqual({k: r[6] for k, r in records.items()}, {"1:2:SNP:G>T": "PASS", "1:4:SNP:A>G": "LowAF",
                                                               "3:0:SNP:T>G": "LowQual", "4:6:SNP:A>C": "PASS",
                                                               "7:0:SNP:C>A": "LowQual"})
        self.assertEqual(records["1:2:SNP:G>T"][5], "3.979")  # -10 log10(1 - 0.6): the rescored p_somatic
        info = parse_info(records["1:2:SNP:G>T"][7])
        self.assertEqual((info["P_SOMATIC"], info["OFFREF_P_SOMATIC"], info["OFFREF_P_NON"]), ("0.00100", "0.60000", "0.30000"))
        self.assertNotIn("OFFREF_P_SOMATIC", parse_info(records["4:6:SNP:A>C"][7]))
        self.assertEqual((stats["offref_rescored"], stats["offref_called"]), (2, 1))
        self.assertIn("offref_call (argmax)", stats["call"])

    def test_threshold_and_min_af_follow_the_checkpoint(self):
        stats, out = self.f.graph("SNV", target_recall=0.8)  # threshold 0.85
        self.assertEqual({r[2]: r[6] for r in columns(out)}, {"1:4:SNP:A>G": "LowAF", "3:0:SNP:T>G": "LowQual",
                                                              "4:6:SNP:A>C": "PASS", "7:0:SNP:C>A": "LowQual"})
        # the floor is the training labels' (0.07), not the labels predict's test used (0.06 in .metrics.json)
        self.assertEqual((stats["threshold"], stats["min_af"]), (0.85, 0.07))
        self.assertIn(f"##pansoma_checkpoint={self.f.checkpoint.resolve()}", read_vcf(out)[0])
        other = write_checkpoint(self.f.root / "other.pth", labels=[dict(snv_min_af=None)])  # --checkpoint wins
        stats, out = self.f.graph("SNV", checkpoint=other)
        self.assertEqual((stats["min_af"], {r[2]: r[6] for r in columns(out)}["1:4:SNP:A>G"]), (None, "PASS"))
        write_checkpoint(self.f.checkpoint, at_recall={})  # a checkpoint from before val.somatic_at_recall
        with self.assertRaisesRegex(SystemExit, "no validation threshold at recall 0.9"):
            self.f.graph("SNV")
        self.assertEqual(self.f.graph("SNV", threshold=0.5)[0]["filters"], {"PASS": 2, "LowQual": 1, "LowAF": 1})
        write_checkpoint(self.f.checkpoint, labels=[dict(snv_min_af=0.07), dict(snv_min_af=None)])
        with self.assertRaisesRegex(SystemExit, "no single snv_min_af"):
            self.f.graph("SNV")
        self.assertEqual(self.f.graph("SNV", min_af=0.07)[0]["min_af"], 0.07)
        report = self.f.predictions["SNV"].with_name("sample.tensors.SNV.metrics.json")
        report.write_text(json.dumps(dict(threshold=0.5)))  # no checkpoint to read
        with self.assertRaisesRegex(SystemExit, "no --checkpoint"):
            self.f.graph("SNV")
        self.assertIsNone(self.f.graph("SNV", threshold=0.5, min_af=0.07)[0]["checkpoint"])

    def test_predictions_must_match_the_summaries(self):
        lines = gzip.open(self.f.predictions["SNV"], "rt").read().splitlines()
        with gzip.open(self.f.predictions["SNV"], "wt") as f:
            f.write("\n".join([lines[1], lines[0]] + lines[2:]) + "\n")
        with self.assertRaisesRegex(ValueError, "line 1"):
            self.f.graph("SNV")

    def test_event_must_fit_the_graph(self):
        with open(self.f.tensors / "SNV" / "chr1_variant_summary.ndjson") as f:
            lines = f.read().replace('"ref": "A", "alt": "C"', '"ref": "G", "alt": "C"')
        (self.f.tensors / "SNV" / "chr1_variant_summary.ndjson").write_text(lines)
        with self.assertRaisesRegex(ValueError, "does not fit the graph"):
            self.f.graph("SNV")


class LinearVcfTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.f = Fixture(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_snv_projection_and_unplaced(self):
        stats, out = self.f.linear("SNV")
        got = [(r[0], r[1], r[2], r[3], r[4], r[6]) for r in columns(out)]
        self.assertEqual(got, [("chr1", "5", "1:4:SNP:A>G", "A", "G", "LowAF"), ("chr1", "18", "4:6:SNP:A>C", "A", "C", "PASS"),
                               ("chr1", "24", "7:0:SNP:C>A", "C", "A", "LowQual")])
        for r in columns(out):
            self.assertEqual(CHR1[int(r[1]) - 1], r[3])
        self.assertEqual([r[2] for r in columns(self.f.root / "SNV.linear.unplaced.vcf.gz")], ["3:0:SNP:T>G"])
        self.assertEqual((stats["unplaced"], stats["unplaced_pass"], stats["records"]), (1, 1, 3))
        self.assertIn("##contig=<ID=chr1,length=41>", read_vcf(out)[0])

    def test_indel_left_alignment_padding_and_merge(self):
        stats, out = self.f.linear("INDEL")
        got = {r[2]: (r[1], r[3], r[4], parse_info(r[7]).get("LEFTSHIFT"), parse_info(r[7]).get("MERGED"))
               for r in columns(out)}
        self.assertEqual(got, {
            "4:0:INS:>A;2:1:INS:>A": ("11", "C", "CA", None, "2"),   # one allele from both sides of the junction
            "4:0:DEL:GG>": ("11", "CGG", "C", None, None),
            "5:0:DEL:TA>@6": ("20", "GAT", "G", "1", None),          # TA at 21 = AT at 20
            "7:5:INS:>G": ("27", "T", "TG", "1", None),
            "7:13:DEL:A>": ("33", "CA", "C", "3", None)})
        for r in columns(out):  # REF is the reference; ALT keeps REF's padding base
            pos, ref, alt = int(r[1]) - 1, r[3], r[4]
            self.assertEqual(CHR1[pos:pos + len(ref)], ref)
            self.assertEqual(alt[0], ref[0])
        self.assertEqual([r[2] for r in columns(self.f.root / "INDEL.linear.unplaced.vcf.gz")], ["5:0:DEL:TC>@7"])
        self.assertEqual((stats["merged"], stats["left_shifted"], stats["unplaced"]), (1, 3, 1))
        merged = next(r for r in columns(out) if r[1] == "11" and r[4] == "CA")
        self.assertEqual(parse_info(merged[7])["P_SOMATIC"], "0.80000")  # the more probable record stays
        self.assertEqual(parse_info(merged[7])["MERGED_ALT"], "10")      # 6 + 4 A1 reads
        self.assertEqual(merged[9], "0/1:20:14,6:0.3")                   # FORMAT: the kept record's

    def test_merge_keeps_a_pass_record(self):
        lines = [json.loads(x) for x in gzip.open(self.f.predictions["INDEL"], "rt")]
        lines[2].update(p_somatic=0.85, pred="somatic")  # 2:1:INS:>A: more probable, but LowAF (AF 0.2)
        with gzip.open(self.f.predictions["INDEL"], "wt") as f:
            f.write("".join(json.dumps(x) + "\n" for x in lines))
        _, out = self.f.linear("INDEL", min_af=0.25)
        merged = next(r for r in columns(out) if r[1] == "11" and r[4] == "CA")
        self.assertEqual((merged[2], merged[6]), ("4:0:INS:>A;2:1:INS:>A", "PASS"))

    def test_projection_equals_the_pipeline(self):
        """Every placed record's event, before left-alignment, is ReferencePath.linear of its summary."""
        from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath
        path = ReferencePath(self.f.reference_path)
        _, graph = self.f.graph("INDEL")
        for r in columns(graph):
            info = parse_info(r[7])
            case = next(c for c in INDEL if c[0] == r[2])
            self.assertEqual(path.linear(int(r[0]), int(info["START"]), info.get("EVREF", ""), info.get("EVALT", ""),
                                         info["KIND"], RP(info.get("PATH"))),
                             path.linear(case[1], case[2], case[3], case[4], case[5], case[6]))

    def test_reference_mismatch_stops(self):
        root = Path(self.tmp.name) / "bad"
        root.mkdir()
        f = Fixture(root, fasta_bases=CHR1[:17] + "T" + CHR1[18:])  # chr1 18 is A in the graph
        with self.assertRaisesRegex(SystemExit, "disagree"):
            f.linear("SNV")


RTG = os.environ.get("RTG") or shutil.which("rtg")


def write_vcf_gz(path, header_contigs, records, sample=False):
    lines = ["##fileformat=VCFv4.2", *header_contigs,
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + ("\tFORMAT\tS" if sample else "")]
    lines += ["\t".join(r) for r in records]
    Path(str(path)[:-3]).write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(path)[:-3], preset="vcf", force=True)
    return path


class VcfevalPartsTest(unittest.TestCase):
    def test_record_kinds(self):
        self.assertEqual(vcfeval.record_kinds("A", "C"), {"SNV"})
        self.assertEqual(vcfeval.record_kinds("ACCCC", "A,ACC"), {"INDEL"})
        self.assertEqual(vcfeval.record_kinds("A", "AT,G"), {"SNV", "INDEL"})
        self.assertEqual(vcfeval.record_kinds("A", "<DEL>"), set())

    def test_summary_without_baseline_and_nan(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "summary.txt").write_text("0 total baseline variants, no summary statistics available\n")
            self.assertEqual(vcfeval.read_summary(tmp), [])
            (Path(tmp) / "summary.txt").write_text("Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  "
                                                   "Precision  Sensitivity  F-measure\n" + "-" * 20 + "\n"
                                                   "     None      0      0      0    702        NaN     0.0000        NaN\n")
            [row] = vcfeval.read_summary(tmp)
            self.assertEqual((row["fn"], row["precision"], row["f1"]), (702, None, None))

    def test_describe_picks_the_highest_threshold_reaching_each_recall(self):
        roc = [dict(score=s, p_somatic=p, tp_baseline=t, fp=fp, tp_call=t, fn=10 - t, precision=t / (t + fp),
                    recall=t / 10, f1=2 * t / (t + fp + 10))
               for s, p, t, fp in ((30, 0.999, 4, 0), (20, 0.99, 5, 1), (10, 0.9, 8, 4), (3, 0.5, 9, 20))]
        summary = [dict(threshold=None, tp_baseline=5, tp_call=5, fp=1, fn=5, precision=5 / 6, recall=0.5, f1=0.625)]
        curve = [dict(threshold=None, tp_baseline=9, tp_call=9, fp=20, fn=1, precision=9 / 29, recall=0.9, f1=0.46)]
        d = vcfeval.describe(dict(summary=summary, roc=[]), dict(summary=curve, roc=roc))
        self.assertEqual((d["baseline"], d["ceiling"], d["curve_records"]), (10, 0.9, 29))
        self.assertEqual(d["best"]["score"], 10)
        self.assertEqual({k: v and v["score"] for k, v in d["precision_at_recall"].items()},
                         {"0.5": 20, "0.8": 10, "0.9": 3, "0.95": None})


@unittest.skipUnless(RTG, "rtg not found (set $RTG)")
class VcfevalTest(unittest.TestCase):
    """PoN tagging and rtg vcfeval on the fixture's SNV calls: chr1 18 A>C PASS (a truth), chr1 24 C>A LowQual
    (a truth, in CoLoRSdb at AF 0.01), chr1 5 A>G LowAF (in 1000G); a third truth has no call."""

    def test_calls_curve_and_pon(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = Fixture(tmp)
            _, calls = f.linear("SNV")
            sdf = Path(tmp) / "g.sdf"
            subprocess.run([RTG, "format", "-o", str(sdf), str(f.fasta)], check=True, capture_output=True)
            contig = ["##contig=<ID=chr1,length=41>"]
            gt = '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">'
            truth = write_vcf_gz(Path(tmp) / "truth.vcf.gz", [gt],   # no ##contig, VAF undeclared (as COLO829T)
                                 [("chr1", "3", ".", "G", "T", ".", ".", "VAF=0.3", "GT", "0|1"),
                                  ("chr1", "18", ".", "A", "C", ".", ".", "VAF=0.4", "GT", "1|0"),
                                  ("chr1", "24", ".", "C", "A", ".", "PASS", ".", "GT", "0|1"),
                                  ("chr1", "30", ".", "A", "T", ".", "Filtered", ".", "GT", "0|1"),  # left out
                                  ("chr1", "33", ".", "CA", "C", ".", ".", ".", "GT", "0|1")], sample=True)
            af = '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">'
            sao = '##INFO=<ID=SAO,Number=1,Type=Integer,Description="SAO">'
            pons = [write_vcf_gz(Path(tmp) / f"{name}.vcf.gz", contig + [info], records) for name, info, records in (
                ("gnomad", af, [("chr1", "18", ".", "A", "G,C", ".", ".", "AF=0.2,0.00005"),  # A>C below AF 0.0001
                                ("chr1", "24", ".", "C", "T", ".", ".", "AF=0.3")]),       # other allele
                ("dbsnp", sao, [("chr1", "18", ".", "A", "C", ".", ".", "SAO=2"),          # somatic: no match
                                ("chr1", "30", ".", "A", "T", ".", ".", "SAO=0")]),
                ("1000g", af, [("chr1", "5", ".", "A", "G", ".", ".", ".")]),              # the LowAF record
                ("colors", af, [("chr1", "24", ".", "C", "A", ".", ".", "AF=0.01")]))]     # the LowQual truth
            args = vcfeval.parse_args(["--calls", str(calls), "--pon", *map(str, pons), "--truth", str(truth),
                                       "--sdf", str(sdf), "--rtg", RTG, "--rtg-mem", "1g", "--threads", "1",
                                       "--output", str(Path(tmp) / "eval")])
            report = vcfeval.evaluate(args)
            self.assertEqual((report["truth_records"], report["truth_not_pass"]), (3, 1))  # the PASS SNV truth records
            self.assertEqual({k: report["truth_pon_tags"][k] for k in ("records", "tagged", "PoN4_CoLoRSdb")},
                             dict(records=3, tagged=1, PoN4_CoLoRSdb=1))  # chr1 24
            self.assertEqual(report["unplaced_pass"], 1)
            text = (Path(tmp) / "eval" / "report.txt").read_text()
            self.assertIn("PoN tags 1 of 3 truth records", text)
            self.assertIn("1 PASS calls without a GRCh38 position are not evaluated", text)
            self.assertEqual({k: report["pon_tags"][k] for k in ("records", "tagged", "tagged_pass", "PoN3_1000G", "PoN4_CoLoRSdb")},
                             dict(records=3, tagged=2, tagged_pass=0, PoN3_1000G=1, PoN4_CoLoRSdb=1))
            raw, pon = report["raw"], report["pon"]
            a = raw["at_threshold"]
            self.assertEqual((a["tp_baseline"], a["fp"], a["fn"], raw["baseline"]), (1, 0, 2, 3))
            self.assertAlmostEqual(raw["ceiling"], 2 / 3, places=3)   # the LowQual truth is on the curve
            self.assertEqual(raw["curve_records"], 2)                # PASS + LowQual, not LowAF
            self.assertAlmostEqual(pon["ceiling"], 1 / 3, places=3)   # ... and tagged by CoLoRSdb
            self.assertEqual(pon["at_threshold"]["tp_baseline"], 1)
            again = vcfeval.parse_args(["--calls", str(calls), "--pon-vcf", report["pon_vcf"], "--truth", str(truth),
                                        "--sdf", str(sdf), "--rtg", RTG, "--rtg-mem", "1g", "--threads", "1",
                                        "--output", str(Path(tmp) / "eval2")])
            second = vcfeval.evaluate(again)
            self.assertEqual((second["pon"]["ceiling"], second["pon_files"], second["truth_pon_tags"]),
                             (pon["ceiling"], report["pon_files"], report["truth_pon_tags"]))
            stale = Path(tmp) / "stale.pon.vcf.gz"   # a PoN VCF of other calls is refused
            header, cols, records = read_vcf(report["pon_vcf"])
            from ..vcfio import write_vcf
            write_vcf(stale, header, cols, list(records)[:1])
            with self.assertRaisesRegex(SystemExit, "is not the PoN-tagged"):
                vcfeval.evaluate(vcfeval.parse_args(["--calls", str(calls), "--pon-vcf", str(stale), "--truth", str(truth),
                                                     "--sdf", str(sdf), "--rtg", RTG, "--output", str(Path(tmp) / "e3")]))
            plain = vcfeval.evaluate(vcfeval.parse_args([   # no PoN (INDELs for now): the raw calls only
                "--calls", str(calls), "--truth", str(truth), "--sdf", str(sdf), "--rtg", RTG, "--rtg-mem", "1g",
                "--threads", "1", "--output", str(Path(tmp) / "eval5")]))
            self.assertEqual(("pon" in plain, plain["pon_vcf"], plain["truth_pon_tags"]), (False, None, None))
            self.assertEqual(plain["raw"]["at_threshold"], raw["at_threshold"])
            text = (Path(tmp) / "eval5" / "report.txt").read_text()
            self.assertNotIn("\npon ", text)
            self.assertNotIn("PoN tags", text)
            self.assertIn("with them as false calls: raw 0.500\n", text + "\n")  # 1 TP, 1 unplaced
            with self.assertRaisesRegex(SystemExit, "no PASS INDEL truth"):
                vcfeval.evaluate(vcfeval.parse_args(["--calls", str(calls), "--pon-vcf", report["pon_vcf"], "--kind", "INDEL",
                                                     "--truth", str(Path(tmp) / "eval" / "truth.SNV.vcf.gz"), "--sdf", str(sdf),
                                                     "--rtg", RTG, "--output", str(Path(tmp) / "e4")]))


if __name__ == "__main__":
    unittest.main()

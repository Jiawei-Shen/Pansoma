import gzip
import json
import tempfile
import unittest
from pathlib import Path

from .. import combine

TRUTH = [(1, "chr1", "SNP", True), (2, "chr1", "INS", True), (3, "chr1", "DEL", True), (4, "chr1", "SNP", True),
         (5, "chr1", "SNP", False), (6, "chr2", "SNP", True)]  # (id, chrom, kind, in BED)


def rec(test_label, pred, truth_ids=(), in_test=True):
    return dict(chrom="chr1", in_test=in_test, test_label=test_label if in_test else None, pred=pred,
                truth_ids=list(truth_ids))


SNV = [rec(1, "somatic", [1]), rec(1, "somatic", [1]),   # two calls of one truth
       rec(1, "somatic", [2]),                            # an SNV tensor that partially matches the INS truth
       rec(0, "somatic"), rec(2, "non"),
       rec(1, "somatic", [5]),                            # its truth does not count (outside the BED)
       rec(1, "somatic", [4], in_test=False)]             # not in the test
INDEL = [rec(1, "somatic", [2]),                          # the same INS truth again, from the INDEL set
         rec(1, "non", [3]), rec(0, "somatic")]


class CombineTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name) / "sample" / "v3_tensors"
        root.mkdir(parents=True)
        with open(root / "somatic.recall.tsv", "w") as f:
            f.write("truth_id\tchrom\tvcf_pos\tkind\tpassed\tin_bed\tstatus\tdetail\n")
            for tid, chrom, kind, in_bed in TRUTH:
                f.write(f"{tid}\t{chrom}\t{100 * tid}\t{kind}\tTrue\t{in_bed}\ttensor\t\n")
        self.files = {}
        for kind, records in (("SNV", SNV), ("INDEL", INDEL)):
            d = Path(self.tmp.name) / f"run_{kind}" / "test_chr1"
            d.mkdir(parents=True)
            (d / f"sample.v3_tensors.{kind}.metrics.json").write_text(json.dumps(dict(directory=str(root / kind))))
            with gzip.open(d / f"sample.v3_tensors.{kind}.predictions.ndjson.gz", "wt") as f:
                f.write("".join(json.dumps(r) + "\n" for r in records))
            self.files[kind] = d.parent

    def tearDown(self):
        self.tmp.cleanup()

    def test_both_sets_count_every_truth_once(self):
        [r] = combine.main([str(self.files["SNV"]), str(self.files["INDEL"])])
        self.assertEqual((r["truth_alleles"], r["with_tensor"], r["tp"], r["fp"], r["fn"]), (4, 3, 2, 2, 2))
        self.assertAlmostEqual(r["precision"], 0.5)
        self.assertAlmostEqual(r["recall"], 0.5)
        self.assertEqual(r["found_by"], {"SNV": 1, "INDEL+SNV": 1})
        self.assertEqual(r["repeated_calls"], 2)   # the second call of truth 1, the INDEL call of truth 2
        self.assertEqual(r["calls"]["SNV"], dict(calls=5, on_truth=3, repeated=1, false=1, on_uncounted_truth=1))
        self.assertEqual({k: (v["truth_alleles"], v["with_tensor"], v["found"]) for k, v in r["per_truth_kind"].items()},
                         {"SNP": (2, 1, 1), "INS": (1, 1, 1), "DEL": (1, 1, 0)})

    def test_one_set_counts_other_kind_truth_as_true_calls(self):
        [r] = combine.main([str(self.files["SNV"])])
        self.assertEqual((r["truth_alleles"], r["tp"], r["other_kind_tp"], r["fp"]), (2, 1, 1, 1))
        self.assertAlmostEqual(r["precision"], 2 / 3)
        self.assertAlmostEqual(r["recall"], 0.5)

    def test_two_predictions_of_one_kind_are_refused(self):
        with self.assertRaises(SystemExit):
            combine.main([str(self.files["SNV"]), str(self.files["SNV"] / "test_chr1")])


if __name__ == "__main__":
    unittest.main()

"""The once-per-graph and once-per-GAM tools: graph_prep gfa / components / audit and gam_prep sort / check.
The cases that run vg need PANSOMA_VG=/path/to/vg (scripts/use_vg.sh sets it); the others always run."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from .fixtures import graph_files, graph_fixture, simple_alignment, tiny_gam, write_gam
from ..tensor_postprocessing.reference_path import ReferencePath
from ..tools import gam_prep, graph_prep

VG = os.environ.get("PANSOMA_VG")
# Paths of the graph_files GFA per node (W lines; chr2 revisits node 8, which counts once).
PATHS = {1: 2, 2: 1, 3: 2, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1}


def segments(gfa):
    return {int(f[1]): f[2] for f in (line.split("\t") for line in Path(gfa).read_text().splitlines()) if f[0] == "S"}


def chr_table(path, firsts):
    path.write_text("chrom\tfirst_node\tlast_node\tnodes\tdataset\tsource\n"
                    + "".join(f"chr{k}\t{n}\t{n}\t1\tautosome\ttest\n" for k, n in enumerate(firsts, 1)))
    return path


class AuditTest(unittest.TestCase):
    def test_counts_and_sequences_against_the_gfa(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gfa, _ = graph_files(tmp)
            seqs = segments(gfa)
            db = graph_fixture(tmp / "g.sqlite", [(n, s, PATHS[n]) for n, s in seqs.items()])
            table = chr_table(tmp / "idx.tsv", [1, 8, 10])
            report = graph_prep.audit(db, gfa, table, tmp / "graph_audit.json", random_nodes=5)
            self.assertTrue(report["passed"])
            self.assertTrue({"1", "8", "10"} <= set(report["nodes"]))
            self.assertEqual({int(n): c for n, c in report["nodes"].items()}, {n: PATHS[n] for n in map(int, report["nodes"])})
            self.assertEqual(json.loads((tmp / "graph_audit.json").read_text())["passed"], True)
            # a wrong count or a wrong sequence fails; the report goes to .failed, not to the output
            bad = graph_fixture(tmp / "bad.sqlite", [(n, s if n != 10 else "A" * len(s), PATHS[n] + (n == 1))
                                                     for n, s in seqs.items()])
            with self.assertRaisesRegex(ValueError, "audit failed"):
                graph_prep.audit(bad, gfa, table, tmp / "bad_audit.json", random_nodes=0)
            self.assertFalse((tmp / "bad_audit.json").exists())
            failed = json.loads((tmp / "bad_audit.json.failed").read_text())
            self.assertEqual(sorted(m["node"] for m in failed["mismatches"]), [1, 10])


class GamCheckTest(unittest.TestCase):
    def test_sorted_gam_passes_and_unsorted_or_mismatched_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gam, rows = tiny_gam(tmp)  # smallest node IDs 10, 10, 10, 20, 20, 30
            report = gam_prep.check(gam)
            self.assertEqual((report["gai_version"], report["first_records"], report["first_records_in_node_order"]),
                             (1, len(rows), True))
            unsorted = write_gam(tmp / "unsorted.gam", list(reversed(rows)))
            with self.assertRaisesRegex(ValueError, "Not a sorted"):
                gam_prep.check(unsorted)
            big = write_gam(tmp / "big.gam", [simple_alignment((n,), name=f"r{n}") for n in range(10, 400)])
            with self.assertRaisesRegex(ValueError, "mismatched GAM/index"):
                gam_prep.check(gam, index=str(big) + ".gai")

    def test_sort_refuses_a_temporary_disk_without_room(self):
        with tempfile.TemporaryDirectory() as tmp:
            gam, _ = tiny_gam(tmp)
            with self.assertRaisesRegex(ValueError, "GB free, vg gamsort needs"):  # before vg runs
                gam_prep.sort(gam, Path(tmp) / "out.sorted.gam", vg="/nonexistent/vg", tmp_dir=tmp, min_free=1e15)
            self.assertFalse(list(Path(tmp).glob("out.sorted.gam*")))

    def test_stats_counts_are_parsed(self):
        text = "Total alignments: 12\nTotal primary: 10\nTotal secondary: 2\nTotal aligned: 11\nTotal perfect: 3\n" \
               "Total gapless (softclips allowed): 9\nMatches: 1500 bp\n"
        self.assertEqual(gam_prep.parse_counts(text), {"Total alignments": 12, "Total primary": 10, "Total secondary": 2,
                                                       "Total aligned": 11, "Total perfect": 3, "Matches": 1500})


@unittest.skipUnless(VG, "PANSOMA_VG not set")
class VgToolsTest(unittest.TestCase):
    """GBZ from the graph_files GFA (GRCh38 as the reference sample), then gfa, ref-path-scan, components; gamsort."""

    def gbz(self, tmp):
        gfa, _ = graph_files(tmp)
        text = gfa.read_text().replace("H\tVN:Z:1.1\n", "H\tVN:Z:1.1\tRS:Z:GRCh38\n", 1)
        (tmp / "r.gfa").write_text(text)
        subprocess.run([VG, "gbwt", "-G", str(tmp / "r.gfa"), "--gbz-format", "-g", str(tmp / "g.gbz")], check=True,
                       capture_output=True)
        return tmp / "g.gbz", gfa

    def test_gfa_keeps_the_node_ids_and_components_cover_each_chromosome(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gbz, source = self.gbz(tmp)
            info = graph_prep.to_gfa(gbz, tmp / "out.gfa", threads=2)
            self.assertTrue(info["vg"].startswith("vg version"))
            self.assertEqual(segments(tmp / "out.gfa"), segments(source))
            with self.assertRaisesRegex(ValueError, "Output exists"):
                graph_prep.to_gfa(gbz, tmp / "out.gfa")
            graph_prep.scan(tmp / "out.gfa", tmp / "rp")
            self.assertEqual(sorted(ReferencePath(tmp / "rp").contigs), ["chr1", "chr2", "chr3"])
            meta = graph_prep.components(gbz, tmp / "rp", tmp / "comp", autosomes=("chr1", "chr2", "chr3"))
            got = {c: sorted(map(int, (tmp / "comp" / c / f"{c}.component.nodes.raw.txt").read_text().split()))
                   for c in ("chr1", "chr2", "chr3")}
            self.assertEqual(got, dict(chr1=list(range(1, 8)), chr2=[8, 9], chr3=[10, 11]))  # node 7 is off-reference
            self.assertEqual((meta["chromosomes"]["chr1"]["path"], meta["chromosomes"]["chr1"]["off_reference_nodes"]),
                             ("GRCh38#0#chr1", 1))
            self.assertTrue((tmp / "comp" / "summary.json").exists() and not (tmp / "comp.tmp").exists())

    def test_sort_publishes_a_checked_gam_with_the_same_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            _, rows = tiny_gam(tmp)
            unsorted = write_gam(tmp / "unsorted.gam", list(reversed(rows)))
            result = gam_prep.sort(unsorted, tmp / "out.sorted.gam", threads=1, tmp_dir=tmp / "sort_tmp")
            self.assertTrue(result["check"]["first_records_in_node_order"])
            self.assertTrue(Path(str(tmp / "out.sorted.gam") + ".gai").exists())
            self.assertFalse(any(p.name.endswith(".tmp") for p in tmp.iterdir()))
            self.assertEqual(list((tmp / "sort_tmp").iterdir()), [])  # gamsort's temporary directory is gone
            checked = gam_prep.check(tmp / "out.sorted.gam", unsorted=unsorted)
            self.assertEqual(checked["counts"]["Total alignments"], len(rows))
            with self.assertRaisesRegex(ValueError, "Output exists"):
                gam_prep.sort(unsorted, tmp / "out.sorted.gam")


if __name__ == "__main__":
    unittest.main()

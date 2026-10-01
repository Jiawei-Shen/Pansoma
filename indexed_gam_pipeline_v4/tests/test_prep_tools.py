"""The once-per-graph and once-per-GAM tools: graph_prep gfa / components / audit and gam_prep sort / check.
The cases that run vg need PANSOMA_VG=/path/to/vg (scripts/use_vg.sh sets it); the others always run."""
import json
import os
import re
from pathlib import Path
import sqlite3
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
    def test_every_node_is_recounted_the_same_for_any_process_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gfa, _ = graph_files(tmp)
            seqs = segments(gfa)
            db = graph_fixture(tmp / "g.sqlite", [(n, s, PATHS[n]) for n, s in seqs.items()])
            table = chr_table(tmp / "idx.tsv", [1, 8, 10])
            reports = [graph_prep.audit(db, gfa, table, tmp / f"audit{k}.json", random_nodes=5, processes=k)
                       for k in (1, 2, 3, 5, gfa.stat().st_size + 3)]  # byte ranges cut lines anywhere, even empty
            for report in reports:
                self.assertTrue(report["passed"], report["checks"])
                self.assertEqual((report["nodes_compared"], report["gfa_totals"]),
                                 (11, dict(s_lines=11, paths=4, path_visits=14, beyond_index=0)))
            self.assertTrue({"1", "8", "10"} <= set(reports[0]["sequence_nodes"]))
            self.assertEqual(json.loads((tmp / "audit1.json").read_text())["passed"], True)

    def test_wrong_counts_sequences_and_outside_nodes_fail_to_a_failed_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gfa, _ = graph_files(tmp)
            seqs = segments(gfa)
            table = chr_table(tmp / "idx.tsv", [1, 8, 10])
            # node 5 is no sampled node: only the full recount sees its count; node 10's sequence is wrong
            bad = graph_fixture(tmp / "bad.sqlite", [(n, s if n != 10 else "A" * len(s), PATHS[n] + (n == 5))
                                                     for n, s in seqs.items()])
            with self.assertRaisesRegex(ValueError, "audit failed"):
                graph_prep.audit(bad, gfa, table, tmp / "bad.json", random_nodes=0, processes=2)
            self.assertFalse((tmp / "bad.json").exists())
            failed = json.loads((tmp / "bad.json.failed").read_text())
            self.assertEqual((failed["checks"]["path_count_mismatches"], failed["checks"]["sequence_mismatches"]), (1, 1))
            self.assertEqual(failed["path_count_examples"], [dict(node=5, index_count=2, gfa_paths=1)])
            # a P line counts like a W line; a path node without an index node fails
            text = gfa.read_text() + "P\tp1\t1+,7+,3-\t*\nW\tHG2\t1\tctg2\t0\t4\t>11>99\n"
            (tmp / "p.gfa").write_text(text)
            extra = {**PATHS, 1: 3, 7: 2, 3: 3, 11: 2}
            db = graph_fixture(tmp / "p.sqlite", [(n, s, extra[n]) for n, s in seqs.items()])
            with self.assertRaisesRegex(ValueError, "audit failed"):
                graph_prep.audit(db, tmp / "p.gfa", table, tmp / "p.json", random_nodes=0, processes=3)
            failed = json.loads((tmp / "p.json.failed").read_text())
            self.assertEqual(failed["checks"], dict(path_count_mismatches=0, path_nodes_outside_index=1,
                                                    s_lines_minus_index_nodes=0, sequence_mismatches=0))
            self.assertEqual(failed["gfa_totals"]["paths"], 6)


class SparseNodeIdTest(unittest.TestCase):
    """A filtered graph (e.g. HPRC v2.1 d46) keeps its nodes' IDs: fewer nodes than the largest ID, and chromosome
    components with ID gaps. ref-path-check compares the index with the GFA's segments, chr-index uses the index."""

    def files(self, tmp):
        gfa, fasta = graph_files(tmp)
        lines = []
        for line in gfa.read_text().splitlines():
            f = line.split("\t")
            if f[0] == "S":
                f[1] = str(int(f[1]) * 10)
            elif f[0] == "W":
                f[6] = re.sub(r"\d+", lambda m: str(int(m.group()) * 10), f[6])
            lines.append("\t".join(f))
        sparse = tmp / "sparse.gfa"
        sparse.write_text("\n".join(lines) + "\n")
        db = graph_fixture(tmp / "sparse.sqlite", [(n * 10, s, PATHS[n]) for n, s in segments(gfa).items()])
        with sqlite3.connect(db) as connection:  # chr-index records the index's GBZ fingerprint
            (value,) = connection.execute("SELECT value FROM graph_metadata").fetchone()
            connection.execute("UPDATE graph_metadata SET value=?",
                               (json.dumps(dict(json.loads(value), source=dict(path="sparse.gbz"))),))
        return sparse, fasta, db

    def components(self, root, lists):
        for chrom, nodes in lists.items():
            (root / chrom).mkdir(parents=True, exist_ok=True)
            (root / chrom / f"{chrom}.component.nodes.raw.txt").write_text("".join(f"{n}\n" for n in nodes))
        return root

    def test_gapped_ids_pass_with_the_graph_index_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gfa, fasta, db = self.files(tmp)
            meta = graph_prep.scan(gfa, tmp / "rp")
            self.assertEqual((meta["segments"], meta["max_node"]), (11, 110))
            report = graph_prep.check(tmp / "rp", db, fasta, samples=50)
            self.assertTrue(report["passed"])
            self.assertEqual((report["graph_index_nodes"], report["gfa_segments"], report["gfa_max_node"]), (11, 11, 110))
            chroms = ("chr1", "chr2", "chr3")
            good = self.components(tmp / "good", dict(chr1=range(10, 80, 10), chr2=[80, 90], chr3=[100, 110]))
            blocks, meta = graph_prep.build(good, tmp / "rp", tmp / "idx", graph_index=db, autosomes=chroms)
            self.assertEqual([(b["chrom"], b["first_node"], b["last_node"], b["nodes"]) for b in blocks],
                             [("chr1", 10, 70, 7), ("chr2", 80, 90, 2), ("chr3", 100, 110, 2)])
            self.assertEqual((meta["covered_nodes"], meta["uncovered_nodes"]), (11, 0))
            with self.assertRaisesRegex(ValueError, "give the graph index"):
                graph_prep.build(good, tmp / "rp", tmp / "idx2", autosomes=chroms)
            # node 30 is a chr1 node left out of chr1's list: another component's node inside chr1's interval
            bad = self.components(tmp / "bad", dict(chr1=[10, 20, 40, 50, 60, 70], chr2=[80, 90], chr3=[100, 110]))
            with self.assertRaisesRegex(ValueError, "1 graph nodes of other components"):
                graph_prep.build(bad, tmp / "rp", tmp / "idx3", graph_index=db, autosomes=chroms)


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

"""Graph index reader contract; native builder tests (tools.graph_index_build) run only when the
tools are configured."""
import os
from pathlib import Path
import sqlite3
import subprocess
import tempfile
import unittest

from .fixtures import graph_fixture  # noqa: F401
from .. import graph_index
from ..graph_index import GraphIndex
from ..tools import graph_index_build
from ..tools.graph_index_build import build_index


class GraphIndexTest(unittest.TestCase):
    def test_chunked_reads_missing_nodes_and_read_only(self):
        with tempfile.TemporaryDirectory() as d:
            p = graph_fixture(Path(d) / "index.sqlite", [(i, "AC", i % 91) for i in range(1, 2001)])
            with GraphIndex(p) as idx:
                records = idx.get_nodes(range(1, 2001))
                self.assertEqual(len(records), 2000)
                self.assertEqual(records[92], dict(sequence="AC", distinct_path_count=1))
                self.assertEqual(idx.get_nodes([91, 91])[91]["distinct_path_count"], 0)
                self.assertEqual(idx.get_nodes([]), {})
                self.assertEqual(idx.performance["queries"], 4)
                with self.assertRaisesRegex(ValueError, "missing"):
                    idx.get_nodes([3000])
                with self.assertRaises(sqlite3.OperationalError):
                    idx.db.execute("DELETE FROM nodes")
            with sqlite3.connect(p) as db:
                db.execute("UPDATE graph_metadata SET value=?", ("{}",))
            with self.assertRaisesRegex(ValueError, "Incomplete"):
                GraphIndex(p)

    def test_the_builder_lives_in_tools(self):
        """The runtime module is the reader only; the one-time builder and its C++ source sit in tools/."""
        self.assertTrue(Path(graph_index_build.__file__).with_name("gbz_graph_index.cpp").is_file())
        self.assertFalse(Path(graph_index.__file__).with_name("gbz_graph_index.cpp").exists())
        for name in ("COUNT_DEFINITION", "compile_builder", "build_index", "main"):
            self.assertFalse(hasattr(graph_index, name), name)
            self.assertTrue(hasattr(graph_index_build, name), name)
        self.assertEqual((graph_index_build.SCHEMA, graph_index_build.METRIC), (graph_index.SCHEMA, graph_index.METRIC))


@unittest.skipUnless(all(os.environ.get(k) for k in ("GBZ_TOOL", "GBZ_GRAPH_INDEX")), "native tools not configured")
class NativeBuilderTest(unittest.TestCase):
    """Set GBZ_TOOL=/path/to/gbztool and GBZ_GRAPH_INDEX=/path/to/compiled/gbz_graph_index."""

    def test_counts_and_atomic_publication(self):
        for second in (2, 2001):
            with self.subTest(second=second), tempfile.TemporaryDirectory() as d:
                root = Path(d)
                gfa, gbz, output = root / "graph.gfa", root / "graph.gbz", root / "index.sqlite"
                gfa.write_text(f"H\tVN:Z:1.1\nS\t1\tA\nS\t{second}\tC\n"
                               f"L\t1\t+\t{second}\t+\t0M\nL\t{second}\t+\t1\t+\t0M\nL\t1\t+\t1\t-\t0M\n"
                               f"W\ts\t1\tchr1\t0\t4\t>1>{second}>1<1\nW\tt\t1\tchr1\t0\t2\t<{second}<1\n")
                subprocess.run([os.environ["GBZ_TOOL"], "convert", str(gfa), str(gbz)], check=True, capture_output=True)
                result = build_index(gbz, output, os.environ["GBZ_GRAPH_INDEX"])
                self.assertEqual(result["nodes"], 2)
                with GraphIndex(output) as idx:
                    records = idx.get_nodes([1, second])
                    self.assertEqual({n: r["sequence"] for n, r in records.items()}, {1: "A", second: "C"})
                    # Path s revisits node 1 in both orientations: still one path. Two paths touch each node.
                    self.assertEqual({n: r["distinct_path_count"] for n, r in records.items()}, {1: 2, second: 2})
                with self.assertRaises(FileExistsError):
                    build_index(gbz, output, os.environ["GBZ_GRAPH_INDEX"])
                gbz.write_bytes(b"invalid")
                with self.assertRaises(subprocess.CalledProcessError):
                    build_index(gbz, root / "bad.sqlite", os.environ["GBZ_GRAPH_INDEX"])
                self.assertFalse((root / "bad.sqlite").exists())


if __name__ == "__main__":
    unittest.main()

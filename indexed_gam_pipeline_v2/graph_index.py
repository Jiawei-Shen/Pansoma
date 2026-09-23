"""Unified GBZ graph index: forward sequence and distinct GBWT path count per node.

One read-only SQLite file, built once per graph by the native helper
`gbz_graph_index.cpp`, serves both the graph-reference channel and the
path-count channel of every tensor. Build it with:

    python -m indexed_gam_pipeline_v2.graph_index compile --deps <gbwtgraph prefix> --output bin/gbz_graph_index
    python -m indexed_gam_pipeline_v2.graph_index build --gbz graph.gbz --builder bin/gbz_graph_index --output graph.sqlite
"""
import argparse
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import tempfile
import time

from indexed_gam_pipeline_v2.common import sha256_file, stamp

SCHEMA = "gbz-all-nodes-v1"
METRIC = "distinct-gbwt-paths-v1"
COUNT_DEFINITION = ("distinct logical GBWT paths visiting the node in either orientation, "
                    "including reference paths; repeated visits count once")


class GraphIndex:
    """Read-only lookup of node sequences and path counts in chunks of 900 IDs."""

    def __init__(self, path):
        self.path = Path(path).resolve()
        self.db = sqlite3.connect(self.path.as_uri() + "?mode=ro", uri=True)
        try:
            rows = self.db.execute("SELECT value FROM graph_metadata").fetchall()
            if len(rows) != 1:
                raise ValueError("Missing/ambiguous graph index metadata")
            self.metadata = json.loads(rows[0][0])
            if (self.metadata.get("schema") != SCHEMA or self.metadata.get("metric") != METRIC
                    or self.metadata.get("status") != "complete"):
                raise ValueError("Incomplete or incompatible graph index")
            self.db.execute("SELECT node_id, seq, distinct_path_count FROM nodes LIMIT 0")
            self.db.execute("BEGIN")  # one read snapshot for the whole build
        except BaseException:
            self.db.close()
            raise
        self.performance = dict(queries=0, nodes_read=0, query_seconds=0.0)

    def get_nodes(self, nodes):
        """{node_id: {'sequence', 'distinct_path_count'}}; a missing node is an error."""
        start = time.perf_counter()
        ids = sorted(set(map(int, nodes)))
        result = {}
        for offset in range(0, len(ids), 900):
            chunk = ids[offset:offset + 900]
            query = ("SELECT node_id, seq, distinct_path_count FROM nodes WHERE node_id IN ("
                     + ",".join("?" for _ in chunk) + ")")
            for node, seq, count in self.db.execute(query, chunk):
                if not seq or not isinstance(count, int) or not 0 <= count <= 2147483647:
                    raise ValueError(f"Invalid graph index record: {node}")
                result[node] = dict(sequence=seq, distinct_path_count=count)
            self.performance["queries"] += 1
        missing = set(ids) - result.keys()
        if missing:
            raise ValueError(f"Nodes missing from GBZ graph index: {sorted(missing)[:10]}")
        self.performance["nodes_read"] += len(result)
        self.performance["query_seconds"] += time.perf_counter() - start
        return result

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.db.close()


def compile_builder(deps, output, cxx="g++"):
    """Compile gbz_graph_index.cpp against an installed gbwtgraph dependency prefix."""
    deps, output = Path(deps).resolve(), Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    libs = [str(deps / "lib" / f"lib{lib}.a")
            for lib in ("gbwtgraph", "gbwt", "sdsl", "handlegraph", "divsufsort", "divsufsort64")]
    subprocess.run([cxx, "-O3", "-std=c++17", "-fopenmp", "-I" + str(deps / "include"),
                    str(Path(__file__).with_name("gbz_graph_index.cpp")),
                    "-Wl,--start-group", *libs, "-Wl,--end-group", "-pthread", "-lsqlite3",
                    "-o", str(output)], check=True)
    return output


def build_index(gbz, output, builder):
    """Build the index in a temporary directory and publish it atomically (never overwrite)."""
    start = time.perf_counter()
    gbz, output = Path(gbz).resolve(), Path(output).absolute()
    if output.exists():
        raise FileExistsError(output)
    source = stamp(gbz)
    metadata = dict(schema=SCHEMA, metric=METRIC, source=dict(source, sha256=sha256_file(gbz)),
                    count_definition=COUNT_DEFINITION)
    hash_seconds = time.perf_counter() - start
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=output.name + ".building-", dir=output.parent) as temp:
        temp = Path(temp)
        meta = temp / "metadata.json"
        meta.write_text(json.dumps(metadata))
        dbpath = temp / "index.sqlite"
        subprocess.run([str(Path(builder).resolve()), str(gbz), str(dbpath), str(meta)], check=True)
        if stamp(gbz) != source:
            raise ValueError("GBZ changed during index build")
        t = time.perf_counter()
        with GraphIndex(dbpath) as index:
            if index.db.execute("SELECT count(*) FROM nodes").fetchone()[0] != index.metadata["nodes"]:
                raise ValueError("Incomplete graph index node table")
            if index.db.execute("PRAGMA quick_check").fetchall() != [("ok",)]:
                raise ValueError("SQLite integrity check failed")
            result = index.metadata
        result["timing"].update(source_hash_seconds=hash_seconds,
                                validation_seconds=time.perf_counter() - t,
                                total_build_seconds=time.perf_counter() - start)
        with sqlite3.connect(dbpath) as db:
            db.execute("UPDATE graph_metadata SET value=?", (json.dumps(result),))
        os.link(dbpath, output)  # same filesystem: atomic and refuses to overwrite
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    c = commands.add_parser("compile", help="compile the native index builder")
    c.add_argument("--deps", required=True, help="prefix containing include/ and lib/ of gbwtgraph")
    c.add_argument("--output", required=True)
    c.add_argument("--cxx", default="g++")
    b = commands.add_parser("build", help="build the SQLite index from a GBZ")
    b.add_argument("--gbz", required=True)
    b.add_argument("--builder", required=True, help="executable produced by `compile`")
    b.add_argument("--output", required=True, help="new .sqlite path")
    args = parser.parse_args()
    if args.command == "compile":
        print(compile_builder(args.deps, args.output, args.cxx))
    else:
        print(json.dumps(build_index(args.gbz, args.output, args.builder), indent=2))


if __name__ == "__main__":
    main()

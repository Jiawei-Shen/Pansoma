"""One-time unified GBZ graph index: forward sequence and distinct GBWT path count per node.

The native helper gbz_graph_index.cpp (next to this file) builds one read-only SQLite file per
graph, read at run time by graph_index.GraphIndex for both the graph-reference channel and the
path-count channel of every tensor:

    python -m indexed_gam_pipeline_v3.tools.graph_index_build compile --deps <gbwtgraph prefix> --output bin/gbz_graph_index
    python -m indexed_gam_pipeline_v3.tools.graph_index_build build --gbz graph.gbz --builder bin/gbz_graph_index --output graph.sqlite
"""
import argparse
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import time

from ..common import sha256_file, stamp
from ..graph_index import METRIC, SCHEMA, GraphIndex

COUNT_DEFINITION = ("distinct logical GBWT paths visiting the node in either orientation, "
                    "including reference paths; repeated visits count once")


def compile_builder(deps, output, cxx=None, portable=True):
    """Compile gbz_graph_index.cpp against an installed gbwtgraph dependency prefix.

    cxx: $CXX, else g++. portable: link libstdc++ and libgcc statically when the toolchain has
    their static archives (else dynamically), so the builder needs little beyond glibc, libgomp
    and libsqlite3. Writes <output>.build.json (compiler, flags, the glibc / libstdc++ symbol
    versions and CPU extensions the binary requires) and warns when it uses AVX/BMI
    instructions: those come from dependencies built with -march=native (sdsl-lite's default
    and inherited by gbwt/gbwtgraph), and then the builder only runs on CPUs like the build
    host's -- rebuild the dependencies with generic flags for a portable builder.
    """
    from .binary_requirements import binary_requirements
    cxx = cxx or os.environ.get("CXX") or "g++"
    deps, output = Path(deps).resolve(), Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_name("gbz_graph_index.cpp")
    libs = [str(deps / "lib" / f"lib{lib}.a")
            for lib in ("gbwtgraph", "gbwt", "sdsl", "handlegraph", "divsufsort", "divsufsort64")]
    command = [cxx, "-O3", "-std=c++17", "-fopenmp", "-I" + str(deps / "include"), str(source),
               "-Wl,--start-group", *libs, "-Wl,--end-group", "-pthread", "-lsqlite3"]
    attempts = ([["-static-libstdc++", "-static-libgcc"]] if portable and sys.platform.startswith("linux") else []) + [[]]
    errors = []
    for extra in attempts:
        result = subprocess.run(command + extra + ["-o", str(output)], capture_output=True, text=True)
        if result.returncode == 0:
            sys.stderr.write(result.stderr)
            break
        errors.append(result.stderr[-4000:])
    else:
        raise RuntimeError("compiling the graph index builder failed:\n" + "\n".join(errors))
    try:
        version = subprocess.run([cxx, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    except (OSError, IndexError):
        version = None
    requirements = binary_requirements(output)
    info = dict(source_sha256=sha256_file(source), deps=str(deps), compiler=cxx, compiler_version=version,
                command=command[1:5] + extra, static_runtime=bool(extra), requirements=requirements)
    output.with_name(output.name + ".build.json").write_text(json.dumps(info, indent=2) + "\n")
    cpu = requirements.get("cpu_extensions", {})
    if any(cpu.values()):
        sys.stderr.write(f"warning: {output.name} uses CPU extensions {cpu}; its dependencies were built for this "
                         "host's CPU (-march=native) -- rebuild them with generic flags for other machines\n")
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
    c.add_argument("--cxx", help="C++17 compiler (default: $CXX, else g++)")
    c.add_argument("--no-portable", action="store_true", help="link libstdc++/libgcc dynamically")
    b = commands.add_parser("build", help="build the SQLite index from a GBZ")
    b.add_argument("--gbz", required=True)
    b.add_argument("--builder", required=True, help="executable produced by `compile`")
    b.add_argument("--output", required=True, help="new .sqlite path")
    args = parser.parse_args()
    if args.command == "compile":
        print(compile_builder(args.deps, args.output, args.cxx, portable=not args.no_portable))
    else:
        print(json.dumps(build_index(args.gbz, args.output, args.builder), indent=2))


if __name__ == "__main__":
    main()

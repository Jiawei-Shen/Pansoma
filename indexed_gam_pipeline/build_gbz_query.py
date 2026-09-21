"""Compile the GBZ query helper against an installed gbwtgraph dependency prefix."""
import argparse
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--deps', required=True, help='Prefix containing include/ and lib/ (e.g. gbz-tool/dependency)')
    parser.add_argument('--output', required=True)
    parser.add_argument('--cxx', default='g++')
    parser.add_argument('--graph-index', action='store_true', help='Compile offline unified SQLite index builder')
    args = parser.parse_args()
    deps = Path(args.deps).resolve()
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([args.cxx, '-O3', '-std=c++17', '-fopenmp',
        '-I' + str(deps / 'include'), str(Path(__file__).with_name('gbz_graph_index.cpp' if args.graph_index else 'gbz_node_counts.cpp')),
        '-Wl,--start-group', *[str(deps / 'lib' / ('lib' + lib + '.a'))
            for lib in ('gbwtgraph', 'gbwt', 'sdsl', 'handlegraph', 'divsufsort', 'divsufsort64')],
        '-Wl,--end-group', '-pthread', *(['-lsqlite3'] if args.graph_index else []), '-o', str(output)], check=True)


if __name__ == '__main__':
    main()

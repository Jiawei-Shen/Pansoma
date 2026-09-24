"""python -m indexed_gam_pipeline_v2.tensor_postprocessing <command> ...  (run from the repository root)"""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m indexed_gam_pipeline_v2.tensor_postprocessing", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    p = commands.add_parser("ref-path-scan", help="one GFA pass: node lengths, GRCh38 path coordinates, walk ranges")
    p.add_argument("--gfa", required=True)
    p.add_argument("--output", required=True, help="new directory")
    p.add_argument("--reference-sample", default="GRCh38")

    p = commands.add_parser("ref-path-check", help="check a ref-path directory against the graph index and a FASTA")
    p.add_argument("--path", required=True)
    p.add_argument("--graph-index", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--samples", type=int, default=100000)

    p = commands.add_parser("chr-index", help="node ID -> chromosome block table (chr1-22 + non-autosomal groups)")
    p.add_argument("--components-dir", required=True, help="directory with chrN/chrN.component.nodes.raw.txt")
    p.add_argument("--reference-path", required=True, help="ref-path-scan output directory")
    p.add_argument("--output", required=True, help="output prefix; writes <prefix>.tsv and <prefix>.json")
    p.add_argument("--graph-index", help="graph index SQLite (records the GBZ fingerprint, checks the node count)")

    p = commands.add_parser("merge", help="task_*/shard_* -> <chrom>_shard_* of --shard-size tensors, byte-verified")
    p.add_argument("--root", required=True, help="run root with outputs.json (orchestrate run)")
    p.add_argument("--chr-index", required=True, help="chr-index .tsv")
    p.add_argument("--shard-size", type=int, default=32768)
    p.add_argument("--keep-sources", action="store_true", help="do not delete task_* directories after verification")
    p.add_argument("--workers", type=int, default=8, help="parallel shard verification processes")
    p.add_argument("--spots", type=int, default=200, help="random tensors per chromosome reloaded from the sources")
    p.add_argument("--reference-path", help="ref-path-scan directory: adds GRCh38 coordinates (grch38) to every summary record")

    args = parser.parse_args(argv)
    if args.command == "ref-path-scan":
        from indexed_gam_pipeline_v2.tensor_postprocessing.reference_path import scan
        result = scan(args.gfa, args.output, args.reference_sample)
        result = {k: v for k, v in result.items() if k != "contigs"}
    elif args.command == "ref-path-check":
        from indexed_gam_pipeline_v2.tensor_postprocessing.reference_path import check
        result = check(args.path, args.graph_index, args.fasta, args.samples)
    elif args.command == "chr-index":
        from indexed_gam_pipeline_v2.tensor_postprocessing.chr_index import build
        blocks, meta = build(args.components_dir, args.reference_path, args.output, args.graph_index)
        result = dict(blocks=blocks, walk_check=meta["walk_check"], uncovered_nodes=meta["uncovered_nodes"])
    elif args.command == "merge":
        from indexed_gam_pipeline_v2.tensor_postprocessing.merge_shards import merge
        result = merge(args.root, args.chr_index, args.shard_size, args.keep_sources, args.workers, args.spots,
                       reference_path=args.reference_path)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

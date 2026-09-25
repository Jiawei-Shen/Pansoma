"""python -m indexed_gam_pipeline_v3.tensor_postprocessing merge|label ...  (run from the repository root)

Manual re-merge or re-label of a run root (orchestrate finalize runs both). The one-time graph
preparation (ref-path-scan, ref-path-check, chr-index) is python -m indexed_gam_pipeline_v3.tools.graph_prep.
"""
import argparse
import json


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m indexed_gam_pipeline_v3.tensor_postprocessing", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    p = commands.add_parser("merge", help="task_*/shard_* -> <chrom>_shard_* of --shard-size tensors, byte-verified")
    p.add_argument("--root", required=True, help="run root with outputs.json (orchestrate run)")
    p.add_argument("--chr-index", required=True, help="chr-index .tsv")
    p.add_argument("--shard-size", type=int, default=32768)
    p.add_argument("--keep-sources", action="store_true", help="do not delete task_* directories after verification")
    p.add_argument("--workers", type=int, default=8, help="parallel shard verification processes")
    p.add_argument("--spots", type=int, default=200, help="random tensors per chromosome reloaded from the sources")
    p.add_argument("--reference-path", help="ref-path-scan directory: adds GRCh38 coordinates (grch38) to every summary record")

    p = commands.add_parser("label", help="germline/somatic VCF -> node keys -> labels for merged directories")
    p.add_argument("--tensors", required=True, help="tensor root holding the merged <kind>/ directories")
    p.add_argument("--kinds", nargs="+", default=["SNV", "INDEL"])
    p.add_argument("--reference-path", required=True)
    p.add_argument("--fasta", required=True, help="GRCh38 FASTA (with .fai) matching the graph's GRCh38 path")
    p.add_argument("--somatic-vcf", required=True)
    p.add_argument("--somatic-bed", required=True)
    p.add_argument("--germline-vcf", required=True)
    p.add_argument("--germline-bed", required=True)
    p.add_argument("--truth-dir", required=True, help="where <set>.graph.tsv truth tables are written")
    p.add_argument("--recall-dir", help="where recall reports go (default: --tensors)")

    args = parser.parse_args(argv)
    if args.command == "merge":
        from .merge_shards import merge
        result = merge(args.root, args.chr_index, args.shard_size, args.keep_sources, args.workers, args.spots,
                       reference_path=args.reference_path)
    elif args.command == "label":
        from .truth_labels import label_run
        result = label_run(args.tensors, args.kinds, args.reference_path, args.fasta, args.somatic_vcf,
                           args.somatic_bed, args.germline_vcf, args.germline_bed, args.truth_dir, args.recall_dir)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

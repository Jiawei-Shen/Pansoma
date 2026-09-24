"""Post-processing of finished tensor runs: graph coordinates, shard merging, truth labels.

Nothing here imports the tensor builder (candidates/build), so the same code serves every
tensor format; formats are read from each output directory's manifest.

    reference_path  one pass over the GFA: node lengths, GRCh38 path coordinates, every walk's node range
    chr_index       node ID -> chromosome block (chr1-22 from connected components, plus non-autosomal groups)
    merge_shards    task_*/shard_* (2,048 each) -> <chrom>_shard_* (32,768 each), byte-verified
    truth_labels    germline/somatic VCF -> node candidate keys -> per-tensor labels

Run as `python -m indexed_gam_pipeline_v2.tensor_postprocessing <command>`.
"""

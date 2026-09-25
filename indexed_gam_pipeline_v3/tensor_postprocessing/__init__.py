"""Post-processing of finished tensor runs: chromosome blocks, shard merging, truth labels.

Nothing here imports the tensor builder (candidates/build), so the same code serves every
tensor format; formats are read from each output directory's manifest.

    reference_path  reader of the GFA reference-path directory: node lengths, GRCh38 path coordinates, walk ranges
    chr_index       node ID -> chromosome block (chr1-22 from connected components, plus non-autosomal groups)
    merge_shards    task_*/shard_* (2,048 each) -> <chrom>_shard_* (32,768 each), byte-verified
    truth_labels    germline/somatic VCF -> node candidate keys -> per-tensor labels

`orchestrate finalize` runs the merge and then the labels; for a manual re-merge or re-label run
`python -m indexed_gam_pipeline_v3.tensor_postprocessing merge|label`. The reference-path directory
and the chromosome block table are made once per graph by tools.graph_prep.
"""

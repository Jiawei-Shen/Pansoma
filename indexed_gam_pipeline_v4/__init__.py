"""Candidate-centered pileup tensors built directly from an indexed, sorted GAM (tensor format indexed-gam-candidate).

For every candidate site (node, start, SNV|INDEL) on the target nodes with a passing allele, one (8, rows, width)
int8 tensor and one summary record, in separate SNV and INDEL outputs; whole-genome runs are merged per chromosome
and truth-labelled by tensor_postprocessing. Every entry point is `python -m`, from the repository root or from a
run's frozen <root>/source. See README.md.

Runtime (frozen into every run by `orchestrate prepare`; data flows top to bottom):
    common          atomic JSON, file stamps, node lists, batching, format names (format_name)
    gam_reader      GAI bins, BGZF reads over large preads (BgzfReader), bounded LRU group cache or the record
                    path (IndexedGam), scan_gam
    gam_record_index
                    per-record index <gam>.gri of a long-read GAM: build | check | info, used by IndexedGam
    graph_index     read-only graph SQLite: node sequence + distinct GBWT path count (GraphIndex)
    candidates      GAM edits -> left-normalized candidates, support counting, (8, rows, width) site tensors
    native          optional C++ record decoder (fastdecode.cpp): compile | check, load-time checks, fallback
    build           per-task builder: fetch -> graph -> decode -> prefilter -> cap -> per-site support -> shard
    run             builder CLI: discover | build
    orchestrate     whole-genome runs: prepare | run [--resume] | task | finalize
    tensor_postprocessing
                    chromosome blocks, per-chromosome merge, truth labels (merge | label)
Tools (offline, never imported by the runtime, not frozen):
    tools.graph_index_build   compile gbz_graph_index.cpp, build the graph SQLite
    tools.graph_prep          ref-path-scan | ref-path-check | chr-index
    tools.validate_examples   independent audit of --debug-rows outputs
    tools.binary_requirements glibc / CPU-extension report of a compiled binary
    tools.compare_runs        byte and normalized comparison of two runs or build directories
Tests: tests/ (python -m unittest discover -s indexed_gam_pipeline_v4/tests -t .), goldens in tests/golden.py.
"""

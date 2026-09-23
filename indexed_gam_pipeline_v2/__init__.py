"""Candidate-centered pileup tensors built directly from an indexed, sorted GAM.

Module map (data flows top to bottom):

    gam_reader     GAI index reading/building, BGZF group seeking, LRU group cache
    graph_index    read-only GBZ node sequence + distinct GBWT path count SQLite
    candidates     decode GAM edits -> candidates, support counting, (8, rows, width) tensors
    build          per-batch builder: fetch -> graph -> decode -> prefilter -> per-site support -> encode -> shard
    run            command line: index / discover / validate / build
    orchestrate    Slurm-scale runs: freeze source, partition nodes, one process per task
    validate_examples  independent audit of debug tensors against the source GAM/graph
    inspect_tensor     human-readable text dump of a tensor folder
"""
__version__ = "2.0.0"

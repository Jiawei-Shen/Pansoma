"""Node ID -> chromosome block, as a table of contiguous node-ID intervals.

Minigraph-Cactus builds one graph per reference chromosome and then gives node IDs
chromosome by chromosome, so each chromosome's graph is one contiguous ID interval
(for HPRC v1.1 d9, in contig-name order: chr1, chr10..chr19, chr2, chr20..chr22, chr3..chr9,
chrEBV, chrM, unplaced, chrX, chrY).

* Autosome blocks (chr1-22) are exactly the node sets of `vg chunk -C -p GRCh38#0#chrN`,
  i.e. the whole connected component of the chromosome's reference path, including every
  node that is not on GRCh38 (insertions, alternative branches). Built from those lists
  (`scripts/build_chr_node_filters.sh`), which must each be one gap-free interval.
* Non-autosomal groups (chrEBV, chrM, chrX, chrY; `*_random`/`chrUn_*`/alt contigs -> "unplaced")
  take the remaining IDs up to the largest node, split at the smallest node of each group's
  reference walks.
* Every walk of every sample in the GFA must stay inside one block (a walk is one
  assembly contig); a walk leaving an autosome block is an error, because it would mean
  chromosome material with IDs outside its block.

Lookup is one np.searchsorted over the block starts.
"""
import csv
import hashlib
from pathlib import Path

import numpy as np

from ..common import read_json

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
FIELDS = ("chrom", "first_node", "last_node", "nodes", "dataset", "source")


def select_nodes(nodes, selection="all", table=None):
    """Keep the nodes of the chosen chromosome blocks: (kept int64 array, report).

    `selection` is "all" (no filtering, no table needed), "autosome" (chr1-22) or a
    comma-separated list of block names from the table (e.g. "chr1,chr2,chrX").
    """
    nodes = np.asarray(nodes, dtype=np.int64)
    if selection == "all":
        return nodes, dict(selection="all", nodes_in=int(nodes.size), nodes_kept=int(nodes.size))
    if not table:
        raise ValueError("--chromosomes other than 'all' needs --chr-index")
    index = ChrIndex(table)
    if selection == "autosome":
        wanted = [k for k, d in enumerate(index.dataset) if d == "autosome"]
    else:
        names = [n.strip() for n in selection.split(",") if n.strip()]
        unknown = sorted(set(names) - set(index.names))
        if not names or unknown:
            raise ValueError(f"Unknown chromosome blocks {unknown}; the index has {index.names}")
        wanted = [index.names.index(n) for n in names]
    blocks = index.lookup(nodes)
    keep = np.isin(blocks, wanted)
    removed = {index.names[k] if k >= 0 else "no_block": int(n)
               for k, n in zip(*np.unique(blocks[~keep], return_counts=True))}
    report = dict(selection=selection, chromosomes=[index.names[k] for k in wanted],
                  chr_index=dict(path=str(index.path), sha256=index.sha256),
                  nodes_in=int(nodes.size), nodes_kept=int(keep.sum()), removed=removed)
    return nodes[keep], report


class ChrIndex:
    """Loaded block table: lookup(node IDs) -> block index (-1 = no block)."""

    def __init__(self, table):
        table = Path(table)
        with table.open() as stream:
            blocks = list(csv.DictReader(stream, delimiter="\t"))
        self.path = table.resolve()
        self.sha256 = hashlib.sha256(table.read_bytes()).hexdigest()
        meta_path = table.with_suffix(".json")
        self.meta = read_json(meta_path) if meta_path.exists() else {}
        if self.meta and self.meta.get("tsv_sha256") != self.sha256:
            raise ValueError(f"{table} does not match its metadata checksum")
        self._set(blocks)

    @classmethod
    def from_blocks(cls, blocks):
        self = cls.__new__(cls)
        self.path, self.sha256, self.meta = None, None, {}
        self._set(blocks)
        return self

    def _set(self, blocks):
        blocks = sorted(blocks, key=lambda b: int(b["first_node"]))
        self.names = [b["chrom"] for b in blocks]
        self.dataset = [b["dataset"] for b in blocks]
        self.first = np.array([int(b["first_node"]) for b in blocks], dtype=np.int64)
        self.last = np.array([int(b["last_node"]) for b in blocks], dtype=np.int64)
        if (self.last < self.first).any() or (self.first[1:] <= self.last[:-1]).any():
            raise ValueError("Chromosome blocks must be non-empty and non-overlapping")
        if len(set(self.names)) != len(self.names):
            raise ValueError("Duplicate chromosome block names")

    def lookup(self, nodes):
        nodes = np.asarray(nodes, dtype=np.int64)
        k = np.searchsorted(self.first, nodes, side="right") - 1
        inside = (k >= 0) & (nodes <= self.last[np.maximum(k, 0)])
        return np.where(inside, k, -1)

    def names_of(self, nodes):
        return [self.names[k] if k >= 0 else None for k in self.lookup(nodes)]

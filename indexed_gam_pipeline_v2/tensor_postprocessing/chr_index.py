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

from indexed_gam_pipeline_v2.common import read_json, stamp, write_json

VERSION = "chr-node-ranges-v1"
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
NAMED_GROUPS = ("chrX", "chrY", "chrM", "chrEBV")
UNPLACED = "unplaced"
FIELDS = ("chrom", "first_node", "last_node", "nodes", "dataset", "source")


def group_of(contig):
    """Non-autosomal group of a reference contig name."""
    return contig if contig in NAMED_GROUPS else UNPLACED


def component_block(path):
    """(first, last, count) of one component node list; it must be one gap-free interval."""
    ids = np.fromfile(path, sep=" ", dtype=np.int64)
    if not ids.size:
        raise ValueError(f"Empty component node list: {path}")
    first, last = int(ids.min()), int(ids.max())
    if np.unique(ids).size != ids.size:
        raise ValueError(f"Duplicate node IDs in {path}")
    if ids.size != last - first + 1:
        raise ValueError(f"Component {path} is not one contiguous node-ID interval "
                         f"({ids.size} nodes in [{first}, {last}])")
    return first, last, int(ids.size)


def build(components_dir, reference_path, output, graph_index=None, autosomes=AUTOSOMES):
    """Write <output>.tsv + <output>.json; raises on any violated invariant."""
    from indexed_gam_pipeline_v2.tensor_postprocessing.reference_path import ReferencePath
    components_dir, output = Path(components_dir), Path(output)
    path = ReferencePath(reference_path)
    max_node = path.meta["max_node"]
    blocks, sources = [], {}
    for chrom in autosomes:
        source = components_dir / chrom / f"{chrom}.component.nodes.raw.txt"
        first, last, count = component_block(source)
        blocks.append(dict(chrom=chrom, first_node=first, last_node=last, nodes=count, dataset="autosome",
                           source=f"vg chunk -C connected component ({source.name})"))
        sources[chrom] = stamp(source)
    blocks.sort(key=lambda b: b["first_node"])
    for a, b in zip(blocks, blocks[1:]):
        if a["last_node"] >= b["first_node"]:
            raise ValueError(f"Overlapping autosome blocks: {a['chrom']} and {b['chrom']}")
    # Non-autosomal groups start at their reference contigs' smallest node.
    starts = {}
    for contig in path.meta["contigs"]:
        if contig["name"] in autosomes:
            continue
        group = group_of(contig["name"])
        starts[group] = min(starts.get(group, contig["min_node"]), contig["min_node"])
    autosome_starts = [b["first_node"] for b in blocks]
    ordered = sorted(starts.items(), key=lambda kv: kv[1])
    for i, (group, first) in enumerate(ordered):
        if any(b["first_node"] <= first <= b["last_node"] for b in blocks):
            raise ValueError(f"Reference contigs of {group} start inside an autosome block")
        following = [s for s in autosome_starts if s > first] + [s for _, s in ordered[i + 1:]]
        last = min(following) - 1 if following else max_node
        blocks.append(dict(chrom=group, first_node=first, last_node=last, nodes=last - first + 1,
                           dataset="non_autosomal", source="reference contigs' smallest node .. next block"))
    blocks.sort(key=lambda b: b["first_node"])
    index = ChrIndex.from_blocks(blocks)
    # Invariants: every reference contig inside its own block, every walk inside one block.
    for contig in path.meta["contigs"]:
        expected = contig["name"] if contig["name"] in autosomes else group_of(contig["name"])
        got = index.names_of([contig["min_node"], contig["max_node"]])
        if list(got) != [expected, expected]:
            raise ValueError(f"Reference contig {contig['name']} is not inside block {expected}: {got}")
    walks = dict(total=0, crossing_autosome=0, crossing_non_autosomal=0, unassigned=0, examples=[])
    for walk in path.walks():
        walks["total"] += 1
        a, b = index.lookup([walk["min"], walk["max"]])
        if a == b and a >= 0:
            continue
        if a < 0 or b < 0:
            walks["unassigned"] += 1
        elif "autosome" in (index.dataset[a], index.dataset[b]):
            walks["crossing_autosome"] += 1
        else:
            walks["crossing_non_autosomal"] += 1
        if len(walks["examples"]) < 20:
            walks["examples"].append(walk)
    if walks["crossing_autosome"] or walks["unassigned"]:
        raise ValueError(f"Walks leave their chromosome block: {walks}")
    covered = sum(b["nodes"] for b in blocks)
    meta = dict(version=VERSION, max_node=max_node, covered_nodes=covered, uncovered_nodes=max_node - covered,
                autosome_components=sources, reference_path=dict(path=str(Path(reference_path).resolve()),
                                                                 source=path.meta["source"]),
                walk_check=walks)
    if graph_index:
        from indexed_gam_pipeline_v2.graph_index import GraphIndex
        with GraphIndex(graph_index) as graph:
            meta["graph_index"] = dict(path=str(Path(graph_index).resolve()), nodes=graph.metadata["nodes"],
                                       gbz=graph.metadata["source"])
            if graph.metadata["nodes"] != max_node:
                raise ValueError("Graph index and GFA disagree on the node count")
    output.parent.mkdir(parents=True, exist_ok=True)
    table = output.with_name(output.name + ".tsv")
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(blocks)
    meta["tsv_sha256"] = hashlib.sha256(table.read_bytes()).hexdigest()
    write_json(output.with_name(output.name + ".json"), meta)
    return blocks, meta


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

    def is_autosome(self, k):
        return k >= 0 and self.dataset[k] == "autosome"

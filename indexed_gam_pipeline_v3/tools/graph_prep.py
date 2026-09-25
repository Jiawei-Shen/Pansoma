"""One-time graph preparation: the GFA reference-path directory and the chromosome block table.

    python -m indexed_gam_pipeline_v3.tools.graph_prep ref-path-scan --gfa G.gfa --output DIR [--reference-sample GRCh38]
    python -m indexed_gam_pipeline_v3.tools.graph_prep ref-path-check --path DIR --graph-index DB --fasta FA [--samples 100000]
    python -m indexed_gam_pipeline_v3.tools.graph_prep chr-index --components-dir D --reference-path DIR --output PREFIX [--graph-index DB]

The formats are unchanged from v2 and documented with their run-time readers:
tensor_postprocessing.reference_path (gfa-reference-path-v1 directory, ReferencePath) and
tensor_postprocessing.chr_index (chr-node-ranges-v1 TSV plus JSON with tsv_sha256, ChrIndex).

Moved verbatim from v2: reference_path SEPARATORS, AWK, parse_walk, scan and check; chr_index
VERSION (here CHR_INDEX_VERSION, because VERSION is the reference-path version that scan writes),
NAMED_GROUPS, UNPLACED, group_of, component_block and build; the three graph-prep subcommands of
tensor_postprocessing/__main__.py.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
import subprocess
import time

import numpy as np

from ..common import read_json, stamp, write_json
from ..tensor_postprocessing.chr_index import AUTOSOMES, FIELDS, ChrIndex
from ..tensor_postprocessing.reference_path import VERSION, ReferencePath, rc

# --- reference path: one pass over a GFA (v2 tensor_postprocessing/reference_path.py) -----------

SEPARATORS = bytes.maketrans(b"<>", b"  ")
# S lines -> "id<TAB>length" into LENGTHS (LN:i: when the sequence is '*'); W lines pass through.
AWK = r'''BEGIN { FS = OFS = "\t" }
$1 == "S" { n = length($3)
            if ($3 == "*") { n = -1; for (i = 4; i <= NF; i++) if ($i ~ /^LN:i:/) n = substr($i, 6) + 0 }
            print $2, n > LENGTHS; next }
$1 == "W" { print }'''


def parse_walk(walk):
    """b'>12<3' -> (ids int64, reverse bool)."""
    raw = np.frombuffer(walk, dtype=np.uint8)
    reverse = raw[(raw == 60) | (raw == 62)] == 60
    ids = np.array(walk.translate(SEPARATORS).split(), dtype=np.int64)
    if ids.size != reverse.size or not ids.size:
        raise ValueError("Malformed GFA walk")
    return ids, reverse


def scan(gfa, output, reference_sample="GRCh38"):
    """Build the output directory from one streaming pass (awk prefilter) over `gfa`."""
    started = time.perf_counter()
    output = Path(output)
    if output.exists():
        raise ValueError(f"Output exists: {output}")
    work = output.with_name(output.name + ".tmp")
    work.mkdir(parents=True)
    lengths_txt = work / "segment_lengths.txt"
    reference, walks = [], []
    process = subprocess.Popen(["awk", "-v", f"LENGTHS={lengths_txt}", AWK, str(gfa)], stdout=subprocess.PIPE)
    with (work / "walks.ndjson").open("w") as log:
        for line in process.stdout:
            fields = line.rstrip(b"\n").split(b"\t")
            if len(fields) != 7:
                raise ValueError("Malformed GFA W line")
            ids, reverse = parse_walk(fields[6])
            sample, hap, contig = (f.decode() for f in fields[1:4])
            record = dict(sample=sample, hap=hap, contig=contig, start=int(fields[4]), end=int(fields[5]),
                          nodes=int(ids.size), min=int(ids.min()), max=int(ids.max()))
            log.write(json.dumps(record) + "\n")
            walks.append(record)
            if sample == reference_sample:
                reference.append((record, ids, reverse))
    process.stdout.close()
    if process.wait():
        raise RuntimeError(f"awk failed on {gfa}")
    if not reference:
        raise ValueError(f"No W lines for sample {reference_sample} in {gfa}")
    pairs = np.fromfile(lengths_txt, sep=" ", dtype=np.int64).reshape(-1, 2)
    if (pairs[:, 1] < 0).any() or (pairs[:, 0] <= 0).any():
        raise ValueError("Segments need a positive integer ID and a sequence or LN tag")
    size = int(pairs[:, 0].max()) + 1
    lengths = np.zeros(size, dtype=np.int32)
    lengths[pairs[:, 0]] = pairs[:, 1]
    if np.count_nonzero(lengths) != len(pairs):
        raise ValueError("Duplicate or empty segments in the GFA")
    chrom = np.full(size, -1, dtype=np.int16)
    start0 = np.zeros(size, dtype=np.int64)
    reverse_of = np.zeros(size, dtype=bool)
    visits = np.zeros(size, dtype=np.int64)
    contigs, offset, order = [], 0, []
    for k, (record, ids, reverse) in enumerate(reference):
        if ids.max() >= size or not lengths[ids].all():
            raise ValueError(f"Reference walk {record['contig']} uses nodes without segments")
        span = lengths[ids].astype(np.int64)
        starts = record["start"] + np.cumsum(span) - span
        if record["start"] + int(span.sum()) != record["end"]:
            raise ValueError(f"Reference walk {record['contig']}: node lengths disagree with the W line interval")
        unique, first, counts = np.unique(ids, return_index=True, return_counts=True)
        visits[unique] += counts
        new = unique[chrom[unique] == -1]
        firsts = first[chrom[unique] == -1]
        chrom[new], start0[new], reverse_of[new] = k, starts[firsts], reverse[firsts]
        contigs.append(dict(name=record["contig"], hap=record["hap"], start=record["start"], end=record["end"],
                            nodes=int(ids.size), min_node=record["min"], max_node=record["max"], offset=offset))
        offset += ids.size
        order.append((ids, starts, reverse))
    arrays = dict(lengths=lengths, chrom=chrom, start0=start0, reverse=reverse_of,
                  visits=np.minimum(visits, np.iinfo(np.uint16).max).astype(np.uint16),
                  path_nodes=np.concatenate([o[0] for o in order]), path_starts=np.concatenate([o[1] for o in order]),
                  path_reverse=np.concatenate([o[2] for o in order]))
    for name, values in arrays.items():
        np.save(work / f"{name}.npy", values)
    lengths_txt.unlink()
    meta = dict(version=VERSION, source=stamp(gfa), reference_sample=reference_sample, max_node=size - 1,
                segments=int(len(pairs)), walks=len(walks), samples=len({w["sample"] for w in walks}),
                contigs=contigs, reference_nodes=int(np.count_nonzero(visits)),
                ambiguous_reference_nodes=int(np.count_nonzero(visits > 1)),
                scan_seconds=time.perf_counter() - started, checks={})
    write_json(work / "meta.json", meta)
    work.rename(output)
    return meta


def check(directory, graph_index, fasta, samples=100000, seed=20260923):
    """GFA lengths vs the GBZ graph index, and reference node sequences vs the FASTA; recorded in meta.json."""
    from ..graph_index import GraphIndex
    import pysam
    path = ReferencePath(directory)
    rng = random.Random(seed)
    lengths = np.asarray(path.lengths)
    present = np.flatnonzero(lengths)
    any_nodes = [int(present[rng.randrange(len(present))]) for _ in range(samples)]
    reference = np.flatnonzero(np.asarray(path.visits) == 1)
    ref_nodes = [int(reference[rng.randrange(len(reference))]) for _ in range(samples)]
    report = dict(graph_index=str(graph_index), fasta=str(fasta), samples=samples, seed=seed)
    with GraphIndex(graph_index) as graph, pysam.FastaFile(str(fasta)) as genome:
        report["graph_index_nodes"] = graph.metadata["nodes"]
        report["gfa_max_node"] = path.meta["max_node"]
        records = graph.get_nodes(set(any_nodes) | set(ref_nodes))
        report["length_mismatches"] = int(sum(len(records[n]["sequence"]) != int(lengths[n]) for n in set(any_nodes)))
        mismatches, skipped = 0, 0
        for n in sorted(set(ref_nodes)):
            contig = path.contigs[path.chrom[n]]
            if contig not in genome.references:
                skipped += 1
                continue
            start = int(path.start0[n])
            linear = genome.fetch(contig, start, start + int(lengths[n])).upper()
            node = records[n]["sequence"].upper()
            mismatches += int((rc(node) if path.reverse[n] else node) != linear)
        report.update(fasta_sequence_mismatches=mismatches, fasta_nodes_checked=len(set(ref_nodes)) - skipped,
                      fasta_contig_missing=skipped)
    report["passed"] = (report["length_mismatches"] == 0 and mismatches == 0
                        and report["graph_index_nodes"] == report["gfa_max_node"])
    meta = read_json(Path(directory) / "meta.json")
    meta["checks"]["graph_index_and_fasta"] = report
    write_json(Path(directory) / "meta.json", meta)
    if not report["passed"]:
        raise ValueError(f"Reference path check failed: {report}")
    return report


# --- chromosome block table (v2 tensor_postprocessing/chr_index.py) -----------------------------

CHR_INDEX_VERSION = "chr-node-ranges-v1"
NAMED_GROUPS = ("chrX", "chrY", "chrM", "chrEBV")
UNPLACED = "unplaced"


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
    from ..tensor_postprocessing.reference_path import ReferencePath
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
    meta = dict(version=CHR_INDEX_VERSION, max_node=max_node, covered_nodes=covered, uncovered_nodes=max_node - covered,
                autosome_components=sources, reference_path=dict(path=str(Path(reference_path).resolve()),
                                                                 source=path.meta["source"]),
                walk_check=walks)
    if graph_index:
        from ..graph_index import GraphIndex
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


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m indexed_gam_pipeline_v3.tools.graph_prep",
                                     description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
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

    args = parser.parse_args(argv)
    if args.command == "ref-path-scan":
        result = scan(args.gfa, args.output, args.reference_sample)
        result = {k: v for k, v in result.items() if k != "contigs"}
    elif args.command == "ref-path-check":
        result = check(args.path, args.graph_index, args.fasta, args.samples)
    elif args.command == "chr-index":
        blocks, meta = build(args.components_dir, args.reference_path, args.output, args.graph_index)
        result = dict(blocks=blocks, walk_check=meta["walk_check"], uncovered_nodes=meta["uncovered_nodes"])
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

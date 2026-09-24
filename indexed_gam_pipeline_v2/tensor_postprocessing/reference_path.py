"""One pass over a GFA: node lengths, GRCh38 path coordinates and the node range of every walk.

Output directory. Arrays indexed by node ID (so a lookup is one vectorized gather):
    lengths.npy        int32   segment length (0 = no such node)
    chrom.npy          int16   index into meta["contigs"] of the first reference contig visiting the node; -1 = off it
    start0.npy         int64   0-based start of that visit on the contig
    reverse.npy        bool    visited as '<' (node forward strand = contig reverse strand)
    visits.npy         uint16  number of reference visits (saturating); > 1 = ambiguous coordinate
Arrays in reference-walk order (contig k occupies [meta.contigs[k].offset, + nodes)):
    path_nodes.npy     int64   node IDs
    path_starts.npy    int64   0-based contig start of each visit (non-decreasing within a contig)
    path_reverse.npy   bool
    walks.ndjson       every W line of every sample: sample, hap, contig, start, end, nodes, min, max
    meta.json          contigs, source stamp, counts, checks (added by `check`)

A walk's node range (min, max) is enough to test chromosome-block membership, because
blocks are contiguous node-ID intervals (see chr_index).
"""
import json
from pathlib import Path
import random
import subprocess
import time

import numpy as np

from indexed_gam_pipeline_v2.common import read_json, stamp, write_json

VERSION = "gfa-reference-path-v1"
ARRAYS = ("lengths", "chrom", "start0", "reverse", "visits")
PATH_ARRAYS = ("path_nodes", "path_starts", "path_reverse")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
SEPARATORS = bytes.maketrans(b"<>", b"  ")
# S lines -> "id<TAB>length" into LENGTHS (LN:i: when the sequence is '*'); W lines pass through.
AWK = r'''BEGIN { FS = OFS = "\t" }
$1 == "S" { n = length($3)
            if ($3 == "*") { n = -1; for (i = 4; i <= NF; i++) if ($i ~ /^LN:i:/) n = substr($i, 6) + 0 }
            print $2, n > LENGTHS; next }
$1 == "W" { print }'''


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


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


class ReferencePath:
    """Memory-mapped view of a scan() directory."""

    def __init__(self, directory):
        self.directory = Path(directory)
        self.meta = read_json(self.directory / "meta.json")
        if self.meta.get("version") != VERSION:
            raise ValueError(f"Unexpected reference path version in {directory}")
        for name in ARRAYS + PATH_ARRAYS:
            setattr(self, name, np.load(self.directory / f"{name}.npy", mmap_mode="r"))
        self.contigs = [c["name"] for c in self.meta["contigs"]]
        self.contig_index = {name: k for k, name in enumerate(self.contigs)}

    def walks(self):
        with (self.directory / "walks.ndjson").open() as stream:
            for line in stream:
                yield json.loads(line)

    def contig_slice(self, name):
        c = self.meta["contigs"][self.contig_index[name]]
        return slice(c["offset"], c["offset"] + c["nodes"])

    def linear(self, node, start, ref, alt, kind, path=None):
        """GRCh38 coordinates of a node-forward candidate, or None when the node has no unique reference visit.

        Returns dict(chrom, pos0, ref, alt, node_reverse): pos0 is the 0-based start of the REF
        interval (SNP/DEL) or the insertion boundary (INS, between pos0 - 1 and pos0); REF/ALT
        are on the contig's forward strand, without anchor bases.
        """
        node = int(node)
        if not 0 < node < len(self.visits) or self.visits[node] != 1:
            return None
        contig, origin, length = int(self.chrom[node]), int(self.start0[node]), int(self.lengths[node])
        reverse = bool(self.reverse[node])
        if kind == "INS":
            pos0 = origin + (length - start if reverse else start)
        elif path:
            # A deletion over several nodes (v6 `path`: further (node, start, end) in forward order)
            # is on the reference only if those are the neighbouring reference nodes, same orientation.
            first = len(ref) - sum(e - s for _, s, e in path)
            spans = []
            for n, s, e in [(node, start, start + first)] + [tuple(p) for p in path]:
                n = int(n)
                if (not 0 < n < len(self.visits) or self.visits[n] != 1 or int(self.chrom[n]) != contig
                        or bool(self.reverse[n]) != reverse):
                    return None
                o, span = int(self.start0[n]), int(self.lengths[n])
                spans.append((o + span - e, o + span - s) if reverse else (o + s, o + e))
            spans.sort()
            if any(a[1] != b[0] for a, b in zip(spans, spans[1:])):
                return None
            pos0 = spans[0][0]
        else:
            pos0 = origin + (length - start - len(ref) if reverse else start)
        if reverse:
            ref, alt = rc(ref), rc(alt)
        return dict(chrom=self.contigs[contig], pos0=pos0, ref=ref, alt=alt, node_reverse=reverse)

    def unique(self, nodes):
        """Nodes with exactly one reference visit (usable as coordinates)."""
        nodes = np.asarray(nodes, dtype=np.int64)
        inside = (nodes > 0) & (nodes < len(self.visits))
        result = np.zeros(nodes.shape, dtype=bool)
        result[inside] = np.asarray(self.visits)[nodes[inside]] == 1
        return result


def check(directory, graph_index, fasta, samples=100000, seed=20260923):
    """GFA lengths vs the GBZ graph index, and reference node sequences vs the FASTA; recorded in meta.json."""
    from indexed_gam_pipeline_v2.graph_index import GraphIndex
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

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

import numpy as np

from ..common import read_json

VERSION = "gfa-reference-path-v1"
ARRAYS = ("lengths", "chrom", "start0", "reverse", "visits")
PATH_ARRAYS = ("path_nodes", "path_starts", "path_reverse")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


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

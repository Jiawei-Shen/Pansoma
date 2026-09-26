"""Print the REF/ALT haplotypes and the local sequence + graph path of reads at one somatic truth allele.

usage: show_locus.py CHROM VCF_POS [max_reads]
"""
import csv, sys
from collections import Counter
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reads as R  # noqa: E402  (v6 reads.py: frozen v6 decoder)

chrom, vcf_pos = sys.argv[1], int(sys.argv[2])
limit = int(sys.argv[3]) if len(sys.argv) > 3 else 12
row = next(r for r in csv.DictReader(open(R.TRUTH), delimiter="\t") if r["chrom"] == chrom and int(r["vcf_pos"]) == vcf_pos)
path = R.ReferencePath(R.REFPATH)
genome = pysam.FastaFile(R.FASTA).fetch(chrom).upper()
pos0 = int(row["pos0"])
places = R.placements(genome, pos0, row["ref"], row["alt"], row["kind"])
lo, hi = min(p[0] for p in places), max(p[0] + len(p[1]) for p in places)
pl, pr = lo - 20, hi + 20
ref_hap = genome[pl:pr]
alt_hap = genome[pl:pos0] + row["alt"] + genome[pos0 + len(row["ref"]):pr]
print(f"{chrom}:{vcf_pos} {row['vcf_ref']}>{row['vcf_alt']} GT {row['gt']}  kind {row['kind']}  span [{lo},{hi})")
print(f"  REF {ref_hap}\n  ALT {alt_hap}")
contig = path.contig_index[chrom]
part = path.contig_slice(chrom)
wn, ws = np.asarray(path.path_nodes[part]), np.asarray(path.path_starts[part])
a = max(0, int(np.searchsorted(ws, pl - 300, side="right")) - 1)
b = int(np.searchsorted(ws, pr + 300, side="right"))
nodes = {int(n) for n in wn[a:b]}
visits, chrom_of, start0, rev_of, lengths = (np.asarray(path.visits), np.asarray(path.chrom), np.asarray(path.start0),
                                             np.asarray(path.reverse), np.asarray(path.lengths))
reader = R.IndexedGam(R.GAM)
with R.GraphIndex(R.GRAPH) as graph:
    records = [x for x in reader.fetch(nodes, {}) if x.mapping_quality > R.MIN_MAPQ]
    subs = [s for s in (R.trim(x, nodes) for x in records) if s]  # >= 600 read bases each side (validate_trim.py)
    seqs = {n: r["sequence"] for n, r in graph.get_nodes({m.position.node_id for s in subs for m in s.path.mapping}).items()}
    shown = Counter()
    for s in subs:
        read, _ = R.v6.decode_alignment(s, seqs, 50, None)
        cols = read.columns
        lin = []
        for c in cols:
            on = visits[c.node] == 1 and chrom_of[c.node] == contig
            lin.append(None if not on or c.boundary else int(start0[c.node] + (lengths[c.node] - 1 - c.pos if rev_of[c.node] else c.pos)))
        idx = {x: i for i, x in enumerate(lin) if x is not None}
        if pl - 1 not in idx or pr not in idx:
            continue
        i0, i1 = sorted((idx[pl - 1], idx[pr]))
        flip = cols[idx[pl - 1]].reverse != bool(rev_of[cols[idx[pl - 1]].node])
        between = cols[i0 + 1:i1]
        seq = "".join(c.read for c in between if c.read != "-")
        if flip:
            seq = R.rc(seq)
        pathnodes = []
        for c in (between[::-1] if flip else between):
            tag = f"{c.node}{'' if visits[c.node] == 1 and chrom_of[c.node] == contig else '*'}"
            if not pathnodes or pathnodes[-1] != tag:
                pathnodes.append(tag)
        ops = "".join(c.op for c in (between[::-1] if flip else between))
        key = "REF" if seq == ref_hap else "ALT" if seq == alt_hap else "other"
        shown[key] += 1
        if shown[key] <= limit:
            print(f"  {key:5s} {seq}  ops {ops.replace('M', '.')}  path {'>'.join(pathnodes)}")
    print("  reads:", dict(shown), "(* = node not on GRCh38)")

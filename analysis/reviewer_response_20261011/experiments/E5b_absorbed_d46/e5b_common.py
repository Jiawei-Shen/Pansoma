"""Shared paths and graph helpers of E5b (analysis code; reads the pipeline's graph indexes, writes nothing).

GRCh38 coordinates of graph nodes come from the pipeline's reference-path index (<graph>.grch38_path/*.npy, the arrays
indexed_gam_pipeline_v4 labels with); node sequences and edges come from the subgraph that s3_extract.sh cuts out of the
graph's GBZ with `vg find -N <window nodes> -c 250 -L` (every node within 250 bp of a window's GRCh38 node, and the
edges among them). A GRCh38 node is a node the GRCh38 path visits once, forward (graph_absorbed_somatic_20261001 s1/s2,
classify.py ref()).
"""
import json
import numpy as np

DATA = '/scratch/jshen/data/reviewer_experiments_20261013/E5b_absorbed_d46'
FASTA = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
GRAPHS = {
    'd9': dict(gbz='/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz',
               grch38_path='/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path',
               sqlite='/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.graph_index.sqlite',
               vg='/scratch/wzhang/bin/vg_v1.65.0', H=90),
    'd46': dict(gbz='/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/hprc-v2.1-mc-grch38.d46.gbz',
                grch38_path='/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/hprc-v2.1-mc-grch38.d46.grch38_path',
                sqlite='/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/hprc-v2.1-mc-grch38.d46.graph_index.sqlite',
                vg='/scratch/jshen/bin/vg', H=464),
}
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
FLIP = {'+': '-', '-': '+'}


class RefPath:
    """GRCh38 placement of the nodes of one graph (pipeline grch38_path arrays, memory-mapped)."""

    def __init__(self, graph, in_memory=True):
        d = GRAPHS[graph]['grch38_path']
        mm = None if in_memory else 'r'      # BeeGFS random reads are slow: load the arrays once (d46 ~5 GB)
        self.meta = json.load(open(f'{d}/meta.json'))
        self.names = [c['name'] for c in self.meta['contigs']]
        self.chrom, self.start, self.length, self.rev, self.visits = (
            np.load(f'{d}/{f}.npy', mmap_mode=mm) for f in ('chrom', 'start0', 'lengths', 'reverse', 'visits'))
        self.p_nodes, self.p_starts = np.load(f'{d}/path_nodes.npy', mmap_mode=mm), np.load(f'{d}/path_starts.npy', mmap_mode=mm)
        self.contig = {c['name']: c for c in self.meta['contigs']}

    def ref(self, n):
        """(chrom, start0, end0 inclusive) if n is a GRCh38 node (visited once, forward), else None."""
        c = int(self.chrom[n])
        if c < 0 or self.visits[n] != 1 or self.rev[n]:
            return None
        s = int(self.start[n])
        return self.names[c], s, s + int(self.length[n]) - 1

    def window_nodes(self, chrom, lo, hi):
        """GRCh38 nodes (node, start0, end0) overlapping [lo, hi] in path order; None if one is walked in reverse."""
        co = self.contig[chrom]
        a, b = co['offset'], co['offset'] + co['nodes']
        starts = self.p_starts[a:b]
        i = max(0, int(np.searchsorted(starts, lo, side='right')) - 1)
        j = int(np.searchsorted(starts, hi, side='right'))
        out = []
        for k in range(a + i, a + j):
            n = int(self.p_nodes[k])
            if self.rev[n]:
                return None
            out.append((n, int(self.p_starts[k]), int(self.p_starts[k]) + int(self.length[n]) - 1))
        return out


def event_window(fa, chrom, pos0, vref, valt):
    """graph_absorbed_somatic_20261001 s1_loci.py event_window: the union of the periodic stretches (period 1-6, length
    >= max(2p, p + 3)) touching the event, +-5 bp; 0-based half-open GRCh38."""
    lref = len(vref)
    if lref == 1 and len(valt) == 1:
        e_lo, e_hi = pos0, pos0 + 1
    else:
        e_lo, e_hi = pos0 + 1, max(pos0 + lref, pos0 + 1)
    s0 = max(0, e_lo - 300)
    seq = fa.fetch(chrom, s0, e_hi + 300).upper()
    lo, hi = e_lo - s0, e_hi - s0
    L, R = lo, hi
    for p in range(1, 7):
        i = 0
        while i < len(seq) - p:
            if seq[i] != seq[i + p]:
                i += 1
                continue
            j = i
            while j < len(seq) - p and seq[j] == seq[j + p]:
                j += 1
            a, b = i, j + p
            if b - a >= max(2 * p, p + 3) and a <= hi + 1 and b >= lo - 1:
                L, R = min(L, a), max(R, b)
            i = j + 1
    return s0 + L - 5, s0 + R + 5


class SubGraph:
    """Nodes and edges of a GFA written by vg (S / L lines only)."""

    def __init__(self, gfa):
        self.seq, self.adj = {}, {}
        with open(gfa) as f:
            for line in f:
                if line[0] == 'S':
                    _, n, s = line.rstrip('\n').split('\t')[:3]
                    self.seq[int(n)] = s
                elif line[0] == 'L':
                    _, a, ao, b, bo = line.split('\t')[:5]
                    a, b = int(a), int(b)
                    self.adj.setdefault((a, ao), set()).add((b, bo))
                    self.adj.setdefault((b, FLIP[bo]), set()).add((a, FLIP[ao]))

    def oseq(self, n, o):
        s = self.seq[n]
        return s if o == '+' else s.translate(RC)[::-1]

    def succ(self, n, o):
        return sorted(self.adj.get((n, o), ()))


def spells(G, R, chrom, first, last, lo, hi, target, max_expand=200_000, max_nodes=800):
    """Does a path (first,+) ... (last,+) whose GRCh38 nodes lie in [lo, hi] of chrom spell `target` exactly?
    Prefix-pruned DFS (the s2_graph_paths 'exact' search, one target). Returns (found, path or None, info)."""
    if first not in G.seq or last not in G.seq:
        return False, None, 'anchor_missing'
    s0 = G.oseq(first, '+')
    if not target.startswith(s0):
        return False, None, 'first_mismatch'
    stack = [((first, '+'), len(s0), ((first, '+'),))]
    nexp = 0
    while stack:
        (n, o), pos, path = stack.pop()
        if n == last:
            if o == '+' and pos == len(target):
                return True, path, ''
            continue
        if len(path) >= max_nodes:
            continue
        nexp += 1
        if nexp > max_expand:
            return False, None, 'expand_cap'
        for m, mo in G.succ(n, o):
            r = R.ref(m)
            if r and (r[0] != chrom or r[1] < lo or r[2] > hi):
                continue
            x = G.oseq(m, mo)
            if target.startswith(x, pos):
                stack.append(((m, mo), pos + len(x), path + ((m, mo),)))
    return False, None, ''


def elements(R, path):
    """Departures of a path from GRCh38: list of ('branch', n_nodes, replaced_bases) / ('skip', bases)."""
    out, last, run = [], None, 0
    for n, o in path:
        r = R.ref(n)
        if r is None:
            run += 1
            continue
        if last is not None:
            gap = r[1] - last[2] - 1
            if run:
                out.append(('branch', run, gap))
            elif gap != 0:
                out.append(('skip', gap))
        last, run = r, 0
    return out

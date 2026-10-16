"""Shared helpers of E9c (graph-native population filter) and E9d (assembly validation) for unplaced calls.

Reads only: pipeline outputs (VCFs, rtg outputs), tensor sets, graph indexes, the HPRC v1.1 raw VCF, assemblies.
Writes only under DATA (large) and EXP (scripts, small tables).
"""
import gzip
import heapq
import json
import os
import re
import sqlite3
import sys
from collections import Counter

import numpy as np

REPO = '/scratch/jshen/Github/Pansoma'
sys.path.insert(0, REPO)
from indexed_gam_pipeline_v4.tensor_postprocessing.reference_path import ReferencePath  # noqa: E402

DATA = '/scratch/jshen/data/reviewer_experiments_20261013/E9c_offref_filter_validation'
EXP = f'{REPO}/analysis/reviewer_response_20261011/experiments/E9c_offref_filter_validation'
RUNS = '/scratch/jshen/data/pansoma_net_v2_runs'
TENSORS = '/scratch/jshen/data/pansoma_v2_tensors'
GRAPHS = {
    'd9': dict(index=f'{TENSORS}/graph_index/hprc-v1.1-mc-grch38.d9.graph_index.sqlite',
               refpath=f'{TENSORS}/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path', edges=f'{DATA}/edges_d9', H=90),
    'd46': dict(index='/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/hprc-v2.1-mc-grch38.d46.graph_index.sqlite',
                refpath='/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/hprc-v2.1-mc-grch38.d46.grch38_path',
                edges=f'{DATA}/edges_d46', H=464),
}
FASTA = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
HPRC_VCF = f'{DATA}/ref/hprc-v1.1-mc-grch38.raw.vcf.gz'   # symlink to /scratch/jshen/data/hprc-v1.1-mc-grch38.raw.vcf.gz
GOLD = f'{RUNS}/assembly_offref_20260929/offref_somatic.tsv'
_Q = '/scratch/qfu/COLO829BL_DSA/ragtag_hg38'
ASSEMBLIES = {   # sample -> [(target label, role, fasta)]
    'HG008T': [('HG008T_v3.2', 'tumor', '/scratch/jshen/data/HG008_GIAB/HG008T_v3.2.fasta.gz'),
               ('HG008N_v6.2', 'normal', '/scratch/jshen/data/HG008_GIAB/HG008N_curatedv6_250714_bothhaps_polished6.2.fasta.gz')],
    'COLO829T': [('COLO829BL_hapX', 'normal', f'{_Q}/ragtag_SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFTH1B2J6-verkko_2.1_hapX_primary_normalized.fa'),
                 ('COLO829BL_hapY', 'normal', f'{_Q}/ragtag_SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFK4V33V8-verkko_2.1_hapY_primary_normalized.fa'),
                 ('COLO829BL_hapX_unloc', 'normal', f'{_Q}/ragtag_SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFTH1B2J6-verkko_2.1_hapX_unlocalized_normalized.fa'),
                 ('COLO829BL_hapY_unloc', 'normal', f'{_Q}/ragtag_SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFK4V33V8-verkko_2.1_hapY_unlocalized_normalized.fa')],
}
BASE = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 5: 'N'}
GAP = 6
COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def info(field):
    return dict(x.split('=', 1) if '=' in x else (x, True) for x in field.split(';'))


def vcf_records(path):
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] != '#':
                yield line.rstrip('\n').split('\t')


def vcf_header(path):
    out = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            if not line.startswith('##'):
                break
            m = re.match(r'##(pansoma_\w+)=(.*)', line.rstrip('\n'))
            if m:
                out[m.group(1)] = m.group(2)
    return out


def graph_of(index_path):
    if 'd46' in index_path:
        return 'd46'
    if '.d9.' in index_path:
        return 'd9'
    return 'linear' if 'linear' in index_path else None


def sample_of(tensor_dir):
    name = tensor_dir.rstrip('/').split('/')[-2]
    for s in ASSEMBLIES:
        if name.startswith(s):
            return s
    return name.split('_')[0]


# ---------------------------------------------------------------- tensor contexts

def _consensus(vals, cov, min_rows=2):
    """Per column: the majority code of the covered rows (> half of them), else None."""
    out = []
    for c in range(vals.shape[1]):
        v = vals[cov[:, c], c]
        if len(v) < min_rows:
            out.append(None)
            continue
        u, n = np.unique(v, return_counts=True)
        k = int(n.argmax())
        out.append(int(u[k]) if 2 * n[k] > len(v) else None)
    return out


def _spell(cons, lo, hi, step, limit):
    """Bases of a consensus walking from column lo towards hi (exclusive) by step; stops at a column without a
    majority; gap columns (6) add nothing. Returns the bases in left-to-right order."""
    out, c, cols = [], lo, 0
    while (c < hi if step > 0 else c > hi) and cols < limit:
        v = cons[c]
        if v is None or v == 0:
            break
        if v in BASE:
            out.append(BASE[v])
        c += step
        cols += 1
    return ''.join(out[::-1] if step < 0 else out)


def infer_layout(x):
    """Row blocks and site columns from the tensor alone (equal to variant_summary's row_groups / site_layout for the
    A1 and REF blocks): channel 2 is non-zero exactly over the site columns, A1 is the first block (row 0 onwards,
    identical channel-2 site pattern), REF rows carry the graph bases there (channel 2 == channel 5)."""
    nz = np.flatnonzero((x[2] != 0).any(axis=0))
    if not len(nz):
        return None
    s0, s1 = int(nz.min()), int(nz.max()) + 1
    site = x[2, :, s0:s1]
    p0 = site[0]
    a1 = 0
    while a1 < x.shape[1] and (site[a1] == p0).all() and p0.any():
        a1 += 1
    ref = [r for r in range(a1, x.shape[1]) if site[r].any() and (site[r] == x[5, r, s0:s1]).all()]
    return dict(row_groups=[dict(allele='A1', start_row=0, end_row=a1)] +
                ([dict(allele='REF', start_row=ref[0], end_row=ref[-1] + 1)] if ref else []),
                anchor_column=s0, site_layout=dict(insertion_slots=0, span=s1 - s0))


def tensor_contexts(x, summary, flank=50):
    """Contexts around the site, node-forward, from one tensor (8, rows, width) and its variant_summary record.

    alt_read  consensus read bases of the A1 rows (the sample's ALT haplotype, with its own nearby alleles)
    ref_read  the same for the REF rows ('' when < 2 REF rows cover the site)
    alt_graph graph bases of the A1 rows' path (channel 5) with the site allele (channel 2) in the site columns
    ref_graph graph bases of the A1 rows' path including the site span (the REF path)
    Each with the site interval [s, e) in the string (e == s for a deletion).
    """
    if 'row_groups' not in summary:
        summary = infer_layout(x)
        if summary is None:
            return None
    groups = {g['allele']: np.arange(g['start_row'], g['end_row']) for g in summary['row_groups']}
    lay = summary['site_layout']
    s0 = summary.get('anchor_column', x.shape[2] // 2)
    s1 = s0 + lay['insertion_slots'] + lay['span']
    out = {}
    a1 = groups.get('A1')
    if a1 is None or len(a1) == 0:
        return None
    cov = x[4, a1] > 0
    read = _consensus(x[0, a1], cov)
    graph = _consensus(x[5, a1], cov)
    allele = _consensus(x[2, a1], cov)
    width = x.shape[2]

    def assemble(cons_flank, site_cons):
        if any(site_cons[c] is None for c in range(s0, s1)):
            return None
        left = _spell(cons_flank, s0 - 1, -1, -1, flank)
        site = ''.join(BASE[site_cons[c]] for c in range(s0, s1) if site_cons[c] in BASE)
        right = _spell(cons_flank, s1, width, 1, flank)
        return left + site + right, len(left), len(left) + len(site)

    out['alt_read'] = assemble(read, read)
    out['alt_graph'] = assemble(graph, allele)
    out['ref_graph'] = assemble(graph, graph) if s1 > s0 else None
    if out['ref_graph'] is None and out['alt_graph']:   # an insertion: the REF path is the flanks only
        s, a, b = out['alt_graph']
        out['ref_graph'] = (s[:a] + s[b:], a, a)
    ref = groups.get('REF')
    out['ref_read'] = None
    if ref is not None and len(ref) >= 2:
        rcov = x[4, ref] > 0
        rread = _consensus(x[0, ref], rcov)
        if all(rread[c] is not None for c in range(s0, s1)):
            out['ref_read'] = assemble(rread, rread)
        elif s1 == s0:
            left = _spell(rread, s0 - 1, -1, -1, flank)
            right = _spell(rread, s0, width, 1, flank)
            out['ref_read'] = (left + right, len(left), len(left))
    return out


class TensorSet:
    """Candidate id -> tensor + variant_summary record, for one tensor set (one kind)."""

    def __init__(self, tensor_dir, kind):
        self.dir = f'{tensor_dir.rstrip("/")}/{kind}'
        self.shards = {}

    def locations(self, chrom, wanted):
        """{candidate_id: dict(shard_file, index_within_shard)} from the (small) labels ndjson."""
        found = {}
        rx = re.compile(r'"candidate_id": ?"([^"]+)".*?"shard_file": ?"([^"]+)", ?"index_within_shard": ?(\d+)')
        with open(f'{self.dir}/{chrom}_labels.ndjson') as f:
            for line in f:
                m = rx.search(line)
                if m and m.group(1) in wanted:
                    found[m.group(1)] = dict(shard_file=m.group(2), index_within_shard=int(m.group(3)))
        return found

    def summaries(self, chrom, wanted):
        """variant_summary records of the wanted candidate ids (A1 candidate_id) on chrom."""
        found = {}
        cid_re = re.compile(r'"candidate_id": ?"([^"]+)"')
        with open(f'{self.dir}/{chrom}_variant_summary.ndjson') as f:
            for line in f:
                m = cid_re.search(line)
                if m and m.group(1) in wanted:
                    found[m.group(1)] = json.loads(line)
        return found

    def tensor(self, rec):
        sh = rec['shard_file']
        if sh not in self.shards:
            self.shards[sh] = np.load(f'{self.dir}/{sh}', mmap_mode='r')
        i = rec['index_within_shard']
        return read_rows(self.shards[sh], i, i + 1)[0]


def query_at(ctx, f):
    """Centered query: f bases each side of the site interval (fewer when the context is shorter)."""
    s, a, b = ctx
    return s[max(0, a - f):b + f], min(a, f), min(len(s) - b, f)


def adaptive_query(alt, ref_strings, f0=10, fmax=40, step=5):
    """Query with flank f0, widened until it is not a substring of any REF context (repeat), up to fmax."""
    f = f0
    while True:
        q, lf, rf = query_at(alt, f)
        amb = any(r and (q in r or rc(q) in r) for r in ref_strings)
        if not amb or f >= fmax or (lf < f and rf < f):
            return q, f, amb
        f += step


def homopolymer_at(s, a, b):
    """Longest homopolymer run touching the site interval [a, b) of s (b == a: the deletion boundary)."""
    best = 0
    for p in set([max(0, a - 1), min(len(s) - 1, a), min(len(s) - 1, max(a, b - 1)), min(len(s) - 1, b)]):
        if not s:
            break
        ch = s[p]
        lo = p
        while lo > 0 and s[lo - 1] == ch:
            lo -= 1
        hi = p
        while hi + 1 < len(s) and s[hi + 1] == ch:
            hi += 1
        best = max(best, hi - lo + 1)
    return best


def str_at(s, a, b, max_unit=6, min_copies=3):
    """(unit, copies) of the longest tandem repeat (unit 2..6, >= 3 copies) overlapping [a - 1, b + 1) of s."""
    best = ('', 0)
    for u in range(2, max_unit + 1):
        for start in range(max(0, a - u * 12), min(len(s), b + 1)):
            unit = s[start:start + u]
            if len(unit) < u:
                break
            k = 1
            while s[start + k * u:start + (k + 1) * u] == unit:
                k += 1
            end = start + k * u
            if k >= min_copies and start <= b and end >= a - 1 and k * u > best[1] * len(best[0] or 'x'):
                best = (unit, k)
    return best


# ---------------------------------------------------------------- graph

def read_rows(mm, lo, hi):
    """Rows [lo, hi) of a memory-mapped .npy (C order) by one sequential read: page faults of memory-mapped files on
    the shared file system are far slower than large reads."""
    row = int(np.prod(mm.shape[1:], dtype=np.int64)) * mm.dtype.itemsize
    with open(mm.filename, 'rb') as f:
        f.seek(mm.offset + lo * row)
        buf = f.read((hi - lo) * row)
    return np.frombuffer(buf, dtype=mm.dtype).reshape((hi - lo,) + mm.shape[1:])


class _Win:
    """An array with the rows [lo, hi) in RAM and the rest memory-mapped."""

    def __init__(self, mm, lo, hi):
        self.mm, self.lo, self.hi = mm, lo, hi
        self.ram = read_rows(mm, lo, hi)

    def __getitem__(self, i):
        return self.ram[i - self.lo] if self.lo <= i < self.hi else self.mm[i]

    def __len__(self):
        return len(self.mm)


def node_range(graph, chrom):
    tsv = GRAPHS[graph]['refpath'].replace('.grch38_path', '.chr_node_ranges.tsv')
    for line in open(tsv):
        f = line.rstrip('\n').split('\t')
        if f[0] == chrom:
            return int(f[1]), int(f[2])
    return None


class Graph:
    def __init__(self, name, in_memory=False):
        """Edge and reference-path arrays memory-mapped; focus(chrom) reads one chromosome component's node range
        into RAM (random page faults on the shared file system are slow); in_memory reads everything (d9 2 GB, d46 6 GB)."""
        g = GRAPHS[name]
        self.name, self.H = name, g['H']
        mode = None if in_memory else 'r'
        self.off = np.load(f'{g["edges"]}/offsets.npy', mmap_mode=mode)
        self.tgt = np.load(f'{g["edges"]}/targets.npy', mmap_mode=mode)
        self.rp = ReferencePath(g['refpath'])
        if in_memory:
            for k in ('lengths', 'visits', 'chrom', 'start0'):
                setattr(self.rp, k, np.load(f'{g["refpath"]}/{k}.npy'))
        self.db = sqlite3.connect(f'file:{g["index"]}?mode=ro', uri=True)
        self._seq = {}
        self._pc = {}
        self._hwin = None

    def focus(self, chrom):
        r = node_range(self.name, chrom)
        if r is None:
            return 0
        lo, hi = r[0], r[1] + 1
        for k in ('lengths', 'visits', 'chrom', 'start0'):
            setattr(self.rp, k, _Win(getattr(self.rp, k), lo, hi))
        off = read_rows(self.off, 2 * lo, 2 * hi + 1)
        self._hwin = (2 * lo, 2 * hi, off, read_rows(self.tgt, int(off[0]), int(off[-1])))
        return hi - lo

    def succ(self, h):
        w = self._hwin
        if w is not None and w[0] <= h < w[1]:
            base = int(w[2][0])
            return w[3][int(w[2][h - w[0]]) - base:int(w[2][h - w[0] + 1]) - base].tolist()
        return self.tgt[self.off[h]:self.off[h + 1]].tolist()

    def length(self, node):
        return int(self.rp.lengths[node])

    def on_ref(self, node):
        return 0 < node < len(self.rp.visits) and int(self.rp.visits[node]) >= 1

    def ref_coord(self, node):
        """(contig, start0, end0) of a node with exactly one GRCh38 visit, else None."""
        if not (0 < node < len(self.rp.visits)) or int(self.rp.visits[node]) != 1:
            return None
        s = int(self.rp.start0[node])
        return self.rp.contigs[int(self.rp.chrom[node])], s, s + self.length(node)

    def preload(self, chrom):
        """Sequences and path counts of a chromosome component's node-ID range: from a compressed cache
        (DATA/cache/<graph>_<chrom>_nodes.npz, written on first use) or one sequential sqlite scan."""
        r = node_range(self.name, chrom)
        if r is None:
            return 0
        lo, hi = r
        cache = f'{DATA}/cache/{self.name}_{chrom}_nodes.npz'
        if os.path.exists(cache):
            z = np.load(cache)
            ids, pcs, off, buf = z['ids'], z['pcs'], z['off'], z['seq'].tobytes().decode()
        else:
            rows = self.db.execute('select node_id, seq, distinct_path_count from nodes where node_id between ? and ? '
                                   'order by node_id', (lo, hi)).fetchall()
            ids = np.array([x[0] for x in rows], np.int64)
            pcs = np.array([x[2] for x in rows], np.int64)
            seqs = [x[1].upper() for x in rows]
            off = np.zeros(len(seqs) + 1, np.int64)
            off[1:] = np.cumsum([len(x) for x in seqs])
            buf = ''.join(seqs)
            os.makedirs(os.path.dirname(cache), exist_ok=True)
            tmp = f'{cache}.{os.getpid()}.tmp.npz'
            np.savez_compressed(tmp, ids=ids, pcs=pcs, off=off, seq=np.frombuffer(buf.encode(), np.uint8))
            os.replace(tmp, cache)
        ids_l, pcs_l, off_l = ids.tolist(), pcs.tolist(), off.tolist()
        for k, n in enumerate(ids_l):
            self._seq[n] = buf[off_l[k]:off_l[k + 1]]
            self._pc[n] = pcs_l[k]
        return len(ids_l)

    def seqs(self, nodes):
        need = [n for n in set(nodes) if n not in self._seq]
        for k in range(0, len(need), 900):
            part = need[k:k + 900]
            q = f'select node_id, seq, distinct_path_count from nodes where node_id in ({",".join("?" * len(part))})'
            for n, s, pc in self.db.execute(q, part):
                self._seq[n] = s.upper()
                self._pc[n] = pc
        return {n: self._seq.get(n, '') for n in nodes}

    def path_count(self, node):
        self.seqs([node])
        return self._pc.get(node)

    def hseq(self, h):
        s = self.seqs([h >> 1])[h >> 1]
        return rc(s) if h & 1 else s

    def nearest_ref(self, h, max_bp=200_000, max_nodes=50_000):
        """Walk forward from handle h (Dijkstra by bp, not through GRCh38 nodes): the first GRCh38 node(s) reached,
        [(bp from the end of h, node)] at the smallest distance; [] if none within max_bp / max_nodes."""
        start = h >> 1
        dist = {h: 0}
        pq = [(0, h)]
        found, best = [], None
        while pq:
            d, x = heapq.heappop(pq)
            if d > dist.get(x, 1 << 62) or (best is not None and d > best):
                continue
            for y in self.succ(x):
                node = y >> 1
                if node == start:
                    continue
                if self.on_ref(node):
                    if best is None or d <= best:
                        best = d
                        found.append((d, node))
                    continue
                nd = d + self.length(node)
                if nd <= max_bp and nd < dist.get(y, 1 << 62):
                    dist[y] = nd
                    heapq.heappush(pq, (nd, y))
            if len(dist) > max_nodes:
                break
        return [(d, n) for d, n in found if d == best]

    def locus(self, node, chrom):
        """GRCh38 interval of the bubble around an off-reference node: its nearest GRCh38 node on each side.
        Returns dict(status, chrom, lo, hi, left_bp, right_bp)."""
        if self.on_ref(node):
            c = self.ref_coord(node)
            if c and c[0] == chrom:
                return dict(status='on_ref', chrom=chrom, lo=c[1], hi=c[2], left_bp=0, right_bp=0)
        sides = {}
        for side, h in (('right', 2 * node), ('left', 2 * node + 1)):
            hits = [(d, self.ref_coord(n)) for d, n in self.nearest_ref(h)]
            hits = [(d, c) for d, c in hits if c and c[0] == chrom]
            if hits:
                sides[side] = hits
        if not sides:
            return dict(status='no_anchor', chrom=chrom, lo=None, hi=None, left_bp=None, right_bp=None)
        coords = [c for hs in sides.values() for _, c in hs]
        lo, hi = min(c[1] for c in coords), max(c[2] for c in coords)
        status = 'bubble' if len(sides) == 2 else 'one_side'
        if hi - lo > 100_000:
            status = 'wide'
        return dict(status=status, chrom=chrom, lo=lo, hi=hi,
                    left_bp=sides['left'][0][0] if 'left' in sides else None,
                    right_bp=sides['right'][0][0] if 'right' in sides else None)

    def neighborhood(self, node, max_bp, max_nodes=20_000):
        """Handles within max_bp of the node (both directions, through any node), as a set of handles."""
        dist = {}
        pq = []
        for h in (2 * node, 2 * node + 1):
            dist[h] = 0
            heapq.heappush(pq, (0, h))
        while pq and len(dist) < max_nodes:
            d, x = heapq.heappop(pq)
            if d > dist.get(x, 1 << 62):
                continue
            for y in self.succ(x):
                nd = d + self.length(y >> 1)
                if nd <= max_bp and nd < dist.get(y, 1 << 62):
                    dist[y] = nd
                    heapq.heappush(pq, (nd, y))
        hs = set(dist)
        hs |= {h ^ 1 for h in hs}
        return hs

    def spell(self, q, handles, max_steps=200_000):
        """A walk inside `handles` spelling q (any orientation, since both handles of every node are in the set):
        returns the list of nodes of the first walk found, None if there is none, 'capped' past max_steps."""
        if not q:  # an empty query (a deletion context without flank) spells nothing (IndexError in jobs 387265_8/_10)
            return None
        nodes = {h >> 1 for h in handles}
        self.seqs(list(nodes))
        first = q[0]
        steps = [0]

        def extend(h, i, j, walk):
            steps[0] += 1
            if steps[0] > max_steps:
                raise TimeoutError
            s = self.hseq(h)
            n = min(len(s) - i, len(q) - j)
            if s[i:i + n] != q[j:j + n]:
                return None
            walk = walk + [h >> 1]
            if j + n == len(q):
                return walk
            for y in self.succ(h):
                if y in handles:
                    r = extend(y, 0, j + n, walk)
                    if r:
                        return r
            return None

        try:
            for h in handles:
                s = self.hseq(h)
                i = s.find(first)
                while i >= 0:
                    r = extend(h, i, 0, [])
                    if r:
                        return r
                    i = s.find(first, i + 1)
        except TimeoutError:
            return 'capped'
        return None


# ---------------------------------------------------------------- HPRC v1.1 haplotypes at a GRCh38 locus

class HprcHaplotypes:
    """Haplotype sequences of the full HPRC v1.1 graph (vg deconstruct raw VCF, top-level sites, 44 diploid
    samples + CHM13 = 89 haplotypes, every allele kept regardless of frequency) over a GRCh38 window."""

    def __init__(self, vcf=HPRC_VCF, fasta=FASTA):
        import pysam
        self.vcf = pysam.TabixFile(vcf)
        self.fa = pysam.FastaFile(fasta)
        header = [h for h in self.vcf.header if h.startswith('#CHROM')][0].split('\t')
        self.samples = header[9:]

    def preload(self, chrom):
        """Every top-level record of a chromosome in memory (sorted starts for bisect), and the chromosome's GRCh38
        sequence: from a compressed cache (DATA/cache/hprc_v11_<chrom>.pkl.gz, written on first use) or one
        sequential read of the VCF and the FASTA."""
        import gzip as _gz
        import pickle
        cache = f'{DATA}/cache/hprc_v11_{chrom}.pkl.gz'
        if os.path.exists(cache):
            with _gz.open(cache, 'rb') as f:
                recs, seq = pickle.load(f)
        else:
            recs = []
            for line in self.vcf.fetch(chrom):
                f = line.rstrip('\n').split('\t')
                if 'LV=0' not in f[7].split(';'):
                    continue
                pos0 = int(f[1]) - 1
                recs.append((pos0, pos0 + len(f[3]), f[3], f[4], '\t'.join(g.split(':')[0] for g in f[9:])))
            seq = self.fa.fetch(chrom).upper()
            os.makedirs(os.path.dirname(cache), exist_ok=True)
            tmp = f'{cache}.{os.getpid()}.tmp'
            with _gz.open(tmp, 'wb', compresslevel=3) as f:
                pickle.dump((recs, seq), f, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(tmp, cache)
        recs.sort()
        self.mem = getattr(self, 'mem', {})
        from itertools import accumulate
        starts = [r[0] for r in recs]
        ends_max = list(accumulate((r[1] for r in recs), max))   # non-decreasing: bisect for the first end > lo
        self.mem[chrom] = (starts, recs, ends_max)
        self.seq_mem = getattr(self, 'seq_mem', {})
        self.seq_mem[chrom] = seq
        return len(recs)

    def _records(self, chrom, lo, hi):
        mem = getattr(self, 'mem', {}).get(chrom)
        if mem is None:
            for line in self.vcf.fetch(chrom, max(0, lo), hi):
                f = line.rstrip('\n').split('\t')
                if 'LV=0' not in f[7].split(';'):
                    continue
                pos0 = int(f[1]) - 1
                yield pos0, pos0 + len(f[3]), f[3], f[4], f[9:]
            return
        import bisect
        starts, recs, ends_max = mem
        k = bisect.bisect_right(ends_max, lo - 1)
        while k < len(recs) and recs[k][0] < hi:
            r = recs[k]
            if r[1] > lo or (r[1] == r[0] and r[0] >= lo):
                yield r[0], r[1], r[2], r[3], r[4].split('\t')
            k += 1

    def window(self, chrom, lo, hi):
        """dict(hap name -> sequence) over [lo', hi') >= [lo, hi) extended to cover every overlapping top-level
        site; 'GRCh38' is the reference itself. Also the number of top-level sites and of overlapping skips."""
        recs = [(s0, e0, r.upper(), [x.upper() for x in a.split(',')], gts)
                for s0, e0, r, a, gts in self._records(chrom, lo, hi)]
        if recs:
            lo = min(lo, min(r[0] for r in recs))
            hi = max(hi, max(r[1] for r in recs))
        lo = max(0, lo)
        mem = getattr(self, 'seq_mem', {}).get(chrom)
        ref = mem[lo:hi] if mem is not None else self.fa.fetch(chrom, lo, hi).upper()
        haps = {'GRCh38': ref}
        recs.sort()
        skipped = 0
        for si, name in enumerate(self.samples):
            for hk in (0, 1):
                parts, cur, ok = [], lo, True
                for s0, e0, r, alts, gts in recs:
                    gt = gts[si].split(':')[0].replace('/', '|').split('|')
                    if hk >= len(gt):
                        ok = False
                        break
                    a = gt[hk]
                    if a in ('.', '0'):
                        continue
                    if s0 < cur:
                        skipped += 1
                        continue
                    alt = alts[int(a) - 1]
                    if alt.startswith('<') or alt == '*':
                        continue
                    parts.append(ref[cur - lo:s0 - lo])
                    parts.append(alt)
                    cur = e0
                if not ok:
                    continue
                parts.append(ref[cur - lo:])
                haps[f'{name}#{hk + 1}'] = ''.join(parts)
        return haps, lo, hi, len(recs), skipped

    @staticmethod
    def carriers(haps, q):
        r = rc(q)
        return [h for h, s in haps.items() if q in s or r in s]

"""How the d9 graph spells each I1 (allele fully in the graph) somatic INDEL truth of HG008T Illumina.

For each truth: enumerate the local graph paths between the GRCh38 nodes ~60 bp around it, keep those whose sequence is
the truth ALT haplotype (GRCh38 with the VCF allele applied), and classify how the best such path leaves GRCh38:
  INS single bubble node  - one branch node V between adjacent GRCh38 nodes X, Y (X->V->Y), seq(V) = the insertion
  INS multi-node branch    - a pure insertion spelled by >= 2 consecutive branch nodes between adjacent X, Y
  branch replaces GRCh38   - a branch run between non-adjacent GRCh38 nodes (it replaces some reference bases)
  DEL skip edge            - only GRCh38 nodes, one edge skips the deleted bases
  several events            - more than one departure from GRCh38
  matches with nearby germline - only GRCh38 + the truth + nearby germline (PASS, up to 6) is a graph path
  closest graph allele differs by 1 bp / >= 2 bp - no path; edit distance of the closest graph haplotype
  no matching path (cap) / no GRCh38 node in the window
Input: analysis/hg008_somatic_miss_20260926/illumina/indel_no_candidate.tsv (final_class I1*). Output: classes.tsv.
"""
import csv, itertools, json, sqlite3, sys
from collections import Counter
import numpy as np
import pysam
IDX = '/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13.path_index.sqlite'
FASTA = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
MISS = '/scratch/jshen/Github/Pansoma/analysis/hg008_somatic_miss_20260926/illumina/indel_no_candidate.tsv'
OUT = '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930'
FLANK, CAP, MAX_NODES = 60, 20000, 400
RC = str.maketrans('ACGTN', 'TGCAN'); FLIP = {'+': '-', '-': '+'}
db = sqlite3.connect(f'file:{IDX}?mode=ro', uri=True)   # node sequences and edges (complete); its path_coords
# lack GRCh38 for chr3-9 and part of chr22, so GRCh38 coordinates come from the pipeline's reference-path index
GP = '/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path'
GMETA = json.load(open(f'{GP}/meta.json'))
NAMES = [c['name'] for c in GMETA['contigs']]
G_CHROM, G_START, G_LEN, G_REV, G_VISITS = (np.load(f'{GP}/{f}.npy', mmap_mode='r') for f in ('chrom', 'start0', 'lengths', 'reverse', 'visits'))
P_NODES, P_STARTS = np.load(f'{GP}/path_nodes.npy', mmap_mode='r'), np.load(f'{GP}/path_starts.npy', mmap_mode='r')
_seq, _succ, _ref = {}, {}, {}

def seq(n):
    if n not in _seq:
        _seq[n] = db.execute('select seq from nodes where node_id=?', (str(n),)).fetchone()[0]
    return _seq[n]

def oseq(n, o):
    return seq(n) if o == '+' else seq(n).translate(RC)[::-1]

def succ(n, o):
    if (n, o) not in _succ:
        out = {(int(m), mo) for m, mo in db.execute('select to_id, to_orient from edges where from_id=? and from_orient=?', (str(n), o))}
        out |= {(int(m), FLIP[mo]) for m, mo in db.execute('select from_id, from_orient from edges where to_id=? and to_orient=?', (str(n), FLIP[o]))}
        _succ[(n, o)] = sorted(out)
    return _succ[(n, o)]

def ref(n):
    """(chrom, start0, end0) if n is a GRCh38 node (visited once, forward), else None."""
    n = int(n)
    if n not in _ref:
        c = int(G_CHROM[n])
        _ref[n] = (NAMES[c], int(G_START[n]), int(G_START[n]) + int(G_LEN[n]) - 1) if c >= 0 and G_VISITS[n] == 1 else None
    return _ref[n]


def window_nodes(chrom, lo, hi):
    """GRCh38 nodes (node, start0, end0) overlapping [lo, hi], in path order."""
    co = GMETA['contigs'][NAMES.index(chrom)]
    a, b = co['offset'], co['offset'] + co['nodes']
    starts = P_STARTS[a:b]
    i = max(0, int(np.searchsorted(starts, lo, side='right')) - 1)
    j = int(np.searchsorted(starts, hi, side='right'))
    out = []
    for k in range(a + i, a + j):
        n = int(P_NODES[k])
        if G_REV[n]:
            return None                                  # GRCh38 traverses a node in reverse: not handled
        out.append((n, int(P_STARTS[k]), int(P_STARTS[k]) + int(G_LEN[n]) - 1))
    return out

def events(path):
    """Departures from GRCh38 along a path: ('skip', bases) or ('branch', n_nodes, replaced_bases, branch_seq)."""
    out, last, run = [], None, []
    for n, o in path:
        r = ref(n)
        if r is None:
            run.append((n, o)); continue
        if last is not None:
            gap = r[1] - last[2] - 1
            if run:
                out.append(('branch', len(run), gap, ''.join(oseq(*x) for x in run)))
            elif gap != 0:
                out.append(('skip', gap))
        last, run = r, []
    return out

def graph_haplotypes(chrom, first, last, targets=None, limit=3000, max_len=None):
    """Sequences spelled from the first to the last GRCh38 node (all, or only those in targets)."""
    out, stack = set(), [((first[0], '+'), oseq(first[0], '+'))]
    while stack and len(out) < limit:
        (n, o), s = stack.pop()
        if n == last[0]:
            if targets is None or s in targets:
                out.add(s)
            continue
        if targets is not None and not any(t.startswith(s) for t in targets):
            continue
        if max_len and len(s) > max_len:
            continue
        for m, mo in succ(n, o):
            r = ref(m)
            if r and (r[0] != chrom or r[1] < first[1] or r[1] > last[2]):
                continue
            stack.append(((m, mo), s + oseq(m, mo)))
    return out


def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def fallback(chrom, first, last, region, off, vref, valt, alt_hap):
    """No path spells GRCh38 + the truth: (1) GRCh38 + the truth + up to 6 nearby germline (PASS) alleles;
    (2) else the edit distance of the closest graph haplotype to GRCh38 + the truth."""
    base = first[1]
    som = (off, vref.upper(), valt.upper(), None)
    g = [(p - 1 - base, a, b, f'{p}:{a}>{b}') for p, a, alts in GERMLINE.get(chrom, ())
         if base <= p - 1 and p - 1 + len(a) <= last[2] + 1 for b in alts][:6]
    targets = {}
    for k in range(1, len(g) + 1):
        for combo in itertools.combinations(g, k):
            vs = sorted(list(combo) + [som], key=lambda v: v[0])
            if any(vs[i][0] + len(vs[i][1]) > vs[i + 1][0] for i in range(len(vs) - 1)):
                continue
            parts, cur = [], 0
            for o, a, b, _ in vs:
                parts += [region[cur:o], b]; cur = o + len(a)
            parts.append(region[cur:])
            targets.setdefault(''.join(parts), [v[3] for v in vs if v[3]])
    if targets:
        hits = graph_haplotypes(chrom, first, last, set(targets), limit=1)
        if hits:
            return 'matches with nearby germline', ';'.join(targets[next(iter(hits))])
    haps = graph_haplotypes(chrom, first, last, limit=3000, max_len=len(alt_hap) + 200)
    if not haps:
        return 'no matching path', ''
    d = min(lev(h, alt_hap) for h in haps)
    return ('closest graph allele differs by 1 bp' if d == 1 else 'closest graph allele differs by >= 2 bp'), f'edit distance {d}'


def classify(chrom, pos1, vref, valt):
    lo, hi = pos1 - 1 - FLANK, pos1 - 1 + len(vref) + FLANK
    nodes = window_nodes(chrom, lo, hi)
    if nodes is None:
        return 'GRCh38 node in reverse (not handled)', ''
    if not nodes:
        return 'no GRCh38 node in the window', ''
    first, last = nodes[0], nodes[-1]
    region = fa.fetch(chrom, first[1], last[2] + 1).upper()
    off = pos1 - 1 - first[1]
    if off < 0 or region[off:off + len(vref)] != vref.upper():
        return 'error: REF mismatch', ''
    alt_hap = region[:off] + valt.upper() + region[off + len(vref):]
    found, n_paths, stack = [], 0, [((first[0], '+'), oseq(first[0], '+'), [(first[0], '+')])]
    while stack and n_paths < CAP:
        (n, o), s, path = stack.pop()
        if not alt_hap.startswith(s) and n != last[0]:
            continue                                  # prune: the ALT haplotype must extend this prefix
        if n == last[0]:
            n_paths += 1
            if s == alt_hap:
                found.append(path)
            continue
        if len(path) > MAX_NODES:
            continue
        for m, mo in succ(n, o):
            r = ref(m)
            if r and (r[0] != chrom or r[1] < first[1] or r[1] > last[2]):
                continue
            stack.append(((m, mo), s + oseq(m, mo), path + [(m, mo)]))
    if not found:
        if n_paths >= CAP:
            return 'no matching path (cap)', ''
        return fallback(chrom, first, last, region, off, vref, valt, alt_hap)
    best = None
    for p in found:
        ev = events(p)
        if len(ev) == 1 and ev[0][0] == 'branch' and ev[0][2] == 0:
            c = 'INS single bubble node' if ev[0][1] == 1 else 'INS multi-node branch'
        elif len(ev) == 1 and ev[0][0] == 'branch':
            c = 'branch replaces GRCh38'
        elif len(ev) == 1 and ev[0][0] == 'skip':
            c = 'DEL skip edge'
        else:
            c = 'several events'
        rank = ['INS single bubble node', 'DEL skip edge', 'INS multi-node branch', 'branch replaces GRCh38', 'several events'].index(c)
        if best is None or rank < best[0]:
            best = (rank, c, ev)
    return best[1], ';'.join(':'.join(map(str, e)) for e in best[2])

fa = pysam.FastaFile(FASTA)
chroms = set(sys.argv[1:])
GERMLINE = {}   # chrom -> [(pos, ref, [alts])], PASS records of the labels' germline truth (no index: one pass)
_man = json.load(open('/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors/INDEL/labels.manifest.json'))
for _rec in pysam.VariantFile(_man['truth']['germline']['vcf']):
    if (chroms and _rec.chrom not in chroms) or (_rec.filter.keys() and 'PASS' not in _rec.filter.keys()):
        continue
    GERMLINE.setdefault(_rec.chrom, []).append((_rec.pos, _rec.ref.upper(), [a.upper() for a in _rec.alts or () if a and a[0] != '<']))
rows = [r for r in csv.DictReader(open(MISS), delimiter='\t') if r['final_class'].startswith('I1') and (not chroms or r['chrom'] in chroms)]
with open(f'{OUT}/classes{"_" + "_".join(sorted(chroms)) if chroms else ""}.tsv', 'w') as w:
    w.write('chrom\tvcf_pos\tvcf_ref\tvcf_alt\tkind\tlength\tin_bed\treads_reason\tgraph_class\tevents\n')
    for k, r in enumerate(rows):
        c, ev = classify(r['chrom'], int(r['vcf_pos']), r['vcf_ref'], r['vcf_alt'])
        w.write(f"{r['chrom']}\t{r['vcf_pos']}\t{r['vcf_ref']}\t{r['vcf_alt']}\t{r['kind']}\t{r['length']}\t{r['in_bed']}\t{r['reason']}\t{c}\t{ev}\n")
        if k % 100 == 0:
            print(k, len(rows), flush=True)

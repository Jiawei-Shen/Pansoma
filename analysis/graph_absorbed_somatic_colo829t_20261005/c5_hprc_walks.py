"""COLO829T step 5: every HPRC haplotype's local walk through each window of loci.tsv, from one parallel pass over the d9 GFA.

Port of analysis/graph_absorbed_somatic_20261001/s4_hprc_walks.py (HG008T), logic unchanged; only D (the COLO829T data
dir) and the validation seed differ.
Input: $D/loci.tsv (c1b; 497 truths, first_node / last_node anchors, win_start0 / win_end0), hprc-v1.1-mc-grch38.d9.gfa,
the pipeline grch38_path arrays (lengths.npy = every node's length, from the same GFA: meta.json source size is checked).
GFA layout (d9): S lines, then P lines = clipped HPRC haplotype fragments 'sample#hap#contig#offset' (segments
'12+,13-'), then W lines (GRCh38, CHM13 and unclipped HPRC contigs, 'W sample hap seqid start end >12<13'), then L lines.
Method (graph_prep.audit_part style): byte ranges aligned to line starts, one pool of processes; each W/P line ->
int64 ids (translate + np.fromstring); a boolean node mask marks every anchor and every interior GRCh38 node of a window;
only lines with a mark are decoded further (orientation from the '<' / '-' bytes, contig offsets from the lengths).
Per locus and line, the anchors are paired in line order: forward +first ... +last, reverse -last ... -first (reversed
and flipped to GRCh38-forward orientation). A repeated first anchor before the last one stays one traversal (a cycle in
the window) unless more than MAXN nodes apart; every pairing is a separate traversal (duplications).
  complete     both anchors, <= MAXN nodes apart; local_path = '>187380>187381<187385...' first..last
  too_long     both anchors in order but more than MAXN nodes apart (no path)
  left_only    first anchor only (the fragment ends, or leaves the window, before the last anchor)
  right_only   last anchor only
  inner_only   interior GRCh38 nodes of the window but neither anchor
For partial rows partial_end = line_end (the fragment ends <= PART_MAX nodes from the anchor without leaving the window:
the d9 frequency filter cut the walk), exits_window (it reaches a GRCh38 node outside the window first, e.g. a germline
deletion of the other anchor) or continues (neither within PART_MAX nodes); partial_path = the walk from the anchor up to
that point (GRCh38-forward orientation), empty for continues. Partial rows inside a complete traversal of the same line
(an anchor revisited in the other orientation inside the window) are dropped. inner_only: partial_path = the walk from
the first to the last interior node hit. A haplotype whose germline deletion removes an anchor node gets no complete
traversal (left_only / right_only exits_window; its partial_path still runs from the other anchor up to the deletion).
hap_lo0 / hap_hi0 = 0-based half-open span of the traversal on the haplotype contig (seq_start = W start, or the 4th
PanSN field of P names = the fragment's offset on its contig; assumption, checked only as non-overlap of fragments).
Outputs ($D): hprc_local_paths.tsv.gz (truth_id sample hap seqid seq_start traversal status strand hap_lo0 hap_hi0
n_nodes local_path partial_end partial_path; traversal numbered per truth_id x sample x hap), hprc_haplotypes.tsv
(every sample/hap: W and P lines, nodes, bases, contigs), hprc_window_coverage.tsv (per truth_id over the HPRC
haplotypes = every sample/hap but GRCh38 / CHM13: complete / partial only / absent, distinct paths, GRCh38 and CHM13
status; grch38_ok = one complete GRCh38 traversal equal to the all-GRCh38 window path at win_start0..win_end0).
Validation printed to the log: grch38_ok for every locus, and for 5 loci a few distinct haplotype walks spelled with
oseq(), edge-checked with succ() and compared with ref_hap / alt_hap (edit distance) and with hap_hi0 - hap_lo0.
Assumptions: HPRC v1.1 = 44 samples x 2 haplotypes (88; no EUR sample); a haplotype absent from a window here means its
d9 walk does not touch the window (clipped fragment, frequency-filtered), not that it lacks the sequence.
Run: sbatch tmp/graph_absorbed_somatic_colo829t_20261005/c5_hprc_walks.sbatch (24 processes, 16G = HG008 s4 MaxRSS
13.4 GB + 20 %).
"""
import csv, gzip, json, os, sys, time
from collections import Counter, defaultdict
from multiprocessing import get_context
import numpy as np
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005'
C = '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/classify.py'
GFA = '/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gfa'
MAXN, PART_MAX = 20000, 500
REFS = ('GRCh38', 'CHM13')
WSEP, PSEP = bytes.maketrans(b'<>', b'  '), bytes.maketrans(b',+-', b'   ')
ns = {}
exec(open(C).read().split("fa = pysam.FastaFile(FASTA)\nchroms")[0], ns)
GP, NAMES = ns['GP'], ns['NAMES']
LEN = np.load(f'{GP}/lengths.npy')                              # in memory, shared by the forked workers
assert json.load(open(f'{GP}/meta.json'))['source']['size'] == os.path.getsize(GFA)
G_CHROM, G_START, G_VISITS = ns['G_CHROM'], ns['G_START'], ns['G_VISITS']
loci = list(csv.DictReader(open(f'{D}/loci.tsv'), delimiter='\t'))
ROLES, REFPATH = defaultdict(list), {}
for i, r in enumerate(loci):
    nodes = ns['window_nodes'](r['chrom'], int(r['win_start0']), int(r['win_end0']))
    assert nodes[0][0] == int(r['first_node']) and nodes[-1][0] == int(r['last_node']), r['truth_id']
    REFPATH[i] = ''.join(f'>{n}' for n, _, _ in nodes)
    for k, (n, _, _) in enumerate(nodes):
        ROLES[n].append((i, 'F' if k == 0 else 'L' if k == len(nodes) - 1 else 'I'))
MARK = np.zeros(len(LEN), dtype=bool)
MARK[list(ROLES)] = True
WIN = [(NAMES.index(r['chrom']), int(r['win_start0']), int(r['win_end0'])) for r in loci]


def pansn(name):
    f = name.split('#')
    if len(f) >= 4 and f[3].isdigit():
        return f[0], f[1], f[2], int(f[3])
    if len(f) == 3:
        seqid, s0 = f[2], 0
        if seqid.endswith(']') and '[' in seqid:
            seqid, rng = seqid[:-1].split('[')
            s0 = int(rng.split('-')[0])
        return f[0], f[1], seqid, s0
    return name, '0', name, 0


def walk(ids, rev, lo, hi, flip):
    """Oriented walk of line positions lo..hi (inclusive), in GRCh38-forward orientation."""
    if hi < lo:
        return ''
    a, o = ids[lo:hi + 1], rev[lo:hi + 1]
    if flip:
        a, o = a[::-1], ~o[::-1]
    return ''.join(('<' if x else '>') + str(n) for n, x in zip(a.tolist(), o.tolist()))


def partial(ids, li, p, step):
    """From anchor position p in direction step (+1/-1): (partial_end, last position kept or None)."""
    c, lo, hi = WIN[li]
    room = ids.size - 1 - p if step == 1 else p                    # nodes left on the line beyond the anchor
    q = p + step * np.arange(1, min(room, PART_MAX) + 1)
    n = ids[q]
    out = (G_CHROM[n] >= 0) & (G_VISITS[n] == 1) & ((G_CHROM[n] != c) | (G_START[n] < lo) | (G_START[n] > hi))
    if out.any():
        k = int(np.argmax(out))
        return f'exits_window:{NAMES[int(G_CHROM[n[k]])]}:{int(G_START[n[k]])}', (int(q[k - 1]) if k else p)
    if room <= PART_MAX:
        return 'line_end', (int(q[-1]) if q.size else p)
    return 'continues', None


def traversals(li, ev, ids, rev):
    """ev: [(pos, role)] of locus li in line order -> [(status, lo, hi, strand, partial_end, path)]."""
    found = []
    for strand, seq in (('+', [(p, r) for p, r in ev if r != 'I' and not rev[p]]),
                        ('-', [(p, r) for p, r in reversed(ev) if r != 'I' and rev[p]])):
        op = None
        for p, r in seq:
            if r == 'F':
                if op is not None and abs(p - op) > MAXN:
                    found.append(('left_only', op, strand)); op = None
                op = p if op is None else op
            elif op is None:
                found.append(('right_only', p, strand))
            else:
                found.append(('complete' if abs(p - op) <= MAXN else 'too_long', (op, p), strand)); op = None
        if op is not None:
            found.append(('left_only', op, strand))
    spans = [(min(x[1]), max(x[1])) for x in found if x[0] in ('complete', 'too_long')]
    out = []
    for st, p, strand in found:
        flip = strand == '-'
        if st in ('complete', 'too_long'):
            lo, hi = min(p), max(p)
            out.append((st, lo, hi, strand, '', walk(ids, rev, lo, hi, flip) if st == 'complete' else ''))
            continue
        if any(a <= p <= b for a, b in spans):
            continue
        step = 1 if (st == 'left_only') != flip else -1          # direction away from the anchor along the line
        end, q = partial(ids, li, p, step)
        lo, hi = (p, p if q is None else q) if step == 1 else (p if q is None else q, p)
        out.append((st, lo, hi, strand, end, '' if q is None else walk(ids, rev, lo, hi, flip)))
    if not out and not found:                                       # interior GRCh38 nodes only
        lo, hi = ev[0][0], ev[-1][0]
        flip = bool(rev[lo])
        out.append(('inner_only', lo, hi, '-' if flip else '+', '',
                    walk(ids, rev, lo, hi, flip) if hi - lo <= PART_MAX else ''))
    return out


def part(task):
    start, end = task
    rows, haps, bad, ext = [], defaultdict(lambda: [0, 0, 0, 0, set()]), [], []
    with open(GFA, 'rb', buffering=1 << 24) as f:
        f.seek(max(start - 1, 0))
        if start:
            f.readline()
        pos = f.tell()
        while pos < end:
            line = f.readline()
            if not line:
                break
            pos += len(line)
            k = line[:2]
            if k == b'W\t':
                fl = line.rstrip(b'\n').split(b'\t')
                sample, hap, seqid, s0, w = fl[1].decode(), fl[2].decode(), fl[3].decode(), int(fl[4]), fl[6]
                ids = np.fromstring(w.translate(WSEP), dtype=np.int64, sep=' ')
                rc, fc, kind = 60, 62, 0
            elif k == b'P\t':
                fl = line.split(b'\t', 3)
                (sample, hap, seqid, s0), w = pansn(fl[1].decode()), fl[2].rstrip(b'\n')
                ids = np.fromstring(w.translate(PSEP), dtype=np.int64, sep=' ')
                rc, fc, kind = 45, 43, 1
            else:
                continue
            lens = LEN[ids].astype(np.int64)
            h = haps[(sample, hap)]
            h[kind] += 1; h[2] += int(ids.size); h[3] += int(lens.sum()); h[4].add(seqid)
            if kind == 0 and int(fl[5]) - s0 != int(lens.sum()):
                bad.append(('W_length', sample, hap, seqid, s0))
            hit = np.flatnonzero(MARK[ids])
            if not hit.size:
                continue
            ext.append((sample, hap, seqid, s0, int(lens.sum())))
            raw = np.frombuffer(w, dtype=np.uint8)
            rev = raw[(raw == rc) | (raw == fc)] == rc
            assert rev.size == ids.size, (sample, hap, seqid)
            cum = s0 + np.cumsum(lens) - lens
            ev = defaultdict(list)
            for p in hit.tolist():
                for li, role in ROLES[int(ids[p])]:
                    ev[li].append((p, role))
            for li, e in ev.items():
                for st, lo, hi, strand, pend, path in traversals(li, e, ids, rev):
                    rows.append((li, sample, hap, seqid, s0, st, strand, int(cum[lo]), int(cum[hi] + lens[hi]),
                                 hi - lo + 1, path if st == 'complete' else '', pend, '' if st == 'complete' else path))
    return rows, {k: v for k, v in haps.items()}, bad, ext


def main():
    t0 = time.time()
    procs = int(os.environ.get('SLURM_CPUS_PER_TASK') or 4)
    total = os.path.getsize(GFA)
    cuts = [total * k // (4 * procs) for k in range(4 * procs + 1)]
    rows, haps, bad, ext = [], {}, [], []
    with get_context('fork').Pool(procs) as pool:
        for r, h, b, e in pool.imap_unordered(part, list(zip(cuts, cuts[1:]))):
            rows += r; bad += b; ext += e
            for key, v in h.items():
                a = haps.setdefault(key, [0, 0, 0, 0, set()])
                for j in range(4):
                    a[j] += v[j]
                a[4] |= v[4]
    print(f'scan {time.time() - t0:.0f}s, {len(rows)} rows, {len(haps)} haplotypes, W length mismatches {len(bad)}', bad[:5])
    rows.sort(key=lambda x: (x[0], x[1], x[2], x[3], x[4], x[7]))
    trav = Counter()
    with gzip.open(f'{D}/hprc_local_paths.tsv.gz', 'wt') as w:
        w.write('truth_id\tsample\thap\tseqid\tseq_start\ttraversal\tstatus\tstrand\thap_lo0\thap_hi0\tn_nodes\t'
                'local_path\tpartial_end\tpartial_path\n')
        for x in rows:
            k = (x[0], x[1], x[2]); t = trav[k]; trav[k] += 1
            w.write('\t'.join(map(str, (loci[x[0]]['truth_id'],) + x[1:5] + (t,) + x[5:])) + '\n')
    with open(f'{D}/hprc_haplotypes.tsv', 'w') as w:
        w.write('sample\thap\treference\tn_W_lines\tn_P_lines\tn_lines\tnodes\tbases\tn_contigs\n')
        for (s, hp), v in sorted(haps.items()):
            w.write(f'{s}\t{hp}\t{s in REFS}\t{v[0]}\t{v[1]}\t{v[0] + v[1]}\t{v[2]}\t{v[3]}\t{len(v[4])}\n')
    hprc = sorted(k for k in haps if k[0] not in REFS)
    frag, n_over = defaultdict(list), 0                    # offsets: hit lines of one contig must not overlap
    for e in ext:
        frag[e[:3]].append(e[3:])
    for v in frag.values():
        v.sort()
        n_over += sum(v[k + 1][0] < v[k][0] + v[k][1] for k in range(len(v) - 1))
    print(f'lines with a window node {len(ext)}, contigs with >1 such line {sum(len(v) > 1 for v in frag.values())}, '
          f'overlapping offsets {n_over}')
    by = defaultdict(list)
    for x in rows:
        by[x[0]].append(x)
    cov_cols = ['truth_id', 'chrom', 'vcf_pos', 'kind2', 'n_hprc_haps', 'n_complete', 'n_partial_only', 'n_absent',
                'n_multi_complete', 'n_ref_path', 'n_nonref_path', 'n_distinct_paths', 'n_line_end', 'n_exits_window',
                'grch38_status', 'grch38_ok', 'chm13_status', 'chm13_ref_path']
    n_ok, cov = 0, []
    for i, r in enumerate(loci):
        st = defaultdict(list)
        for x in by[i]:
            st[(x[1], x[2])].append(x)
        comp = {k: [x[10] for x in v if x[5] == 'complete'] for k, v in st.items()}
        g = st.get(('GRCh38', '0'), [])
        gok = (len(g) == 1 and g[0][5] == 'complete' and g[0][10] == REFPATH[i] and g[0][7] == int(r['win_start0'])
               and g[0][8] == int(r['win_end0']) + 1 and g[0][3] == r['chrom'])
        n_ok += gok
        h = [k for k in hprc if k in st]
        c = [k for k in h if comp[k]]
        ch = st.get(('CHM13', '0'), [])
        cov.append(dict(truth_id=r['truth_id'], chrom=r['chrom'], vcf_pos=r['vcf_pos'], kind2=r['kind2'],
                        n_hprc_haps=len(hprc), n_complete=len(c), n_partial_only=len(h) - len(c),
                        n_absent=len(hprc) - len(h), n_multi_complete=sum(len(comp[k]) > 1 for k in c),
                        n_ref_path=sum(all(p == REFPATH[i] for p in comp[k]) for k in c),
                        n_nonref_path=sum(any(p != REFPATH[i] for p in comp[k]) for k in c),
                        n_distinct_paths=len({p for k in c for p in comp[k]}),
                        n_line_end=sum(any(x[11] == 'line_end' for x in st[k]) for k in h),
                        n_exits_window=sum(any(x[11].startswith('exits') for x in st[k]) for k in h),
                        grch38_status=','.join(sorted({x[5] for x in g})) or 'absent', grch38_ok=gok,
                        chm13_status=','.join(sorted({x[5] for x in ch})) or 'absent',
                        chm13_ref_path=any(x[10] == REFPATH[i] for x in ch if x[5] == 'complete')))
    with open(f'{D}/hprc_window_coverage.tsv', 'w') as w:
        w.write('\t'.join(cov_cols) + '\n')
        for c in cov:
            w.write('\t'.join(str(c[k]) for k in cov_cols) + '\n')
    nc = np.array([c['n_complete'] for c in cov])
    print(f'loci {len(loci)} coverage rows {len(cov)}; GRCh38 ok {n_ok}/{len(loci)}; HPRC haplotypes {len(hprc)}; '
          f'complete per locus median {np.median(nc)} min {nc.min()} max {nc.max()}; loci with < 80 complete '
          f'{int((nc < 80).sum())}; status counts {Counter(x[5] for x in rows)}')
    print('partial_end', Counter(x[11].split(':')[0] for x in rows if x[5] != 'complete'))
    validate(by, hprc)
    print(f'total {time.time() - t0:.0f}s')


def validate(by, hprc):
    """Spell distinct haplotype walks of 5 loci (most distinct paths SNV / INDEL, 3 random); edges must exist, lengths match the contig span, few edits vs ref."""
    oseq, succ, lev = ns['oseq'], ns['succ'], ns['lev']
    pick = sorted(by, key=lambda i: -len({x[10] for x in by[i] if x[5] == 'complete'}))
    top = [next(i for i in pick if loci[i]['kind2'] == k) for k in ('SNV', 'INDEL')]
    for i in top + [int(k) for k in np.random.default_rng(20261005).choice(len(loci), 3, replace=False)]:
        r = loci[i]
        seen = {}
        for x in by[i]:
            if x[5] == 'complete' and x[1] not in REFS and x[10] not in seen and len(seen) < 4:
                seen[x[10]] = x
        print(f"check truth {r['truth_id']} {r['chrom']}:{r['vcf_pos']} {r['vcf_ref']}>{r['vcf_alt']} {r['kind2']} "
              f"window {len(r['ref_hap'])} bp, {len({x[10] for x in by[i] if x[5] == 'complete'})} distinct paths")
        for path, x in seen.items():
            steps = [(int(s[1:]), '+' if s[0] == '>' else '-') for s in path.replace('<', ' <').replace('>', ' >').split()]
            s = ''.join(oseq(n, o) for n, o in steps)
            edges = all(steps[k + 1] in succ(*steps[k]) for k in range(len(steps) - 1))
            print(f"  {x[1]}#{x[2]} {x[3]}:{x[7]}-{x[8]} {x[6]} nodes {len(steps)} len {len(s)} span_ok "
                  f"{len(s) == x[8] - x[7]} edges_ok {edges} lev_ref {lev(s, r['ref_hap'])} lev_alt {lev(s, r['alt_hap'])}")


if __name__ == '__main__':
    main()

# Result (2026-10-05, Slurm 380997 on guinness: 24 processes, 8 min 6 s, MaxRSS 10.1 GB; scan 346 s)
#   51,142 rows (complete 38,143, left_only 6,092, right_only 5,839, inner_only 1,068; no too_long); partial_end line_end
#   11,642, exits_window 289; 90 haplotypes (88 HPRC + GRCh38 + CHM13; hprc_haplotypes.tsv identical to the HG008 one);
#   W length mismatches 0; 49,641 lines touch a window, overlapping P-fragment offsets 0.
#   GRCh38 ok 497 / 497. HPRC haplotypes with a complete traversal per locus: SNV median 85 (min 15), INDEL median 78
#   (10th pct 50, min 0); 256 loci < 80, 37 < 44 (36 INDEL, 1 SNV); 31 loci with a haplotype traversing twice.
#   Zero complete: truth 28014 (chr10:12247133 CCT>C): 85 HPRC haplotypes leave the window before the last anchor
#   4662778 (exits_window: GRCh38 carries a rare allele there; CHM13 left_only too); their partial_path still spans the
#   event window. CHM13: complete 488, absent 5, right_only 3, left_only 1; CHM13 walks the GRCh38 path at 214 loci.
#   5 checked loci: every spelled walk has existing edges and length == hap_hi0 - hap_lo0.
#   Format check vs c2: truths where >= 1 HPRC haplotype's complete local_path equals the c2 primary_path: INDEL exact
#   328 / 394, with_germline 6 / 15, closest 11 / 21; SNV exact 53 / 64, closest 2 / 3 (whole-window equality is strict:
#   other germline variants in the window change the path; carriers should be counted by elements, as in HG008 s9).

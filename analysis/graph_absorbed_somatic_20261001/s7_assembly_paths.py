"""Step 7: the d9 node walk of every assembly haplotype through each locus (giraffe), all loci of loci.tsv
(generalizes indel_graph_paths_20260930/assembly/b2_windows.py + b4_nodes.py, Illumina INDELs only).
  python s7_assembly_paths.py windows    -> asm_windows.fq (asm_seq_status.tsv rows = (locus, assembly, contig) with a
                                            minimap2 alignment; window = graph_win_t0..t1 (the GRCh38 graph window
                                            win_start0..win_end0 projected through the REF/ALT context alignment, s6)
                                            +- (80 + extrapolated bases / 4); read name truth_id|assembly|contig|start0)
  (job s7_giraffe.sh, run from tmp/graph_absorbed_somatic_20261001/: vg 1.65 giraffe on d9 + vg view -a
   -> asm_windows.gam.json; s7_parse.sh runs the next two)
  python s7_assembly_paths.py parse      -> asm_local_paths.tsv
  python s7_assembly_paths.py crosscheck -> asm_crosscheck.log (stdout): minimap2 status vs giraffe asm_event_allele, and
                                            vs the earlier Illumina INDEL results (haplotypes_normal.tsv, nodes_status.tsv
                                            recomputed on the new walks with altpaths.tsv elements)
Per window alignment, normalized to GRCh38-forward orientation (reverse-complemented if it passes first_node reversed):
  local_path      - oriented node walk from the first occurrence of first_node to the last occurrence of last_node
                    ('>n' forward, '<n' reverse); status ok / not covered (an anchor missing: the haplotype bypasses it
                    through a branch, or the window aligned elsewhere; note says which) / unmapped / inconsistent (anchors
                    in opposite orientations); for not covered, local_path_ext = the walk between the nearest GRCh38
                    nodes of the locus at or outside the missing anchor(s), and the columns below are computed on it
  n_edits_window / n_indel_edits_window - non-match edits on the walk (first..last mapping) / those changing length
  anchors_full    - the alignment covers both anchor nodes completely (then the walk spans the whole graph window)
  walk_spells     - the walk's node sequence (full nodes) = loci alt_hap / ref_hap / other ('' if anchors not full)
  walk_event_allele / asm_event_allele - over the event window (ev_lo0..ev_hi0, widened outwards to bases the walk
                    places on GRCh38 nodes): graph sequence of the walk (no edits) / assembly sequence as aligned by giraffe
                    (walk + edits), each compared with GRCh38 (REF) and GRCh38 + truth allele (ALT): ALT / REF / other
                    (walk_ext / asm_ext: 'event', else the first extension of the span towards the graph-window ends
                    (start / end / whole window) at which it is REF or ALT - the graph can place its edges outside the
                    event window; 'other' keeps 'event')
  n_edits_event / event_span0 - non-match edits on the event-window bases / their GRCh38 span; read_check - walk + edits
                    reproduce the read window
Assumptions: GRCh38 coordinates via ref() of classify.py (GRCh38 nodes visited once by GRCh38, forward); a node repeated
on the walk (repeat loop) keeps every visit, the event window takes the first / last visit of its boundary bases.
"""
import csv, json, sys
from collections import Counter, defaultdict
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
ns = {}
exec(open('/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/classify.py').read().split("fa = pysam.FastaFile(FASTA)\nchroms")[0], ns)
seq, oseq, ref = ns['seq'], ns['oseq'], ns['ref']
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
ASM = {'normal': '/scratch/jshen/data/HG008_GIAB/HG008N_curatedv6_250714_bothhaps_polished6.2.fasta.gz',
       'tumor': '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/assembly/HG008T_v3.2.fasta'}
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))


def windows():
    fa = {a: pysam.FastaFile(p) for a, p in ASM.items()}
    n = Counter()
    with open(f'{D}/asm_windows.fq', 'w') as out:
        for r in rd(f'{D}/asm_seq_status.tsv'):
            pad = 80 + int(r['ext_bases']) // 4
            a = max(0, int(r['graph_win_t0']) - pad)
            b = min(fa[r['assembly']].get_reference_length(r['contig']), int(r['graph_win_t1']) + 1 + pad)
            s = fa[r['assembly']].fetch(r['contig'], a, b).upper()
            out.write(f"@{r['truth_id']}|{r['assembly']}|{r['contig']}|{a}\n{s}\n+\n{'I' * len(s)}\n")
            n[r['assembly']] += 1
    print(dict(n), 'windows')


def revcomp_aln(maps, read):
    out = []
    for m in reversed(maps):
        n = int(m['position']['node_id']); fl = sum(int(e.get('from_length', 0)) for e in m.get('edit', []))
        out.append({'position': {'node_id': n, 'is_reverse': not m['position'].get('is_reverse', False),
                                 'offset': len(seq(n)) - int(m['position'].get('offset', 0)) - fl},
                    'edit': [dict(e, sequence=e['sequence'].translate(RC)[::-1]) if 'sequence' in e else e for e in reversed(m.get('edit', []))]})
    return out, read.translate(RC)[::-1]


def columns(maps):
    """Per walk graph base: [gcoord or None, graph base, read bases aligned to it (+ a following insertion), edit indices]."""
    C, k, pre = [], 0, ''
    for m in maps:
        n = int(m['position']['node_id']); o = '-' if m['position'].get('is_reverse') else '+'
        s = oseq(n, o); off = int(m['position'].get('offset', 0)); r = ref(n)
        for e in m.get('edit', []):
            f, t, sq = int(e.get('from_length', 0)), int(e.get('to_length', 0)), e.get('sequence')
            edit = sq is not None or f != t
            if f == 0:                                 # insertion: attach to the previous graph base
                if C: C[-1][2] += sq; C[-1][3].add(k)
                else: pre += sq
            else:
                for i in range(f):
                    g = r[1] + off + i if r and o == '+' else None
                    rb = (s[off + i] if sq is None else sq[i]) if i < t else ''
                    if i == f - 1 and t > f: rb += sq[f:]
                    C.append([g, s[off + i], rb, {k} if edit else set()])
                off += f
            k += 1
    if C: C[0][2] = pre + C[0][2]
    return C


def event_alleles(C, L0, H0, ws, ref_hap, alt_hap, lt):
    """Walk / assembly allele over the event window [L0, H0); where it is neither REF nor ALT there, retried over the span
    extended to the graph-window start, to its end, then both (the graph can spell the same haplotype with its edges outside
    the event window). -> (walk allele, its extension, assembly allele, its extension, n edits and span of the event window)."""
    pos = defaultdict(list)
    for i, c in enumerate(C):
        if c[0] is not None: pos[c[0]].append(i)
    we = ws + len(ref_hap) - 1
    left = lambda x: next((y for y in range(x, max(ws, x - 300) - 1, -1) if y in pos), None)
    right = lambda x: next((y for y in range(x, min(we, x + 300) + 1) if y in pos), None)
    out = {}
    for ext, (L, H) in (('event', (left(L0), right(H0 - 1))), ('to window start', (left(ws), right(H0 - 1))),
                        ('to window end', (left(L0), right(we))), ('whole window', (left(ws), right(we)))):
        if L is None or H is None or pos[L][0] > pos[H][-1]:
            continue
        seg = C[pos[L][0]:pos[H][-1] + 1]
        R, A = ref_hap[L - ws:H + 1 - ws], alt_hap[L - ws:H + 1 - ws + lt]
        cls = lambda x: 'ALT' if x == A else 'REF' if x == R else 'other'
        if ext == 'event':
            out['n'], out['span'] = len(set().union(*(c[3] for c in seg))), f'{L}-{H + 1}'
        for k, i in (('walk', 1), ('asm', 2)):
            a = cls(''.join(c[i] for c in seg))
            if k not in out or out[k][0] == 'other' and a != 'other':
                out[k] = (a, ext)
    if 'walk' not in out:
        return '', '', '', '', '', ''
    return out['walk'][0], out['walk'][1], out['asm'][0], out['asm'][1], out.get('n', ''), out.get('span', '')


def locate(q, maps, read):
    """Normalize the alignment to GRCh38-forward at the locus and find the walk ends: first_node .. last_node (status ok),
    else the nearest GRCh38 nodes of the locus outside the missing anchor(s) (status not covered, local_path_ext)."""
    first, last, ws, we = int(q['first_node']), int(q['last_node']), int(q['win_start0']), int(q['win_end0'])
    nodes = [int(m['position']['node_id']) for m in maps]
    rv = {n: m['position'].get('is_reverse', False) for n, m in zip(nodes, maps)}
    near = [rv[n] for n in nodes if ref(n) and ref(n)[0] == q['chrom'] and ws - 2000 <= ref(n)[1] <= we + 2000]
    if first in rv and last in rv and rv[first] != rv[last]:
        return 'inconsistent', maps, read, None, None, 'anchors in opposite orientations'
    if first not in rv and last not in rv and not near:
        return 'not covered', maps, read, None, None, 'off target (no GRCh38 node of the locus)'
    if rv.get(first, rv.get(last, Counter(near).most_common(1)[0][0] if near else False)):
        maps, read = revcomp_aln(maps, read); nodes = nodes[::-1]
    g = lambda k: ref(nodes[k]) if ref(nodes[k]) and ref(nodes[k])[0] == q['chrom'] and not maps[k]['position'].get('is_reverse') else None
    i = nodes.index(first) if first in nodes else next((k for k in range(len(nodes) - 1, -1, -1) if g(k) and g(k)[1] <= ws), None)
    if i is None: i = next((k for k in range(len(nodes)) if g(k) and ws <= g(k)[1] <= we), None)
    j = len(nodes) - 1 - nodes[::-1].index(last) if last in nodes else \
        next((k for k in range(i or 0, len(nodes)) if g(k) and g(k)[2] >= we), None) if i is not None else None
    if j is None and i is not None: j = next((k for k in range(len(nodes) - 1, i - 1, -1) if g(k) and g(k)[2] <= we), None)
    miss = ','.join(x for x, n in (('first', first), ('last', last)) if n not in nodes)
    if i is None or j is None or j < i:
        return ('inconsistent' if not miss else 'not covered'), maps, read, None, None, ('last anchor before first' if not miss else f'missing {miss}, no GRCh38 ends')
    note = '' if not miss else f'missing {miss}; ends {nodes[i]},{nodes[j]}'
    if not miss and (nodes.count(first) > 1 or nodes.count(last) > 1): note = 'anchor repeated'
    return ('ok' if not miss else 'not covered'), maps, read, i, j, note


def parse():
    LOC = {r['truth_id']: r for r in rd(f'{D}/loci.tsv')}
    cols = ['truth_id', 'assembly', 'contig', 'window_start0', 'status', 'local_path', 'n_edits_window', 'n_indel_edits_window', 'mapq',
            'anchors_full', 'walk_spells', 'walk_event_allele', 'walk_ext', 'asm_event_allele', 'asm_ext', 'n_edits_event', 'event_span0',
            'identity', 'read_check', 'local_path_ext', 'note']
    out = open(f'{D}/asm_local_paths.tsv', 'w'); out.write('\t'.join(cols) + '\n')
    for line in open(f'{D}/asm_windows.gam.json'):
        a = json.loads(line)
        tid, asm, contig, start = a['name'].split('|'); q = LOC[tid]
        o = dict(truth_id=tid, assembly=asm, contig=contig, window_start0=start, mapq=a.get('mapping_quality', 0),
                 identity=round(a.get('identity', 0), 4), note='')
        maps, read = a.get('path', {}).get('mapping', []), a.get('sequence', '')
        if not maps:
            o['status'] = 'unmapped'
            out.write('\t'.join(str(o.get(c, '')) for c in cols) + '\n'); continue
        st, maps, read, i, j, o['note'] = locate(q, maps, read)
        o['status'] = st
        if i is not None:
            W = maps[i:j + 1]
            walk = [(int(m['position']['node_id']), '-' if m['position'].get('is_reverse') else '+') for m in W]
            eds = [e for m in W for e in m.get('edit', [])]
            full = int(W[0]['position'].get('offset', 0)) == 0 and \
                int(W[-1]['position'].get('offset', 0)) + sum(int(e.get('from_length', 0)) for e in W[-1].get('edit', [])) == len(seq(walk[-1][0]))
            spelled = ''.join(oseq(*x) for x in walk)
            C = columns(W)
            r0 = sum(int(e.get('to_length', 0)) for m in maps[:i] for e in m.get('edit', []))
            r1 = r0 + sum(int(e.get('to_length', 0)) for e in eds)
            lt = len(q['vcf_alt']) - len(q['vcf_ref'])
            wa, wx, ae, ax, ne, sp = event_alleles(C, int(q['ev_lo0']), int(q['ev_hi0']), int(q['win_start0']), q['ref_hap'], q['alt_hap'], lt)
            path = ''.join(('>' if x[1] == '+' else '<') + str(x[0]) for x in walk)
            o.update(local_path=path if st == 'ok' else '', local_path_ext=path if st != 'ok' else '',
                     n_edits_window=sum(1 for e in eds if 'sequence' in e or int(e.get('from_length', 0)) != int(e.get('to_length', 0))),
                     n_indel_edits_window=sum(1 for e in eds if int(e.get('from_length', 0)) != int(e.get('to_length', 0))),
                     anchors_full=full, walk_spells=('alt_hap' if spelled == q['alt_hap'] else 'ref_hap' if spelled == q['ref_hap'] else 'other') if full and st == 'ok' else '',
                     walk_event_allele=wa, walk_ext=wx, asm_event_allele=ae, asm_ext=ax, n_edits_event=ne, event_span0=sp,
                     read_check=''.join(c[2] for c in C) == read[r0:r1])
        out.write('\t'.join(str(o.get(c, '')) for c in cols) + '\n')
    out.close()
    P = rd(f'{D}/asm_local_paths.tsv')
    print(len(P), 'window alignments', dict(Counter((p['assembly'], p['status']) for p in P)), '| read_check failures',
          sum(p['read_check'] == 'False' for p in P), '| not covered with local_path_ext', sum(p['local_path_ext'] != '' for p in P))


def crosscheck():
    LOC = {r['truth_id']: r for r in rd(f'{D}/loci.tsv')}
    S = {(r['truth_id'], r['assembly'], r['contig']): r for r in rd(f'{D}/asm_seq_status.tsv')}
    P = {(r['truth_id'], r['assembly'], r['contig']): r for r in rd(f'{D}/asm_local_paths.tsv')}
    print('== minimap2 status vs giraffe asm_event_allele (status ok rows; normal / tumor)')
    for asm in ('normal', 'tumor'):
        c = Counter((S[k]['status'], P[k]['asm_event_allele'] if P[k]['status'] == 'ok' else P[k]['status']) for k in P if k[1] == asm)
        agree = sum(v for (m, g), v in c.items() if m == g)
        print(f'  {asm}: {sum(c.values())} windows, agree (ALT=ALT, REF=REF, other=other) {agree}')
        for (m, g), v in c.most_common(): print(f'    {v:5d}  minimap2 {m:9s} giraffe {g}')
        c = Counter((S[k]['status'], P[k]['asm_event_allele']) for k in P if k[1] == asm and P[k]['status'] == 'not covered')
        print(f'    not covered rows (local_path_ext): minimap2 vs giraffe', dict(c.most_common()))
    print('== giraffe walk (graph only) vs assembly over the event window, normal ok rows:',
          dict(Counter((P[k]['walk_event_allele'], P[k]['asm_event_allele']) for k in P if k[1] == 'normal' and P[k]['status'] == 'ok').most_common()))
    E = '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/assembly'
    key = lambda r: (r['chrom'], r['vcf_pos'], r['vcf_ref'].upper(), r['vcf_alt'].upper())
    T = {key(r): r['truth_id'] for r in LOC.values()}
    Qk = {r['k']: T.get(key(r)) for r in rd(f'{E}/queries.tsv')}
    c = Counter()
    for r in rd(f'{E}/haplotypes_normal.tsv'):
        tid = Qk.get(r['k'])
        if tid: c[(r['status'], S.get((tid, 'normal', r['contig']), {}).get('status', 'no row'))] += 1
    print(f"== earlier haplotypes_normal.tsv (Illumina INDEL) vs asm_seq_status: {sum(c.values())} (truth, contig), agree "
          f"{sum(v for (a, b), v in c.items() if a == b)}")
    for (a, b), v in c.most_common(): print(f'    {v:5d}  earlier {a:9s} now {b}')
    AP = {r['k']: json.loads(r['elements']) for r in rd(f'{E}/altpaths.tsv')}
    c = Counter()
    for r in rd(f'{E}/nodes_status.tsv'):
        tid = Qk.get(r['k'])
        if not tid or not r['event_hap'].startswith('chr'): continue
        oth = r['event_hap'][:-1] + ('2' if r['event_hap'].endswith('1') else '1')
        els = [e for e in AP[r['k']] if e['role'] == 'indel']
        for hap, old in ((r['event_hap'], r['event_hap_walks']), (oth, r['other_hap_walks'])):
            p = P.get((tid, 'normal', hap))
            if not p or p['status'] != 'ok':
                new = 'not covered' if p and p['status'] == 'not covered' else p['status'] if p else 'no window'
            else:
                path = [int(x) for x in p['local_path'].replace('<', '>').split('>')[1:]]
                on, pairs = set(path), set(zip(path, path[1:])) | set(zip(path[1:], path))
                if not all(e['x'] in on and e['y'] in on for e in els): new = 'not covered'
                else:
                    w = [(all(n in on for n, _ in e['nodes']) if e['type'] == 'branch' else (e['x'], e['y']) in pairs) for e in els]
                    new = 'all' if all(w) else 'some' if any(w) else 'none'
            c[(old, new)] += 1
    print(f"== earlier nodes_status.tsv event/other-hap walks the ALT-path elements vs the new walks: {sum(c.values())}, agree "
          f"{sum(v for (a, b), v in c.items() if a == b)}")
    for (a, b), v in c.most_common(): print(f'    {v:5d}  earlier {a:11s} now {b}')


if __name__ == '__main__':
    {'windows': windows, 'parse': parse, 'crosscheck': crosscheck}[sys.argv[1]]()

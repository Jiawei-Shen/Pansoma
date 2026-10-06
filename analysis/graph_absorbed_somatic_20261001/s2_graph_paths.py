"""Step 2: the d9 graph paths that spell each truth ALT haplotype, and the branch elements they use.

Input: loci.tsv (s1). For each locus, oriented paths from (first_node, +) to (last_node, +) whose GRCh38 nodes lie in
win_start0..win_end0 of the same chromosome (non-GRCh38 nodes free; at most 600 nodes; a GRCh38 node walked in reverse
counts as a branch node):
  exact          - paths spelling alt_hap (DFS with prefix pruning, b1_altpaths find_path / classify.py; ALL distinct
                   matching paths up to 50, enumeration cap 50,000 paths reaching last_node, 3 M node expansions)
  with_germline  - none spells alt_hap, but alt_hap + k of the <= 6 PASS HG008-N dipcall germline alleles inside the
                   window nearest the truth (all combinations, smallest k that matches) is spelled (b1 'matches with
                   nearby germline'); germline_used / germline_gt (dipcall GT, hap1|hap2) of the primary path
  closest        - else the graph haplotypes at minimum edit distance to alt_hap (branch-and-bound DFS with a
                   vectorised Levenshtein row per path prefix: prune when min(row) > best, best starts at the GRCh38
                   path's distance = grch38_edit_distance; ties kept up to 50; cap 50,000 completed / 300,000
                   expansions / 900 s); closest_edit_distance == grch38_edit_distance means the graph has nothing
                   closer to the ALT than GRCh38 itself
  none           - no path completed
Elements of a path (b1 elements(), extended): 'branch' = run of non-GRCh38 nodes between GRCh38 nodes x and y (seq,
replaced = GRCh38 bases between x and y, start0 = x end + 1, end0 = y start), 'skip' = edge x -> y with replaced != 0
(> 0 deletion, < 0 back edge = repeat copy re-walked). role snv = same-length branch, else indel. subtype: snv_branch
(1 bp) / mnv_branch / ins_branch (replaced 0) / replacing_branch / back_branch (branch with replaced < 0) / del_skip /
back_edge. in_event: [min(start0,end0), max(start0,end0)) overlaps the event window [ev_lo0, ev_hi0); a point counts if
ev_lo0 <= point <= ev_hi0. id 'B:x:y:>n1>n2' or 'S:x:y'.
local_alt / local_ref: the path spells the ALT / GRCh38 sequence of the event window +-10 bp (alt_hap[ev_lo0 - 10 ..
ev_hi0 + 10 + len(ALT) - len(REF)), ref_hap[ev_lo0 - 10 .. ev_hi0 + 10)) as a substring (always local_alt for exact).
event_net = net length (len(seq) - replaced) of the primary's in_event elements vs truth_net = len(ALT) - len(REF):
they differ where the graph puts part of the event in an adjacent repeat stretch that the s1 event window misses
(compound repeats such as (TC)n(TG)n, period > 6 VNTRs) or where a germline indel is in the event window.
Primary path = fewest in_event elements, then fewest nodes (then DFS order); for closest, first a local_alt path, then
the tied path farthest from GRCh38 (largest d_grch38 = edit distance to ref_hap; otherwise a tie with the plain GRCh38
path would always win). Paths in local format (>n<n, GRCh38-forward).
Outputs: graph_paths.tsv (one row per locus, loci.tsv order) and graph_elements.tsv (one row per truth x element over
all kept paths). Run: python s2_graph_paths.py [n_workers]  (Slurm: tmp/graph_absorbed_somatic_20261001/s2.sbatch);
python s2_graph_paths.py report  re-prints the summary and the cross-check with indel_graph_paths_20260930 (Illumina).
Assumptions: the nearest-6 germline rule replaces b1's first-6-in-window; germline phase is not used to filter combos.
"""
import csv, itertools, json, sqlite3, sys, time
from bisect import bisect_left
from collections import Counter
from multiprocessing import Pool
import numpy as np
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
P = '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930'
ns = {}
exec(open(f'{P}/classify.py').read().split("fa = pysam.FastaFile(FASTA)\nchroms")[0], ns)
oseq, succ, ref = ns['oseq'], ns['succ'], ns['ref']
GERM_VCF = '/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz'  # labels manifest germline
MAX_KEEP, CAP, MAX_NODES, MAX_EXPAND, MAX_EXPAND_CLOSEST, MAX_SECONDS, N_GERM = 50, 50000, 600, 3_000_000, 300_000, 900, 6
tag = lambda path: ''.join(('>' if o == '+' else '<') + str(n) for n, o in path)
spell = lambda path: ''.join(oseq(n, o) for n, o in path)


def lev(a, b):
    """Levenshtein distance, one numpy row per character of a."""
    B = np.frombuffer(b.encode(), dtype=np.uint8); ar = np.arange(len(B) + 1, dtype=np.int32); row = ar
    for ch in a.encode():
        t = np.empty_like(row); t[0] = row[0] + 1
        np.minimum(row[1:] + 1, row[:-1] + (B != ch), out=t[1:])
        row = np.minimum.accumulate(t - ar) + ar
    return int(row[-1])


def walk(L, targets=None, A=None):
    """DFS over the window. targets: list of sequences (prefix-pruned exact search) or None with A (numpy alt bytes):
    branch-and-bound closest search. Returns (hits [(path, target or dist)], n_completed, cap_info)."""
    first, last, chrom, lo, hi = L['first'], L['last'], L['chrom'], L['lo'], L['hi']
    hits, done, nexp, info, best, t0 = [], 0, 0, set(), None, time.time()
    max_exp = MAX_EXPAND if A is None else MAX_EXPAND_CLOSEST
    if A is not None:
        ar = np.arange(len(A) + 1, dtype=np.int32); mm = {}
        def step(row, s):
            for ch in s.encode():
                if ch not in mm:
                    mm[ch] = (A != ch).astype(np.int32)
                t = np.empty_like(row); t[0] = row[0] + 1
                np.minimum(row[1:] + 1, row[:-1] + mm[ch], out=t[1:])
                row = np.minimum.accumulate(t - ar) + ar
            return row
        best = L['d_ref']                                          # the GRCh38 path's distance bounds the search
        stack = [((first, '+'), step(ar, oseq(first, '+')), ((first, '+'),))]
    else:
        s0 = oseq(first, '+')
        stack = [((first, '+'), len(s0), tuple(t for t in targets if t.startswith(s0)), ((first, '+'),))]
    while stack:
        if done >= CAP or nexp >= max_exp or (nexp % 1000 == 0 and time.time() - t0 > MAX_SECONDS):
            info.add('enum_cap' if done >= CAP else 'expand_cap' if nexp >= max_exp else 'time_cap'); break
        if A is None:
            (n, o), pos, cands, path = stack.pop()
        else:
            (n, o), row, path = stack.pop()
            if int(row.min()) > best:
                continue
        if n == last:
            if o == '+':
                done += 1
                if A is None:
                    for t in cands:
                        if len(t) == pos:
                            hits.append((path, t)); break
                    if len(hits) > MAX_KEEP:
                        hits.pop(); info.add('more_than_50'); break
                else:
                    d = int(row[-1])
                    if d < best:
                        best, hits = d, []
                    if d == best:
                        if len(hits) < MAX_KEEP:
                            hits.append((path, d))
                        else:
                            info.add('more_than_50')
            continue
        if len(path) >= MAX_NODES:
            info.add('max_nodes'); continue
        nexp += 1
        kids = []
        for m, mo in succ(n, o):
            r = ref(m)
            if r and (r[0] != chrom or r[1] < lo or r[2] > hi):
                continue
            x = oseq(m, mo)
            if A is None:
                c = tuple(t for t in cands if t[pos:pos + len(x)] == x)
                if c:
                    kids.append(((m, mo), pos + len(x), c, path + ((m, mo),)))
            else:
                rw = step(row, x)
                if int(rw.min()) <= best:
                    kids.append((-int(rw.min()), ((m, mo), rw, path + ((m, mo),))))
        if A is not None:                                        # lowest bound popped first
            kids = [k for _, k in sorted(kids, key=lambda z: z[0])]
        stack += kids
    return hits, done, ';'.join(sorted(info))


def elements(path, ev_lo, ev_hi):
    out, last, lastnode, run = [], None, None, []
    for n, o in path:
        r = ref(n) if o == '+' else None
        if r is None:
            run.append((n, o)); continue
        if last is not None and (run or r[1] != last[2] + 1):
            rep, s = r[1] - last[2] - 1, ''.join(oseq(*x) for x in run)
            e = dict(type='branch' if run else 'skip', x=lastnode, y=n, nodes=tag(run), seq=s, replaced=rep,
                     start0=last[2] + 1, end0=r[1], role='snv' if run and len(s) == rep else 'indel')
            e['id'] = f"B:{lastnode}:{n}:{e['nodes']}" if run else f'S:{lastnode}:{n}'
            a, b = min(e['start0'], e['end0']), max(e['start0'], e['end0'])
            e['in_event'] = (ev_lo <= a <= ev_hi) if a == b else (a < ev_hi and b > ev_lo)
            e['subtype'] = ('del_skip' if rep > 0 else 'back_edge') if not run else 'ins_branch' if rep == 0 else \
                ('snv_branch' if rep == 1 else 'mnv_branch') if len(s) == rep else 'back_branch' if rep < 0 else 'replacing_branch'
            out.append(e)
        last, lastnode, run = r, n, []
    return out


def apply(region, base, edits):
    s, cur = [], 0
    for p0, a, b in sorted(edits):
        off = p0 - base
        if off < cur or off + len(a) > len(region):
            return None
        s += [region[cur:off], b]; cur = off + len(a)
    return ''.join(s) + region[cur:]


def init():
    ns['db'] = sqlite3.connect(f"file:{ns['IDX']}?mode=ro", uri=True)   # one connection per worker


def locus(L):
    t0 = time.time()
    out = dict(truth_id=L['truth_id'], chrom=L['chrom'], vcf_pos=L['vcf_pos'], kind2=L['kind2'], match='none', germline_used='',
               germline_gt='', closest_edit_distance='', n_alt_paths=0, cap_hit=False, cap_info='', n_completed=0)
    L = dict(L, first=int(L['first_node']), last=int(L['last_node']), lo=int(L['win_start0']), hi=int(L['win_end0']))
    hits, done, info = walk(L, [L['alt_hap']])
    infos = [f'exact:{info}'] if info else []
    labels = {L['alt_hap']: ('', '')}
    if hits:
        out['match'] = 'exact'
    else:
        base, pos0 = L['lo'], int(L['vcf_pos']) - 1
        som = (pos0, L['vcf_ref'].upper(), L['vcf_alt'].upper())
        g = sorted(L['germ'], key=lambda v: abs(v[0] - pos0))[:N_GERM]
        for k in range(1, len(g) + 1):
            labels = {}
            for combo in itertools.combinations(g, k):
                h = apply(L['ref_hap'], base, [som] + [c[:3] for c in combo])
                if h:
                    labels.setdefault(h, (';'.join(f'{c[0] + 1}:{c[1]}>{c[2]}' for c in combo), ';'.join(c[3] for c in combo)))
            if not labels:
                continue
            hits, done, info = walk(L, sorted(labels))
            if info:
                infos.append(f'germline{k}:{info}')
            if hits:
                out['match'] = 'with_germline'; break
        if not hits:
            L['d_ref'] = out['grch38_edit_distance'] = lev(L['ref_hap'], L['alt_hap'])
            hits, done, info = walk(L, A=np.frombuffer(L['alt_hap'].encode(), dtype=np.uint8))
            if info:
                infos.append(f'closest:{info}')
            if hits:
                out['match'] = 'closest'; out['closest_edit_distance'] = hits[0][1]
    ev_lo, ev_hi = int(L['ev_lo0']), int(L['ev_hi0'])
    a, b, dl = max(0, ev_lo - 10 - L['lo']), ev_hi + 10 - L['lo'], len(L['vcf_alt']) - len(L['vcf_ref'])
    loc_alt, loc_ref = L['alt_hap'][a:b + dl], L['ref_hap'][a:b]
    paths = []
    for path, t in hits:
        el, sp = elements(path, ev_lo, ev_hi), spell(path)
        dg = lev(sp, L['ref_hap']) if out['match'] == 'closest' else 0
        paths.append((loc_alt not in sp, -dg, sum(e['in_event'] for e in el), len(path), path, el, t, loc_ref in sp))
    out.update(n_alt_paths=len(paths), n_completed=done, cap_info='|'.join(infos), cap_hit=bool(infos))
    elem, rows = {}, []
    if paths:
        prim = min(range(len(paths)), key=lambda i: paths[i][:4])
        nla, mdg, _, nn, path, el, t, lr = paths[prim]
        out.update(primary_local_alt=not nla, primary_local_ref=lr, truth_net=dl,
                   event_net=sum(len(e['seq']) - e['replaced'] for e in el if e['in_event']))
        if out['match'] == 'with_germline':
            out['germline_used'], out['germline_gt'] = labels[t]
        if out['match'] == 'closest':
            out['primary_d_grch38'] = -mdg
        sets, dgs = Counter(), {}
        for p in paths:
            s = tuple(sorted(e['id'] for e in p[5] if e['in_event'])); sets[s] += 1
            dgs[s] = (max(dgs.get(s, (0, False))[0], -p[1]), dgs.get(s, (0, False))[1] or not p[0])
        es = [dict(elements=list(s), n_paths=c, **(dict(d_grch38=dgs[s][0], local_alt=dgs[s][1]) if out['match'] == 'closest' else {}))
              for s, c in sets.most_common()]
        out.update(primary_path=tag(path), primary_n_nodes=nn, primary_event_ids=';'.join(sorted(e['id'] for e in el if e['in_event'])),
                   n_event_sets=len(sets), event_element_sets=json.dumps(es),
                   primary_elements=json.dumps([{k: e[k] for k in ('id', 'type', 'x', 'y', 'nodes', 'seq', 'replaced', 'start0', 'end0', 'role', 'in_event', 'subtype')} for e in el]))
        for i, p in enumerate(paths):
            for e in {e['id']: e for e in p[5]}.values():
                r = elem.setdefault(e['id'], dict(e, truth_id=L['truth_id'], element_id=e['id'], n_paths_with_element=0, in_primary=False))
                r['n_paths_with_element'] += 1; r['in_primary'] |= i == prim
        rows = list(elem.values())
    out['seconds'] = round(time.time() - t0, 2)
    return out, rows


PCOLS = ['truth_id', 'match', 'germline_used', 'closest_edit_distance', 'n_alt_paths', 'cap_hit', 'primary_path', 'primary_elements',
         'event_element_sets', 'chrom', 'vcf_pos', 'kind2', 'germline_gt', 'cap_info', 'n_completed', 'primary_n_nodes',
         'n_event_sets', 'primary_event_ids', 'grch38_edit_distance', 'primary_d_grch38', 'primary_local_alt', 'primary_local_ref',
         'truth_net', 'event_net', 'seconds']
ECOLS = ['truth_id', 'element_id', 'type', 'x', 'y', 'nodes', 'seq', 'replaced', 'start0', 'end0', 'role', 'in_event',
         'n_paths_with_element', 'in_primary', 'subtype']


def run(nproc):
    loci = list(csv.DictReader(open(f'{D}/loci.tsv'), delimiter='\t'))
    iv = {}
    for L in loci:
        L['germ'] = []
        iv.setdefault(L['chrom'], []).append((int(L['win_start0']), int(L['win_end0']), L))
    for c in iv:
        iv[c].sort(key=lambda z: z[0])
    starts = {c: [z[0] for z in v] for c, v in iv.items()}
    span = max(b - a for v in iv.values() for a, b, _ in v)
    n_g = 0
    for rec in pysam.VariantFile(GERM_VCF):                       # PASS germline alleles inside each window (one pass)
        v = iv.get(rec.chrom)
        if not v:
            continue
        p0, a = rec.pos - 1, rec.ref.upper()
        j = bisect_left(starts[rec.chrom], p0 - span)
        hit = [L for s, e, L in v[j:bisect_left(starts[rec.chrom], p0 + 1)] if s <= p0 and p0 + len(a) <= e + 1]
        if not hit or (rec.filter.keys() and 'PASS' not in rec.filter.keys()):
            continue
        gt = '|'.join('.' if x is None else str(x) for x in rec.samples[0]['GT'])
        for b in rec.alts or ():
            if b and b[0] not in '<*':
                n_g += 1
                for L in hit:
                    L['germ'].append((p0, a, b.upper(), gt))
    print(len(loci), 'loci,', n_g, 'germline allele x window hits', flush=True)
    t0 = time.time()
    with Pool(nproc, initializer=init) as pool, open(f'{D}/graph_paths.tsv', 'w') as wp, open(f'{D}/graph_elements.tsv', 'w') as we:
        wp.write('\t'.join(PCOLS) + '\n'); we.write('\t'.join(ECOLS) + '\n')
        for k, (o, rows) in enumerate(pool.imap(locus, loci, chunksize=1)):
            wp.write('\t'.join(str(o.get(c, '')) for c in PCOLS) + '\n')
            for r in rows:
                we.write('\t'.join(str(r[c]) for c in ECOLS) + '\n')
            if k % 100 == 0:
                print(k, len(loci), round(time.time() - t0), 's', flush=True)
    report()


def report():
    G = list(csv.DictReader(open(f'{D}/graph_paths.tsv'), delimiter='\t'))
    E = list(csv.DictReader(open(f'{D}/graph_elements.tsv'), delimiter='\t'))
    V = {r['truth_id']: r for r in csv.DictReader(open(f'{D}/variant_set.tsv'), delimiter='\t')}
    print('rows', len(G), 'elements', len(E))
    print('match x kind', sorted(Counter((g['kind2'], g['match']) for g in G).items()))
    print('cap_hit', sorted(Counter((g['kind2'], g['match'], g['cap_info']) for g in G if g['cap_hit'] == 'True').items()))
    print('n_alt_paths', sorted(Counter((g['kind2'], min(int(g['n_alt_paths']), 10)) for g in G).items()))
    print('ambiguous event sets (>1)', sorted(Counter((g['kind2'], g['match'], int(g['n_event_sets'] or 0) > 1) for g in G).items()))
    print('primary in-event element count', sorted(Counter((g['kind2'], len(g['primary_event_ids'].split(';')) if g['primary_event_ids'] else 0) for g in G).items()))
    kd = {g['truth_id']: g['kind2'] for g in G}
    print('element subtypes, primary in_event', sorted(Counter((kd[e['truth_id']], e['subtype']) for e in E if e['in_primary'] == 'True' and e['in_event'] == 'True').items()))
    print('element subtypes, primary out of event', sorted(Counter((kd[e['truth_id']], e['subtype']) for e in E if e['in_primary'] == 'True' and e['in_event'] != 'True').items()))
    print('event_net == truth_net', sorted(Counter((g['kind2'], g['match'], g['event_net'] == g['truth_net']) for g in G).items()))
    print('closest distances', sorted(Counter((g['kind2'], g['closest_edit_distance']) for g in G if g['match'] == 'closest').items()))
    print('closest: (no better than GRCh38, primary local_alt, local_ref)', sorted(Counter((g['kind2'], g['closest_edit_distance'] == g['grch38_edit_distance'],
          g['primary_local_alt'], g['primary_local_ref']) for g in G if g['match'] == 'closest').items()))
    print('with_germline: (primary local_alt, local_ref)', sorted(Counter((g['kind2'], g['primary_local_alt'], g['primary_local_ref']) for g in G if g['match'] == 'with_germline').items()))
    for p in ('PacBio', 'ONT', 'Illumina'):
        print(p, sorted(Counter((g['kind2'], g['match']) for g in G if V[g['truth_id']][f'{p}_perfect'] == 'True').items()))
    # cross-check with the Illumina-only INDEL analysis (indel_graph_paths_20260930): classes_all.tsv and assembly/altpaths.tsv
    key = lambda r: (r['chrom'], r['vcf_pos'], r['vcf_ref'], r['vcf_alt'])
    byk = {key(V[g['truth_id']]): g for g in G}
    C = list(csv.DictReader(open(f'{P}/classes_all.tsv'), delimiter='\t'))
    print('classes_all overlap', sum(key(c) in byk for c in C), 'of', len(C))
    print('graph_class x match', sorted(Counter((c['graph_class'], byk[key(c)]['match']) for c in C if key(c) in byk).items()))
    sub = {g['truth_id']: Counter() for g in G}
    for e in E:
        if e['in_primary'] == 'True' and e['in_event'] == 'True':
            sub[e['truth_id']][e['subtype']] += 1
    print('graph_class x primary in-event subtypes', sorted(Counter((c['graph_class'], tuple(sorted(sub[byk[key(c)]['truth_id']].elements())))
                                                                    for c in C if key(c) in byk).items(), key=lambda z: -z[1])[:25])
    Q = {q['k']: q for q in csv.DictReader(open(f'{P}/assembly/queries.tsv'), delimiter='\t')}
    eid = {}
    for e in E:
        eid.setdefault(e['truth_id'], set()).add((e['type'], e['x'], e['y'], e['seq']))
    cmp = Counter()
    for a in csv.DictReader(open(f'{P}/assembly/altpaths.tsv'), delimiter='\t'):
        g = byk.get(key(Q[a['k']]))
        if g is None or a['path_found'] != 'True':
            cmp['not in set / no path'] += 1; continue
        theirs = {(e['type'], str(e['x']), str(e['y']), e['seq']) for e in json.loads(a['elements']) if e['role'] == 'indel'}
        mine = eid.get(g['truth_id'], set())
        cmp[g['match'] + ': ' + ('all b1 indel elements among ours' if theirs <= mine else 'b1 element missing')] += 1
    print('altpaths.tsv indel elements vs ours', sorted(cmp.items()))


if __name__ == '__main__':
    if sys.argv[1:] == ['report']:
        report()
    else:
        run(int(sys.argv[1]) if sys.argv[1:] else 8)

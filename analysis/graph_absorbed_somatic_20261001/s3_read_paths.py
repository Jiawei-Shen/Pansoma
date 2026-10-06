"""Step 3: the graph nodes the reads actually traverse at each locus of variant_set.tsv (per platform, per spanning read).

Inputs: $D/variant_set.tsv + loci.tsv (ev_lo0 / ev_hi0); the platform's read-level miss analysis reads.py and report.py
(tmp/somatic_miss_analysis_{v6,ont,illumina}, = analysis/hg008_somatic_miss_20260926/<platform>/scripts); the frozen run
source, truth table and target-node parts it names, which moved on 2026-10-01 to
pansoma_v2_tensors/backup_ch6_linear100_20261001/HG008T_<P>/ (run/source, run/parts, truth/somatic.graph.tsv); the
platform's sorted GAM (path unchanged in reads.py).
Method (reads.py's code, not a re-implementation):
  * reads.py is exec'd as a module with its stale path constants remapped (REMAP; every old string is asserted to occur);
  * its analyse() runs unchanged except for one inserted call _hook(locals(), call) after each read call it makes
    (ref_exact, ref_like / unclear, alt_exact / alt_like); records that are not spanning are not recorded;
  * the hook execs reads.py's own site block (from `ordered = between[::-1] ...` to the `kind = ...` statement, cut out
    of its source) on the read's locals, so REF reads get the has-edit flag and representation kind that reads.py gives
    ALT reads (for ALT reads the block's kind / local_edit are asserted equal to reads.py's own);
    has_edit_in_window = local_edit or local_obs (the edit / candidate test of reads.py's ':no_edit' suffix, +-5 bp);
  * local_path: the trimmed GAM record's mappings (reads.py trim(), +-600 read bases), consecutive mappings that continue
    one node visit merged, turned to GRCh38-forward order when the read's left anchor is reversed relative to GRCh38
    (reads.py flip_l: list reversed, every orientation flipped), cut from the last GRCh38 node (visited once by GRCh38,
    same chromosome) holding a base before the event window (start0 < ev_lo0) to the first GRCh38 node from there on
    holding a base after it (end0 >= ev_hi0), both inclusive: the shortest walk that covers the event window (d9 GRCh38
    nodes are up to 1 kb long, so an anchor node often reaches into the window; a read that stays inside one GRCh38 node
    gives a one-node walk). spans_event = both anchors found; otherwise the cut runs to that end of the record;
    path_type = branch (a non-GRCh38 node, or a GRCh38 node traversed against GRCh38, inside the cut), skip_edge
    (GRCh38 nodes only, two consecutive ones not adjacent in GRCh38), grch38;
  * majority_path = the commonest local_path among ALT reads with no edit in the window and spans_event;
    majority_seq = ALT / REF / other: its sequence (classify.py oseq) against GRCh38 between its two anchor nodes with /
    without the truth allele ('other' = the reads also follow a nearby non-GRCh38 branch, e.g. a germline SNP);
  * reason_recomputed = report.py primary() (majority representation of the exact ALT reads, else all ALT reads) on the
    locus record analyse() returns; compared with P_reason of variant_set.tsv (set where P_status = no_candidate).
Modes:  run P K N [LIMIT]  -> $D/read_paths_parts/P/part{K}of{N}.{rows.tsv,loci.jsonl} (LIMIT: a seeded sample of LIMIT
                              no_candidate loci, files test.*; env SUB=j/m: sub-part j of m of part K, files
                              part{K}of{N}s{j}of{m}.*); workers = SLURM_CPUS_PER_TASK, loci in chunks of 5 (reads.py)
        rest-list P TAG    -> $D/read_paths_parts/P/rest{TAG}.ids: the platform's loci with no complete loci.jsonl line in
                              any part*/rest* file yet (fixed once, so every array task splits the same list)
        rest P TAG K N     -> loci K::N of rest{TAG}.ids -> rest{TAG}_{K}of{N}.{rows.tsv,loci.jsonl} (re-run of what a
                              timed-out job left; resubmit rest-list + rest until nothing is left)
        Every run writes a chunk (2 loci) as soon as it completes, rows first, then its loci lines, so a job killed at
        its time limit keeps all finished loci (a locus counts only once its loci line is complete).
        merge P            -> $D/read_paths_P.tsv.gz (one row per spanning record: truth_id read call has_edit_in_window
                              kind mapq spans_event local_path path_type reverse) and $D/read_path_summary_P.tsv
Assumptions: a record is one GAM alignment (two mates of a pair can give two rows with one read name); MAPQ > 10 and the
anchor rule are reads.py's; PacBio SNV reasons in variant_set came from the v5 analysis (v5 decoder), recomputed here with
the v6 decoder of the run.
"""
import ast, collections, csv, gzip, json, os, random, re, sys, textwrap, types
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
B = '/scratch/jshen/data/pansoma_v2_tensors/backup_ch6_linear100_20261001'
T = '/scratch/jshen/Github/Pansoma/tmp'
OLD = '/scratch/jshen/data/pansoma_v2_tensors/'
CLASSIFY = '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/classify.py'
PLAT = {'PacBio': ('somatic_miss_analysis_v6', 'Liss_lab_PacBio_Revio_20240125', 'v6_run/'),
        'ONT': ('somatic_miss_analysis_ont', 'Liss_lab_Northeastern-ONT-UL-20241216', 'v3_run/'),
        'Illumina': ('somatic_miss_analysis_illumina', 'Liss_lab_BCM_Illumina-WGS_20240313', 'v3_run/')}
csv.field_size_limit(sys.maxsize)


def load(p):
    """reads.py of platform p as module R (paths remapped, analyse() instrumented) and report.py primary()."""
    d, old, run = PLAT[p]
    f = f'{T}/{d}/reads.py'
    src = open(f).read()
    remap = [(OLD + old, f'{B}/HG008T_{p}'), (run, 'run/')]
    if p == 'PacBio':
        remap.append((OLD + 'truth/somatic.graph.tsv', f'{B}/HG008T_PacBio/truth/somatic.graph.tsv'))
    for a, b in remap:
        assert a in src, a
        src = src.replace(a, b)
    R = types.ModuleType(f'reads_{p}'); R.__file__ = f; sys.modules[R.__name__] = R
    exec(compile(src, f, 'exec'), R.__dict__)
    tree = ast.parse(src)
    fn = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    lines = ast.get_source_segment(src, fn['analyse']).splitlines()
    at = lambda s: [i for i, l in enumerate(lines) if l.strip().startswith(s)]
    i0, = at('ordered = between[::-1] if flip_l else between')
    i1, = at('info["alt_representation"][kind')
    R._SEG = compile(textwrap.dedent('\n'.join(lines[i0:i1])), f + ':site_block', 'exec')
    for s, call in [('info["ref_exact"] += 1', '"ref_exact"'), ('info["ref_like" if dr < da else "unclear"] += 1',
                    '"ref_like" if dr < da else "unclear"'), ('info["alt_representation"][kind', '"alt_exact" if exact else "alt_like"')]:
        i, = at(s)
        ind = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
        lines.insert(i + 1, f'{ind}_hook(locals(), {call})')
    exec(compile('\n'.join(lines), f + ':analyse_instrumented', 'exec'), R.__dict__)
    main = [n for n in ast.walk(fn['main']) if isinstance(n, ast.Assign) and getattr(n.targets[0], 'id', '') == 'targets']
    R._TARGETS_EXPR = ast.get_source_segment(src, main[0].value)
    rep = open(f'{T}/{d}/report.py').read()
    ns = {'Counter': collections.Counter}
    exec(ast.get_source_segment(rep, [n for n in ast.parse(rep).body if getattr(n, 'name', '') == 'primary'][0]), ns)
    R._hook, R._ROWS = hook, []
    return R, ns['primary']


def hook(L, call):
    R = MOD
    m, x = L['m'], L['x']
    ns = dict(R.__dict__); ns.update(L)
    exec(R._SEG, ns)
    if call.startswith('alt'):
        assert ns['kind'] == L['kind'] and ns['local_edit'] == L['local_edit'], (m['truth_id'], x.name)
    has_edit = bool(ns['local_edit'] or ns['local_obs'])
    visits, chrom_of, start0, rev_of, lengths, contig = (L[k] for k in ('visits', 'chrom_of', 'start0', 'rev_of', 'lengths', 'contig'))
    flip = bool(L['flip_l'])
    vis = []
    for mp in x.path.mapping:
        n, r, off = int(mp.position.node_id), bool(mp.position.is_reverse), int(mp.position.offset)
        fl = sum(e.from_length for e in mp.edit)
        if vis and vis[-1][0] == n and vis[-1][1] == r and off == vis[-1][2] and off > 0:
            vis[-1][2] = off + fl; continue                      # same node visit split over two mappings
        vis.append([n, r, off + fl])
    v = [(n, r != flip) for n, r, _ in (vis[::-1] if flip else vis)]
    co = [(int(start0[n]), int(start0[n] + lengths[n] - 1)) if visits[n] == 1 and chrom_of[n] == contig else None for n, _ in v]
    lo0, hi0 = m['ev_lo0'], m['ev_hi0']
    e0 = min((i for i, c in enumerate(co) if c and c[1] >= hi0), default=len(co) - 1)
    b = max((i for i, c in enumerate(co[:e0 + 1]) if c and c[0] < lo0), default=None)
    e = min((i for i, c in enumerate(co) if c and c[1] >= hi0 and (b is None or i >= b)), default=None)
    i0, i1 = (0 if b is None else b), (len(v) - 1 if e is None else e)
    sub = list(range(i0, i1 + 1))
    if any(co[i] is None or v[i][1] != bool(rev_of[v[i][0]]) for i in sub):
        ptype = 'branch'
    elif any(co[i + 1][0] != co[i][1] + 1 for i in sub[:-1]):
        ptype = 'skip_edge'
    else:
        ptype = 'grch38'
    walk = ''.join(('<' if v[i][1] else '>') + str(v[i][0]) for i in sub)
    R._ROWS.append((m['truth_id'], x.name, call, has_edit, ns['kind'], x.mapping_quality, b is not None and e is not None,
                    walk, ptype, int(flip)))


def work(chunk):
    MOD._ROWS = []
    res = MOD.analyse(chunk)
    return res, MOD._ROWS


def loci(p):
    """variant_set rows joined to the truth table (prepared as reads.py loci() does) and to the event window of loci.tsv."""
    vs = {r['truth_id']: r for r in csv.DictReader(open(f'{D}/variant_set.tsv'), delimiter='\t')}
    lc = {r['truth_id']: r for r in csv.DictReader(open(f'{D}/loci.tsv'), delimiter='\t')}
    out = []
    for row in csv.DictReader(open(MOD.TRUTH), delimiter='\t'):
        t = vs.get(row['truth_id'])
        if t is None:
            continue
        assert all(row[k] == t[k] for k in ('chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt')), row['truth_id']
        row.update(status=t[f'{p}_status'], pos0=int(row['pos0']), keys=row['keys'].split(',') if row['keys'] else [],
                   ev_lo0=int(lc[row['truth_id']]['ev_lo0']), ev_hi0=int(lc[row['truth_id']]['ev_hi0']))
        out.append(row)
    assert len(out) == len(vs)
    return out


def run(p, k, n, limit=None):
    work_loci = sorted(loci(p), key=lambda m: (m['chrom'], m['pos0']))
    if limit:
        work_loci = sorted(random.Random(3).sample([m for m in work_loci if m['status'] == 'no_candidate'], limit),
                           key=lambda m: (m['chrom'], m['pos0']))
        stem = 'test'
    else:
        work_loci, stem = work_loci[k::n], f'part{k}of{n}'
        if os.environ.get('SUB'):                       # 'j/m': sub-part j of m of part k (smaller Slurm jobs)
            j, m_ = map(int, os.environ['SUB'].split('/'))
            work_loci, stem = work_loci[j::m_], f'{stem}s{j}of{m_}'
    process(p, work_loci, stem)


def process(p, work_loci, stem):
    targets = eval(MOD._TARGETS_EXPR, MOD.__dict__)
    print(p, stem, len(work_loci), 'loci', len(targets), 'target nodes', flush=True)
    out = f'{D}/read_paths_parts/{p}'; os.makedirs(out, exist_ok=True)
    chunks = [work_loci[i:i + 2] for i in range(0, len(work_loci), 2)]
    done = 0
    with ProcessPoolExecutor(int(os.environ.get('SLURM_CPUS_PER_TASK', 8)), initializer=MOD.init, initargs=(targets, {})) as pool, \
            open(f'{out}/{stem}.rows.tsv', 'w') as wr, open(f'{out}/{stem}.loci.jsonl', 'w') as wl:
        for fut in as_completed([pool.submit(work, c) for c in chunks]):
            res, rows = fut.result()
            for r in rows:
                wr.write('\t'.join(map(str, r)) + '\n')
            wr.flush()
            for r in res:
                wl.write(json.dumps(r) + '\n')
            wl.flush()
            done += len(res)
            if done % 50 < 2:
                print(f'{done}/{len(work_loci)} loci', flush=True)


def read_loci_lines(path):
    """Complete loci.jsonl lines of one file (a job killed mid-write leaves a truncated last line)."""
    out = []
    for line in open(path):
        try:
            out.append(json.loads(line))
        except ValueError:
            pass
    return out


def stems_of(p):
    out = f'{D}/read_paths_parts/{p}'
    st = {f[:-len('.loci.jsonl')] for f in os.listdir(out) if f.endswith('.loci.jsonl') and f.startswith(('part', 'rest'))}
    return sorted(st, key=lambda s: (not s.startswith('part'), s))     # part* first, then rest* in order


def rest_list(p, tag):
    out = f'{D}/read_paths_parts/{p}'
    done = {r['truth_id'] for s in stems_of(p) for r in read_loci_lines(f'{out}/{s}.loci.jsonl')}
    left = [m['truth_id'] for m in sorted(loci(p), key=lambda m: (m['chrom'], m['pos0'])) if m['truth_id'] not in done]
    open(f'{out}/rest{tag}.ids', 'w').write('\n'.join(left) + '\n')
    print(p, 'done', len(done), 'left', len(left), '->', f'{out}/rest{tag}.ids')


def run_rest(p, tag, k, n):
    ids = open(f'{D}/read_paths_parts/{p}/rest{tag}.ids').read().split()[k::n]
    by = {m['truth_id']: m for m in loci(p)}
    process(p, sorted((by[i] for i in ids), key=lambda m: (m['chrom'], m['pos0'])), f'rest{tag}_{k}of{n}')


COLS = ['truth_id', 'read', 'call', 'has_edit_in_window', 'kind', 'mapq', 'spans_event', 'local_path', 'path_type', 'reverse']


def merge(p, stems=None):
    out = f'{D}/read_paths_parts/{p}'
    stems = stems or stems_of(p)
    info, src = {}, {}                                   # a locus is taken from the first file with its loci line
    for s in stems:
        for r in read_loci_lines(f'{out}/{s}.loci.jsonl'):
            if r['truth_id'] not in info:
                info[r['truth_id']], src[r['truth_id']] = r, s
    rows = collections.defaultdict(list)
    for s in stems:
        for line in open(f'{out}/{s}.rows.tsv'):
            v = line.rstrip('\n').split('\t')
            if len(v) == len(COLS) and line.endswith('\n') and src.get(v[0]) == s:
                rows[v[0]].append(dict(zip(COLS, v)))
    vs = list(csv.DictReader(open(f'{D}/variant_set.tsv'), delimiter='\t'))
    tag = '' if stems[0].startswith(('part', 'rest')) else '_' + stems[0]
    agree, n_cmp, summ = collections.Counter(), 0, []
    with gzip.open(f'{D}/read_paths_{p}{tag}.tsv.gz', 'wt') as w:
        w.write('\t'.join(COLS) + '\n')
        for t in vs:
            for r in rows.get(t['truth_id'], []):
                w.write('\t'.join(r[c] for c in COLS) + '\n')
    scol = ['truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', f'{p}_perfect', f'{p}_status', f'{p}_reason',
            'reason_recomputed', 'reason_agree', 'records', 'low_mapq', 'n_spanning', 'n_ref', 'n_alt', 'n_alt_exact',
            'n_alt_no_edit', 'n_alt_no_edit_spanning', 'top_paths', 'majority_path', 'majority_fraction', 'majority_type',
            'majority_nonref_nodes', 'majority_seq', 'alt_representation']
    with open(f'{D}/read_path_summary_{p}{tag}.tsv', 'w') as w:
        w.write('\t'.join(scol) + '\n')
        for t in vs:
            i = info.get(t['truth_id'])
            if i is None:
                continue
            rs = rows.get(t['truth_id'], [])
            assert len(rs) == i['spanning'], t['truth_id']
            alt = [r for r in rs if r['call'].startswith('alt')]
            good = [r for r in alt if r['has_edit_in_window'] == 'False' and r['spans_event'] == 'True']
            c = collections.Counter(r['local_path'] for r in good)
            top = c.most_common(5)
            maj = top[0][0] if top else ''
            mtype = next((r['path_type'] for r in good if r['local_path'] == maj), '')
            rec = PRIMARY(i)[0]
            old = t[f'{p}_reason']
            ag = '' if not old else str(rec == old)
            if old:
                n_cmp += 1; agree[(old, rec)] += 1
            out_r = dict(t, reason_recomputed=rec, reason_agree=ag, records=i['records'], low_mapq=i['low_mapq'],
                         n_spanning=i['spanning'], n_ref=i['ref_exact'] + i['ref_like'], n_alt=i['alt_exact'] + i['alt_like'],
                         n_alt_exact=i['alt_exact'], n_alt_no_edit=sum(r['has_edit_in_window'] == 'False' for r in alt),
                         n_alt_no_edit_spanning=len(good), top_paths=json.dumps([list(x) for x in top]), majority_path=maj,
                         majority_fraction=round(top[0][1] / len(good), 3) if top else '', majority_type=mtype,
                         majority_nonref_nodes=','.join(NONREF(maj)) if maj else '', majority_seq=SPELL(maj, t) if maj else '',
                         alt_representation=json.dumps(i['alt_representation']))
            w.write('\t'.join(str(out_r.get(k, '')) for k in scol) + '\n')
            summ.append(out_r)
    same = sum(v for (a, b), v in agree.items() if a == b)
    print(p, 'loci', len(info), 'of', len(vs), 'rows', sum(map(len, rows.values())))
    print(p, 'reason agreement', same, '/', n_cmp, f'{same / max(n_cmp, 1):.4f}')
    for (a, b), v in sorted(agree.items(), key=lambda x: -x[1]):
        if a != b:
            print('  disagree', v, a, '->', b)
    for sub in ('all', 'perfect'):
        for kd in ('SNV', 'INDEL'):
            x = [r for r in summ if r['kind2'] == kd and (sub == 'all' or r[f'{p}_perfect'] == 'True')]
            y = [r for r in x if r['majority_path']]
            print(p, sub, kd, 'loci', len(x), 'with a majority no-edit ALT path', len(y), 'median majority_fraction',
                  round(float(np.median([r['majority_fraction'] for r in y])), 3) if y else '',
                  'type', dict(collections.Counter(r['majority_type'] for r in y)), 'spells', dict(collections.Counter(r['majority_seq'] for r in y)))


def NONREF(walk):
    """Node ids of a walk that are not GRCh38 nodes of the d9 reference path (visited once)."""
    ids = [int(x) for x in walk.replace('<', ' ').replace('>', ' ').split()]
    return [str(n) for n in ids if not (_RP.visits[n] == 1 and _RP.chrom[n] >= 0)]


def SPELL(walk, t):
    """ALT / REF / other: the walk's sequence (classify.py oseq) against GRCh38 between its two anchor nodes, with and without
    the truth allele (a no-edit ALT read that also follows a nearby germline branch spells 'other')."""
    w = [(int(x[1:]), '+' if x[0] == '>' else '-') for x in re.findall(r'[<>]\d+', walk)]
    a, z = G['ref'](w[0][0]), G['ref'](w[-1][0])
    if not a or not z or w[0][1] == '-' or w[-1][1] == '-':
        return 'other'
    s = ''.join(G['oseq'](n, o) for n, o in w).upper()
    reg = FA.fetch(t['chrom'], a[1], z[2] + 1).upper()
    off, vref = int(t['vcf_pos']) - 1 - a[1], t['vcf_ref'].upper()
    alt = reg[:off] + t['vcf_alt'].upper() + reg[off + len(vref):] if 0 <= off and off + len(vref) <= len(reg) else None
    return 'ALT' if s == alt else 'REF' if s == reg else 'other'


if __name__ == '__main__':
    mode, p = sys.argv[1], sys.argv[2]
    MOD, PRIMARY = load(p)
    if mode == 'run':
        run(p, int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]) if len(sys.argv) > 5 else None)
    elif mode == 'rest-list':
        rest_list(p, sys.argv[3])
    elif mode == 'rest':
        run_rest(p, sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    else:
        _RP = MOD.ReferencePath(MOD.REFPATH)
        G = {}
        exec(open(CLASSIFY).read().split("fa = pysam.FastaFile(FASTA)\nchroms")[0], G)   # as s1_loci.py
        FA = MOD.pysam.FastaFile(G['FASTA'])
        merge(p, sys.argv[3:] or None)

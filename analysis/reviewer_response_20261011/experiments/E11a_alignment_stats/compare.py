"""E11a tables: the same reads on the GRCh38 BAM, the GRCh38-only graph (linear), HPRC v1.1 d9 and HPRC v2.1 d46.

Two read samples, one analysis each:
  G  genome-wide name-hash sample (long reads; scan_gam.py full scans of the three graph GAMs):
     zlib.crc32(name) % 100 < 10; every sampled read is in every graph arm's scan (mapped or not).
  W  window sample (window_scan.py): reads placed in the fixed GRCh38 windows (windows.tsv; long reads also
     windows_extra.tsv) by at least one arm. Short reads: crc32(name) % 10 < 1 within the windows, all four arms
     from window scans + the GAMs' unmapped heads. Long reads: the G sample's graph-arm rows + the BAM window scan
     (same name rule).
A read is (name, mate). Per read and arm (derived()):
  unmapped     no path (GAM) / flag 0x4 (BAM)
  perfect      mapped, full length, edit-free: no mismatch, no indel, no soft clip (BAM: no S/H either)
  mismatch / indel / softclip   >= 1 mismatched / inserted-or-deleted / soft-clipped base
  identity     match / (match + mismatch + inserted + deleted bases), soft clips excluded
  placement    GRCh38 contig + span (GAM: anchored graph nodes, node granularity; BAM: primary record)
  concordant   both mapped, same contig, spans overlap or lie within 1 kb
W status per arm: `in` = mapped and placed overlapping the read's window (+/- flank: 1 kb short, 25 kb long);
`elsewhere` = mapped elsewhere (window scans: not seen near the window and not unmapped); `unmapped`.
Strata (one anchor placement per read, the same for every arm): G: linear, else d9, else d46; W: the GRCh38 BAM
when it places the read in the window, else linear, d9, d46. BED strata take >= 50 % of the anchor span:
  SMaHT v2 easy / difficult / extreme (disjoint, tile chr1-22,X,Y; the one covering most of the span)
  segdup (GRCh38_segdups.bed); HG008T somatic benchmark v0.2 all.bed (HG008T only)
Run: python compare.py [SET ...]  -> results/<SET>.json, results/tables.md
"""
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = Path('/scratch/jshen/data/reviewer_experiments_20261013/E11a_alignment_stats')
SCANS, WINS = DATA / 'scans', DATA / 'windows'
RES = HERE / 'results'
ARMS = ['grch38_bam', 'linear', 'd9', 'd46']
GRAPH = ['linear', 'd9', 'd46']
ARM_LABEL = {'grch38_bam': 'GRCh38 BAM', 'linear': 'GRCh38-only graph', 'd9': 'HPRC v1.1 d9', 'd46': 'HPRC v2.1 d46'}
SETS = ['HG008T_PacBio', 'HG008T_ONT', 'HG008T_Illumina', 'COLO829T_fiberseq', 'COLO829T_ONT', 'COLO829T_Illumina']
TOTAL_READS = {  # records of every GAM of the set (gam stats / earlier counts; the G scans re-count long reads)
    'HG008T_Illumina': 2974694440, 'COLO829T_Illumina': 2479603422}
BEDS = {
    'smaht_easy': '/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_easy_v2.bed',
    'smaht_difficult': '/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_difficult_v2.bed',
    'smaht_extreme': '/scratch/qfu/reference/hg38_DAC_stratification/SMaHT_extreme_v2.bed',
    'segdup': '/scratch/qfu/real_mixture/repeatmasker/GRCh38_segdups.bed',
    'hg008t_benchmark': '/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed',
}
MAPQ_BINS = [(0, 0), (1, 9), (10, 29), (30, 59), (60, 60)]
NEAR = 1000
PAIRS = [('linear', 'd9'), ('linear', 'd46'), ('d9', 'd46'), ('grch38_bam', 'linear'), ('grch38_bam', 'd9'),
         ('grch38_bam', 'd46')]


# --- BED coverage -------------------------------------------------------------------------------

class Bed:
    def __init__(self, path):
        per = {}
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                c, s, e = line.split('\t', 3)[:3]
                per.setdefault(c, []).append((int(s), int(e)))
        self.iv = {}
        for c, ivs in per.items():
            ivs.sort()
            ms, me = [], []
            for s, e in ivs:
                if ms and s <= me[-1]:
                    me[-1] = max(me[-1], e)
                else:
                    ms.append(s)
                    me.append(e)
            ms, me = np.array(ms, dtype=np.int64), np.array(me, dtype=np.int64)
            self.iv[c] = (ms, me, np.cumsum(me - ms))
        self.total = int(sum(v[2][-1] for v in self.iv.values()))

    def below(self, c, x):
        ms, me, cum = self.iv[c]
        k = np.searchsorted(ms, x, side='left')
        out = np.where(k > 0, cum[np.maximum(k - 1, 0)], 0)
        over = np.where(k > 0, np.maximum(0, me[np.maximum(k - 1, 0)] - x), 0)
        return out - over

    def fraction(self, chrom, start, end):
        frac = np.zeros(len(chrom))
        for c in np.unique(chrom):
            if c not in self.iv:
                continue
            m = chrom == c
            s, e = start[m], end[m]
            frac[m] = (self.below(c, e) - self.below(c, s)) / np.maximum(e - s, 1)
        return frac


# --- per-read metrics ---------------------------------------------------------------------------

def pct(a, b):
    return 100.0 * float(a) / float(b) if b else float('nan')


def derived(r):
    mapped = r['mapped'] == 1
    ind = (r['ins'] > 0) | (r['dl'] > 0)
    clip = (r['sc_l'].astype(np.int64) + r['sc_r']) > 0
    perfect = mapped & (r['mm'] == 0) & ~ind & ~clip & (r['cx'] == 0)
    aligned = r['match'].astype(np.int64) + r['mm'] + r['ins'] + r['dl']
    ident = np.where(aligned > 0, r['match'] / np.maximum(aligned, 1), np.nan)
    return dict(mapped=mapped, perfect=perfect, mismatch=mapped & (r['mm'] > 0), indel=mapped & ind,
                softclip=mapped & clip, identity=ident, aligned=aligned)


def arm_metrics(rows):
    n = len(rows)
    d = derived(rows)
    m = d['mapped']
    nm = int(m.sum())
    out = dict(reads=n, mapped=nm, unmapped_pct=pct(n - nm, n), perfect_pct=pct(d['perfect'].sum(), n),
               perfect_of_mapped_pct=pct(d['perfect'].sum(), nm), mismatch_pct=pct(d['mismatch'].sum(), nm),
               indel_pct=pct(d['indel'].sum(), nm), softclip_pct=pct(d['softclip'].sum(), nm))
    mq = rows['mapq'][m]
    for lo, hi in MAPQ_BINS:
        out[f'mapq_{lo}_{hi}_pct'] = pct(((mq >= lo) & (mq <= hi)).sum(), nm)
    out['mapq_mean'] = float(mq.mean()) if nm else float('nan')
    al = int(d['aligned'][m].sum())
    out['identity_median'] = float(np.nanmedian(d['identity'][m])) if nm else float('nan')
    out['mm_per_kb'] = 1000.0 * float(rows['mm'][m].sum()) / al if al else float('nan')
    out['indel_ev_per_kb'] = (1000.0 * float(rows['ins_ev'][m].astype(np.int64).sum() + rows['del_ev'][m].sum()) / al
                              if al else float('nan'))
    rl = int(rows['read_len'][m].astype(np.int64).sum())
    out['softclip_bp_pct'] = pct((rows['sc_l'][m].astype(np.int64) + rows['sc_r'][m]).sum(), rl)
    off = rows['off_visits'][m].astype(np.int64)
    out['offref_read_pct'] = pct((off > 0).sum(), nm)
    out['offref_visit_pct'] = pct(off.sum(), off.sum() + rows['ref_visits'][m].astype(np.int64).sum())
    out['offref_visits_edit_free_pct'] = pct(off.sum() - rows['off_edited'][m].astype(np.int64).sum(), off.sum())
    out['multi_contig_pct'] = pct((rows['nchrom'][m] > 1).sum(), nm)
    if 'supp' in rows.dtype.names:
        out['supplementary_pct'] = pct((rows['supp'][m] > 0).sum(), nm)
    return out


def concordant(ca, sa, ea, cb, sb, eb):
    return (ca == cb) & (ca != '') & (sa <= eb + NEAR) & (sb <= ea + NEAR)


def pair_metrics(ra, ca, rb, cb):
    """B relative to A for aligned row arrays (the same reads)."""
    da, db = derived(ra), derived(rb)
    n = len(ra)
    both = da['mapped'] & db['mapped']
    nb = int(both.sum())
    qa, qb = ra['mapq'].astype(int), rb['mapq'].astype(int)
    conc = concordant(ca, ra['gstart'], ra['gend'], cb, rb['gstart'], rb['gend']) & both
    hq = both & (qa >= 30) & (qb >= 30)
    out = dict(reads=n, mapped_both=nb,
               unmapped_A_only=int((~da['mapped'] & db['mapped']).sum()),
               unmapped_B_only=int((da['mapped'] & ~db['mapped']).sum()),
               mapq_up_pct=pct((both & (qb > qa)).sum(), nb), mapq_down_pct=pct((both & (qb < qa)).sum(), nb),
               mapq30_gained_pct=pct((both & (qa < 30) & (qb >= 30)).sum(), nb),
               mapq30_lost_pct=pct((both & (qa >= 30) & (qb < 30)).sum(), nb),
               mapq0_A_pct=pct((both & (qa == 0)).sum(), nb), mapq0_B_pct=pct((both & (qb == 0)).sum(), nb),
               perfect_A_pct=pct(da['perfect'].sum(), n), perfect_B_pct=pct(db['perfect'].sum(), n),
               perfect_gained_pct=pct((~da['perfect'] & db['perfect']).sum(), n),
               perfect_lost_pct=pct((da['perfect'] & ~db['perfect']).sum(), n),
               softclip_A_pct=pct(da['softclip'].sum(), n), softclip_B_pct=pct(db['softclip'].sum(), n),
               concordant_pct=pct(conc.sum(), nb), mapq30_both=int(hq.sum()),
               concordant_mapq30_both_pct=pct((conc & hq).sum(), hq.sum()),
               discordant_mapq30_both=int((hq & ~conc).sum()))
    gained = ~da['perfect'] & db['perfect']
    out['gained_perfect'] = int(gained.sum())
    if gained.any():
        g_off = rb['off_visits'][gained] > 0
        gb = gained & both
        out.update(gained_perfect_offref_pct=pct(g_off.sum(), gained.sum()),
                   gained_perfect_offref_visits_median=float(np.median(rb['off_visits'][gained][g_off])) if g_off.any() else 0.0,
                   gained_perfect_A_mismatch_only_pct=pct((da['mismatch'] & ~da['indel'] & ~da['softclip'])[gained].sum(), gained.sum()),
                   gained_perfect_A_indel_pct=pct(da['indel'][gained].sum(), gained.sum()),
                   gained_perfect_A_softclip_pct=pct(da['softclip'][gained].sum(), gained.sum()),
                   gained_perfect_A_unmapped_pct=pct((~da['mapped'])[gained].sum(), gained.sum()),
                   gained_perfect_A_mm_median=float(np.median(ra['mm'][gained])),
                   gained_perfect_concordant_pct=pct((conc & gained).sum(), gb.sum()),
                   gained_perfect_mapq_up_pct=pct((gb & (qb > qa)).sum(), gb.sum()),
                   gained_perfect_mapq_down_pct=pct((gb & (qb < qa)).sum(), gb.sum()))
    lost = da['perfect'] & ~db['perfect']
    out['lost_perfect'] = int(lost.sum())
    if lost.any():
        out['lost_perfect_B_offref_pct'] = pct((rb['off_visits'][lost] > 0).sum(), lost.sum())
        out['lost_perfect_B_unmapped_pct'] = pct((~db['mapped'])[lost].sum(), lost.sum())
    return out


# --- loading ------------------------------------------------------------------------------------

def keys_of(rows):
    return np.char.add(np.char.add(rows['name'], b'\t'), rows['mate'].astype('S1'))


def load(path):
    if not Path(str(path) + '.npz').exists():
        return None
    d = np.load(str(path) + '.npz')
    rows = d['rows']
    contigs = json.loads(str(d['contigs']))
    chrom = np.array([contigs[c] if c >= 0 else '' for c in rows['chrom']], dtype='U32')
    out = dict(rows=rows, key=keys_of(rows), chrom=chrom, summary=json.load(open(str(path) + '.summary.json')))
    if 'unmapped' in d.files:
        out['unmapped'] = d['unmapped']
    return out


def index_of(keys, sub):
    """positions of `sub` in `keys` (-1 when absent)."""
    order = np.argsort(keys)
    pos = np.searchsorted(keys[order], sub)
    pos = np.minimum(pos, len(keys) - 1)
    hit = keys[order][pos] == sub
    return np.where(hit, order[pos], -1)


def read_windows(extra):
    out = []
    for name in ['windows.tsv'] + (['windows_extra.tsv'] if extra else []):
        with open(HERE / name) as f:
            next(f)
            for line in f:
                w, c, s, e, k = line.rstrip('\n').split('\t')[:5]
                out.append((w, c, int(s), int(e), k))
    return out


def window_of(chrom, gs, ge, wins, flank):
    """index of the window the placement overlaps (+/- flank), -1 if none."""
    w = np.full(len(chrom), -1, dtype=np.int64)
    for i, (_, c, s, e, _) in enumerate(wins):
        m = (chrom == c) & (gs < e + flank) & (ge > s - flank) & (w < 0)
        w[m] = i
    return w


def strata_masks(s, chrom, gs, ge, beds):
    fr = {k: b.fraction(chrom, gs, ge) for k, b in beds.items()}
    sm = np.stack([fr['smaht_easy'], fr['smaht_difficult'], fr['smaht_extreme']])
    best = np.argmax(sm, axis=0)
    placed = (chrom != '') & (sm.max(axis=0) > 0)
    out = {f'SMaHT {n}': placed & (best == i) for i, n in enumerate(['easy', 'difficult', 'extreme'])}
    out['segdup'] = fr['segdup'] >= 0.5
    out['not segdup'] = placed & (fr['segdup'] < 0.5)
    if s.startswith('HG008T'):
        out['in HG008T benchmark BED'] = fr['hg008t_benchmark'] >= 0.5
        out['outside HG008T benchmark BED'] = placed & (fr['hg008t_benchmark'] < 0.5)
    return out


# --- analysis G (genome-wide name-hash sample, graph arms) --------------------------------------

def analysis_g(s, beds):
    arms = {a: load(SCANS / f'{s}.{a}') for a in GRAPH}
    arms = {a: v for a, v in arms.items() if v is not None and v['summary'].get('complete')}
    if len(arms) < 2:
        return None
    out = dict(sample='G', arms=list(arms), whole_file={}, per_arm={}, pairs={}, strata={})
    for a, A in arms.items():
        sm = A['summary']
        h = np.array(sm['mapq_hist_with_path'])
        u, cnt = np.unique(A['key'], return_counts=True)
        out['whole_file'][a] = dict(records=sm['records'], unmapped=sm['records'] - sm['with_path'],
                                    unmapped_pct=pct(sm['records'] - sm['with_path'], sm['records']),
                                    mapq60_of_mapped_pct=pct(h[60], h.sum()), mapq0_of_mapped_pct=pct(h[0], h.sum()),
                                    mapq_le5_of_mapped_pct=pct(h[:6].sum(), h.sum()),
                                    sampled=sm['selected'], duplicate_keys=int((cnt > 1).sum()),
                                    gam=sm['gam'], gam_bytes=sm['gam_bytes'])
    # union of the sampled reads; a read absent from an arm's GAM (giraffe did not write it) counts as unmapped there
    keys = np.unique(np.concatenate([A['key'] for A in arms.values()]))
    out['common_reads'] = int(len(keys))
    out['read_presence'] = {a: int(len(A['key'])) for a, A in arms.items()}
    idx = {a: index_of(A['key'], keys) for a, A in arms.items()}
    R, C = {}, {}
    for a in arms:
        miss = idx[a] < 0
        R[a] = arms[a]['rows'][np.maximum(idx[a], 0)].copy()
        C[a] = arms[a]['chrom'][np.maximum(idx[a], 0)].copy()
        if miss.any():
            R[a][miss] = np.zeros(1, dtype=R[a].dtype)
            C[a][miss] = ''
        out['whole_file'][a]['sampled_reads_absent'] = int(miss.sum())
    for a in arms:
        out['per_arm'][a] = arm_metrics(R[a])
    for a, b in PAIRS:
        if a in arms and b in arms:
            out['pairs'][f'{b}_vs_{a}'] = pair_metrics(R[a], C[a], R[b], C[b])
    n = len(keys)
    ach, ast, aen = np.full(n, '', dtype='U32'), np.zeros(n, np.int64), np.zeros(n, np.int64)
    src = np.full(n, '', dtype='U8')
    for a in GRAPH:
        if a in arms:
            take = (ach == '') & (R[a]['mapped'] == 1) & (C[a] != '')
            ach[take], ast[take], aen[take], src[take] = C[a][take], R[a]['gstart'][take], R[a]['gend'][take], a
    out['anchor_source'] = {a: int((src == a).sum()) for a in GRAPH}
    out['anchor_source']['none'] = int((src == '').sum())
    for name, m in strata_masks(s, ach, ast, aen, beds).items():
        e = dict(reads=int(m.sum()), per_arm={}, pairs={})
        if m.any():
            for a in arms:
                e['per_arm'][a] = arm_metrics(R[a][m])
            for a, b in PAIRS:
                if a in arms and b in arms:
                    e['pairs'][f'{b}_vs_{a}'] = pair_metrics(R[a][m], C[a][m], R[b][m], C[b][m])
        out['strata'][name] = e
    return out


# --- analysis W (window sample, all arms) ---------------------------------------------------------

def analysis_w(s, beds):
    long = not s.endswith('Illumina')
    wins = read_windows(extra=long)
    flank = 25_000 if long else 1_000
    arms = {}
    for a in ARMS:
        if long and a in GRAPH:
            A = load(SCANS / f'{s}.{a}')
            if A is not None and A['summary'].get('complete'):
                A['scope'] = 'genome'
                arms[a] = A
        else:
            A = load(WINS / f'{s}.{a}')
            if A is not None and A['summary'].get('complete'):
                A['scope'] = 'window'
                arms[a] = A
    if len(arms) < 2:
        return None
    # per arm: window of each row's placement
    for a, A in arms.items():
        r = A['rows']
        A['win'] = np.where(r['mapped'] == 1, window_of(A['chrom'], r['gstart'], r['gend'], wins, flank), -1)
    universe = np.unique(np.concatenate([A['key'][A['win'] >= 0] for A in arms.values()]))
    n = len(universe)
    out = dict(sample='W', arms=list(arms), scope={a: A['scope'] for a, A in arms.items()}, windows=len(wins),
               flank=flank, reads=n, status={}, per_arm={}, pairs={}, strata={}, whole_file={})
    st, R, C, W = {}, {}, {}, {}
    for a, A in arms.items():
        i = index_of(A['key'], universe)
        found = i >= 0
        rows = A['rows'][np.maximum(i, 0)]
        win = np.where(found, A['win'][np.maximum(i, 0)], -1)
        mapped = found & (rows['mapped'] == 1)
        unm = found & (rows['mapped'] == 0)
        if 'unmapped' in A:
            unm |= ~found & np.isin(universe, A['unmapped'])
        if A['scope'] == 'genome':  # a full scan has every read giraffe wrote: absent = not written = unmapped
            unm |= ~found
        status = np.where(win >= 0, 'in', np.where(unm, 'unmapped', 'elsewhere'))
        st[a], R[a], C[a], W[a] = status, rows, np.where(found, A['chrom'][np.maximum(i, 0)], ''), win
        out['status'][a] = {k: int((status == k).sum()) for k in ('in', 'elsewhere', 'unmapped')}
        out['status'][a].update({f'{k}_pct': pct((status == k).sum(), n) for k in ('in', 'elsewhere', 'unmapped')})
        out['per_arm'][a] = arm_metrics(rows[status == 'in'])
        sm = A['summary']
        if A['scope'] == 'window':
            hd = sm.get('head', [])
            out['whole_file'][a] = dict(head_records=sum(h['records'] for h in hd),
                                        unmapped_records=sum(h['unmapped'] for h in hd))
            if a in GRAPH and s in TOTAL_READS:
                out['whole_file'][a]['unmapped_pct_of_all_reads'] = pct(out['whole_file'][a]['unmapped_records'], TOTAL_READS[s])
    # the read's window: first arm that places it in a window (BAM first)
    wid = np.full(n, -1, dtype=np.int64)
    ach, ast, aen = np.full(n, '', dtype='U32'), np.zeros(n, np.int64), np.zeros(n, np.int64)
    for a in ARMS:
        if a in arms:
            take = (wid < 0) & (st[a] == 'in')
            wid[take] = W[a][take]
            ach[take], ast[take], aen[take] = C[a][take], R[a]['gstart'][take], R[a]['gend'][take]
    kinds = np.array([wins[i][4] for i in wid])
    out['reads_by_window_kind'] = {k: int((kinds == k).sum()) for k in sorted(set(kinds.tolist()))}
    for a, b in PAIRS:
        if a in arms and b in arms:
            both_in = (st[a] == 'in') & (st[b] == 'in')
            p = pair_metrics(R[a][both_in], C[a][both_in], R[b][both_in], C[b][both_in])
            p.update(A_in=int((st[a] == 'in').sum()), B_in=int((st[b] == 'in').sum()), both_in=int(both_in.sum()),
                     A_in_B_elsewhere=int(((st[a] == 'in') & (st[b] == 'elsewhere')).sum()),
                     A_in_B_unmapped=int(((st[a] == 'in') & (st[b] == 'unmapped')).sum()),
                     B_in_A_elsewhere=int(((st[b] == 'in') & (st[a] == 'elsewhere')).sum()),
                     B_in_A_unmapped=int(((st[b] == 'in') & (st[a] == 'unmapped')).sum()))
            out['pairs'][f'{b}_vs_{a}'] = p
    masks = {f'windows: {k}': kinds == k for k in sorted(set(kinds.tolist()))}
    masks.update(strata_masks(s, ach, ast, aen, beds))
    for name, m in masks.items():
        e = dict(reads=int(m.sum()), status={}, per_arm={}, pairs={})
        if m.any():
            for a in arms:
                e['status'][a] = {f'{k}_pct': pct(((st[a] == k) & m).sum(), m.sum()) for k in ('in', 'elsewhere', 'unmapped')}
                e['per_arm'][a] = arm_metrics(R[a][m & (st[a] == 'in')])
            for a, b in PAIRS:
                if a in arms and b in arms:
                    both_in = m & (st[a] == 'in') & (st[b] == 'in')
                    p = pair_metrics(R[a][both_in], C[a][both_in], R[b][both_in], C[b][both_in])
                    p.update(A_in_B_elsewhere_pct=pct((m & (st[a] == 'in') & (st[b] == 'elsewhere')).sum(), (m & (st[a] == 'in')).sum()),
                             A_in_B_unmapped_pct=pct((m & (st[a] == 'in') & (st[b] == 'unmapped')).sum(), (m & (st[a] == 'in')).sum()))
                    e['pairs'][f'{b}_vs_{a}'] = p
        out['strata'][name] = e
    return out


# --- markdown -----------------------------------------------------------------------------------

def f(x, nd=2):
    if x is None:
        return 'n/a'
    if isinstance(x, (int, np.integer)) and not isinstance(x, bool):
        return f'{x:,}'
    if isinstance(x, float) and np.isnan(x):
        return 'n/a'
    return f'{x:.{nd}f}'


METRIC_HEAD = ('| arm | reads | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 1-9 % | MAPQ 10-29 % | '
               'MAPQ 30-59 % | MAPQ 60 % | median identity | mismatch bp / kb | indel events / kb | soft-clipped bp % | '
               'reads on off-ref nodes % | off-ref share of node visits % | off-ref visits edit-free % | multi-contig % |')
METRIC_SEP = '|---|' + '---:|' * 18


def metric_row(label, m, a):
    off = a in ('d9', 'd46')
    return (f"| {label} | {f(m['reads'])} | {f(m['perfect_pct'])} | {f(m['mismatch_pct'])} | {f(m['indel_pct'])} | "
            f"{f(m['softclip_pct'])} | {f(m['mapq_0_0_pct'])} | {f(m['mapq_1_9_pct'])} | {f(m['mapq_10_29_pct'])} | "
            f"{f(m['mapq_30_59_pct'])} | {f(m['mapq_60_60_pct'])} | {f(m['identity_median'], 4)} | {f(m['mm_per_kb'])} | "
            f"{f(m['indel_ev_per_kb'])} | {f(m['softclip_bp_pct'])} | {f(m['offref_read_pct']) if off else '-'} | "
            f"{f(m['offref_visit_pct']) if off else '-'} | {f(m['offref_visits_edit_free_pct']) if off else '-'} | "
            f"{f(m['multi_contig_pct'])} |")


PAIR_HEAD = ('| B vs A | reads | both mapped | MAPQ up % | MAPQ down % | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect A % | '
             'perfect B % | perfect gained % | perfect lost % | soft clip A % | soft clip B % | placement concordant % | '
             'concordant if MAPQ>=30 in both % | discordant with MAPQ>=30 in both |')
PAIR_SEP = '|---|' + '---:|' * 15


def pair_row(k, p):
    b, a = k.split('_vs_')
    return (f"| {ARM_LABEL[b]} vs {ARM_LABEL[a]} | {f(p['reads'])} | {f(p['mapped_both'])} | {f(p['mapq_up_pct'])} | "
            f"{f(p['mapq_down_pct'])} | {f(p['mapq30_gained_pct'])} | {f(p['mapq30_lost_pct'])} | {f(p['perfect_A_pct'])} | "
            f"{f(p['perfect_B_pct'])} | {f(p['perfect_gained_pct'])} | {f(p['perfect_lost_pct'])} | {f(p['softclip_A_pct'])} | "
            f"{f(p['softclip_B_pct'])} | {f(p['concordant_pct'])} | {f(p['concordant_mapq30_both_pct'])} | "
            f"{f(p['discordant_mapq30_both'])} |")


GAIN_HEAD = ('| B vs A | reads perfect in B, not in A | of these on off-ref nodes in B % | median off-ref visits | '
             'A had mismatches only % | A had an indel % | A had a soft clip % | A unmapped % | median A mismatched bp | '
             'placement concordant % | MAPQ up % | MAPQ down % | reads perfect in A, not in B | of these on off-ref in B % |')
GAIN_SEP = '|---|' + '---:|' * 13


def gain_row(k, p):
    b, a = k.split('_vs_')
    if not p.get('gained_perfect'):
        return f"| {ARM_LABEL[b]} vs {ARM_LABEL[a]} | 0 |" + ' |' * 12
    return (f"| {ARM_LABEL[b]} vs {ARM_LABEL[a]} | {f(p['gained_perfect'])} | {f(p['gained_perfect_offref_pct'])} | "
            f"{f(p['gained_perfect_offref_visits_median'], 1)} | {f(p['gained_perfect_A_mismatch_only_pct'])} | "
            f"{f(p['gained_perfect_A_indel_pct'])} | {f(p['gained_perfect_A_softclip_pct'])} | "
            f"{f(p['gained_perfect_A_unmapped_pct'])} | {f(p['gained_perfect_A_mm_median'], 1)} | "
            f"{f(p['gained_perfect_concordant_pct'])} | {f(p['gained_perfect_mapq_up_pct'])} | "
            f"{f(p['gained_perfect_mapq_down_pct'])} | {f(p['lost_perfect'])} | {f(p.get('lost_perfect_B_offref_pct'))} |")


def write_md(results):
    L = ['# E11a tables: the same reads on four references', '',
         'Written by compare.py; definitions in README.md. Percentages: perfect of all reads in the row; mismatch, '
         'indel, soft clip and MAPQ bins of mapped reads.', '']
    for s, res in results:
        L += [f'## {s}', '']
        g = res.get('G')
        if g:
            L += ['### G. Genome-wide name-hash sample (crc32(name) % 100 < 10), graph arms', '',
                  'Whole GAM (every record, exact):', '',
                  '| arm | records | unmapped | unmapped % | MAPQ 0 % of mapped | MAPQ <= 5 % of mapped | MAPQ 60 % of mapped | sampled reads |',
                  '|---|---:|---:|---:|---:|---:|---:|---:|']
            for a in g['arms']:
                w = g['whole_file'][a]
                L.append(f"| {ARM_LABEL[a]} | {f(w['records'])} | {f(w['unmapped'])} | {f(w['unmapped_pct'])} | "
                         f"{f(w['mapq0_of_mapped_pct'])} | {f(w['mapq_le5_of_mapped_pct'])} | {f(w['mapq60_of_mapped_pct'])} | {f(w['sampled'])} |")
            L += ['', f"Sampled reads (union over graph arms): {f(g['common_reads'])}; absent from an arm (counted unmapped there): " + ', '.join(f"{ARM_LABEL[a]} {f(g['whole_file'][a]['sampled_reads_absent'])}" for a in g['arms']) + '.', '', METRIC_HEAD, METRIC_SEP]
            for a in g['arms']:
                m = g['per_arm'][a]
                L.append(metric_row(f"{ARM_LABEL[a]} (unmapped {f(m['unmapped_pct'])} %)", m, a))
            L += ['', 'Same reads, B relative to A:', '', PAIR_HEAD, PAIR_SEP]
            L += [pair_row(k, p) for k, p in g['pairs'].items()]
            L += ['', 'Reads that became perfect / stopped being perfect:', '', GAIN_HEAD, GAIN_SEP]
            L += [gain_row(k, p) for k, p in g['pairs'].items()]
            L += ['', 'By region (anchor: ' + ', '.join(f'{k} {f(v)}' for k, v in g['anchor_source'].items()) + '):', '',
                  '| stratum | reads | arm | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity | reads on off-ref nodes % |',
                  '|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
            for name, e in g['strata'].items():
                for a in g['arms']:
                    m = e['per_arm'].get(a)
                    if m:
                        L.append(f"| {name} | {f(e['reads'])} | {ARM_LABEL[a]} | {f(m['unmapped_pct'])} | {f(m['perfect_pct'])} | "
                                 f"{f(m['mismatch_pct'])} | {f(m['indel_pct'])} | {f(m['softclip_pct'])} | {f(m['mapq_0_0_pct'])} | "
                                 f"{f(m['mapq_60_60_pct'])} | {f(m['identity_median'], 4)} | "
                                 f"{f(m['offref_read_pct']) if a in ('d9', 'd46') else '-'} |")
            L += ['', '| stratum | B vs A | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | placement concordant % | discordant with MAPQ>=30 in both | became perfect: on off-ref nodes % |',
                  '|---|---|---:|---:|---:|---:|---:|---:|---:|']
            for name, e in g['strata'].items():
                for k, p in e['pairs'].items():
                    b, a = k.split('_vs_')
                    L.append(f"| {name} | {ARM_LABEL[b]} vs {ARM_LABEL[a]} | {f(p['mapq30_gained_pct'])} | {f(p['mapq30_lost_pct'])} | "
                             f"{f(p['perfect_gained_pct'])} | {f(p['perfect_lost_pct'])} | {f(p['concordant_pct'])} | "
                             f"{f(p['discordant_mapq30_both'])} | {f(p.get('gained_perfect_offref_pct'))} |")
            L.append('')
        w = res.get('W')
        if w:
            L += [f"### W. Window sample ({w['windows']} GRCh38 windows, flank {f(w['flank'])} bp; reads placed in a window by at least one arm: {f(w['reads'])}; "
                  + ', '.join(f'{k} {f(v)}' for k, v in w['reads_by_window_kind'].items()) + ')', '',
                  '| arm | scope | in the window % | placed elsewhere % | unmapped % |', '|---|---|---:|---:|---:|']
            for a in w['arms']:
                x = w['status'][a]
                L.append(f"| {ARM_LABEL[a]} | {w['scope'][a]} | {f(x['in_pct'])} | {f(x['elsewhere_pct'])} | {f(x['unmapped_pct'])} |")
            if w['whole_file']:
                L += ['', 'Unmapped records of the whole file (GAM head / BAM unplaced): ' + '; '.join(
                    f"{ARM_LABEL[a]} {f(v['unmapped_records'])}" + (f" ({f(v['unmapped_pct_of_all_reads'])} % of all reads)" if 'unmapped_pct_of_all_reads' in v else '')
                    for a, v in w['whole_file'].items())]
            L += ['', 'Reads placed in the window by the arm:', '', METRIC_HEAD, METRIC_SEP]
            for a in w['arms']:
                L.append(metric_row(ARM_LABEL[a], w['per_arm'][a], a))
            L += ['', 'Same reads, both arms place them in the window (B relative to A); last columns: reads A places in the window that B places elsewhere / leaves unmapped:', '',
                  PAIR_HEAD[:-1] + ' A in, B elsewhere | A in, B unmapped | B in, A elsewhere | B in, A unmapped |', PAIR_SEP + '---:|' * 4]
            for k, p in w['pairs'].items():
                L.append(pair_row(k, p) + f" {f(p['A_in_B_elsewhere'])} | {f(p['A_in_B_unmapped'])} | {f(p['B_in_A_elsewhere'])} | {f(p['B_in_A_unmapped'])} |")
            L += ['', GAIN_HEAD, GAIN_SEP]
            L += [gain_row(k, p) for k, p in w['pairs'].items()]
            L += ['', 'By window kind and region:', '',
                  '| stratum | reads | arm | in window % | elsewhere % | unmapped % | perfect % | mismatch % | indel % | soft clip % | MAPQ 0 % | MAPQ 60 % | median identity |',
                  '|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
            for name, e in w['strata'].items():
                for a in w['arms']:
                    m, x = e['per_arm'].get(a), e['status'].get(a)
                    if m:
                        L.append(f"| {name} | {f(e['reads'])} | {ARM_LABEL[a]} | {f(x['in_pct'])} | {f(x['elsewhere_pct'])} | "
                                 f"{f(x['unmapped_pct'])} | {f(m['perfect_pct'])} | {f(m['mismatch_pct'])} | {f(m['indel_pct'])} | "
                                 f"{f(m['softclip_pct'])} | {f(m['mapq_0_0_pct'])} | {f(m['mapq_60_60_pct'])} | {f(m['identity_median'], 4)} |")
            L += ['', '| stratum | B vs A | both in window | MAPQ>=30 gained % | MAPQ>=30 lost % | perfect gained % | perfect lost % | concordant (1 kb) % | A in, B elsewhere % | A in, B unmapped % |',
                  '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
            for name, e in w['strata'].items():
                for k, p in e['pairs'].items():
                    b, a = k.split('_vs_')
                    L.append(f"| {name} | {ARM_LABEL[b]} vs {ARM_LABEL[a]} | {f(p['reads'])} | {f(p['mapq30_gained_pct'])} | "
                             f"{f(p['mapq30_lost_pct'])} | {f(p['perfect_gained_pct'])} | {f(p['perfect_lost_pct'])} | "
                             f"{f(p['concordant_pct'])} | {f(p['A_in_B_elsewhere_pct'])} | {f(p['A_in_B_unmapped_pct'])} |")
            L.append('')
    (RES / 'tables.md').write_text('\n'.join(L) + '\n')


def main():
    sets = sys.argv[1:] or SETS
    RES.mkdir(exist_ok=True)
    beds = {k: Bed(v) for k, v in BEDS.items()}
    for s in sets:
        res = {}
        g = analysis_g(s, beds) if not s.endswith('Illumina') else None
        if g:
            res['G'] = g
        w = analysis_w(s, beds)
        if w:
            res['W'] = w
        if res:
            json.dump(res, open(RES / f'{s}.json', 'w'), indent=1, default=float)
        print(s, {k: (v.get('common_reads') or v.get('reads')) for k, v in res.items()} or 'nothing complete yet')
    allres = [(s, json.load(open(RES / f'{s}.json'))) for s in SETS if (RES / f'{s}.json').exists()]
    write_md(allres)


if __name__ == '__main__':
    main()

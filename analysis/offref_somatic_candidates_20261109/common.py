"""Shared helpers of the off-reference somatic candidate scripts (s2-s5).

Coordinates are 0-based, half-open. A query is an s1 context: alt_q (A1 read consensus, site [a, b)) or ref_q (the same
read haplotype with the graph REF at the site, site from ref_site()).
"""
import hashlib
import os

DATA = '/scratch/jshen/data/pansoma_net_v2_runs/offref_somatic_candidates_20261109'
TARGETS = '/scratch/jshen/Github/Pansoma/tmp/offref_somatic_candidates_20261109/targets.tsv'
GRCH38 = '/scratch/jshen/data/reference/GRCh38_no_alt_analysis_set.fasta'
MIN_LEN = 30          # minimap2 -x sr needs queries of >= 30 bp
MIN_FLANK = 10        # read-consensus bases on each side of the site
MIN_RUN_FLANK = 5     # bases beyond the event window (site, its equivalent placements and the tandem run touching it)
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(s):
    return s.translate(COMP)[::-1]


def qid(seq, site):
    return 'q' + hashlib.sha1(f'{seq}|{site}'.encode()).hexdigest()[:16]


def site_of(r):
    a, b = map(int, r['site'].split('-'))
    return a, b


def ref_site(r):
    """site interval of ref_q: same left flank; the right flank has the same length as alt_q's."""
    a, b = site_of(r)
    right = len(r['alt_q']) - b
    return f'{a}-{len(r["ref_q"]) - right}'


# s1 tables per tumor: the three d46 sets x SNV / INDEL x chr1-22; COLO829T also the new off-reference SNV keys of the
# Illumina / Fiber-seq d46_snvaf05 sets (SNV AF floor 0.05; their INDEL tables and their AF >= 0.08 SNV keys are the
# d46 ones)
S1_SETS = {'HG008T': [(f'HG008T_{p}_d46', k) for p in ('Illumina', 'PacBio', 'ONT') for k in ('SNV', 'INDEL')],
           'COLO829T': [(f'COLO829T_{p}_d46', k) for p in ('Illumina', 'fiberseq', 'ONT') for k in ('SNV', 'INDEL')] +
                       [(f'COLO829T_{p}_d46_snvaf05', 'SNV') for p in ('Illumina', 'fiberseq')]}


def s1_files(tumor):
    files = [f'{DATA}/s1/{s}.{k}.chr{c}.tsv' for s, k in S1_SETS[tumor] for c in range(1, 23)]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        raise SystemExit(f'{tumor}: {len(missing)} of {len(files)} s1 tables missing, e.g. {missing[0]}')
    return files


def site_window(q, a, b, max_period=6):
    """[a - 1, max(a, b) + 1) widened to every maximal periodic stretch (period 1..max_period, at least two copies and
    3 bases) that overlaps it."""
    n = len(q)
    lo, hi = max(0, a - 1), min(n, max(a, b) + 1)
    wl, wh = lo, hi
    for p in range(1, max_period + 1):
        k = max(0, lo - p)
        while k < min(n - p, hi):
            if q[k] != q[k + p]:
                k += 1
                continue
            s = k
            while s > 0 and q[s - 1] == q[s - 1 + p]:
                s -= 1
            e = k
            while e + 1 + p < n and q[e + 1] == q[e + 1 + p]:
                e += 1
            t = e + 1 + p                       # stretch = bases [s, t)
            if t - s >= max(3, 2 * p) and s < hi and t > lo:
                wl, wh = min(wl, s), max(wh, t)
            k = e + 1
    return wl, wh


def placements(q, other):
    """equivalent placements of a pure length difference between q and other: when q is longer, the intervals [x, x + L)
    of q whose removal gives other; when q is shorter, the junctions (x, x) of q where inserting gives other."""
    L = len(q) - len(other)
    if L == 0:
        return []
    long_, short = (q, other) if L > 0 else (other, q)
    L = abs(L)
    p = 0
    while p < len(short) and long_[p] == short[p]:
        p += 1
    s = 0
    while s < len(short) and long_[-1 - s] == short[-1 - s]:
        s += 1
    xs = range(max(0, len(short) - s), min(p, len(short)) + 1)
    if len(q) > len(other):
        return [(x, x + L) for x in xs]
    return [(x, x) for x in xs]


def event_window(q, a, b, other):
    """site window of q widened to every equivalent placement of the event (other = the swapped context)."""
    wl, wh = site_window(q, a, b)
    for lo, hi in placements(q, other):
        wl, wh = min(wl, max(0, lo - 1)), max(wh, min(len(q), hi + 1))
    return wl, wh


def windows(r):
    """(alt window, ref window) of an s1 row, each in its own query's coordinates."""
    a, b = site_of(r)
    ra, rb = map(int, ref_site(r).split('-'))
    return event_window(r['alt_q'], a, b, r['ref_q']), event_window(r['ref_q'], ra, rb, r['alt_q'])


def usable(r):
    """'ok' when the row gives an anchored ALT query; else the reason (no_query / short_flank)."""
    if r['status'] != 'ok' or len(r['alt_q']) < MIN_LEN:
        return 'no_query'
    a, b = site_of(r)
    if min(a, len(r['alt_q']) - b) < MIN_FLANK:
        return 'short_flank'
    (wl, wh), (rl, rh) = windows(r)
    if min(wl, len(r['alt_q']) - wh, rl, len(r['ref_q']) - rh) < MIN_RUN_FLANK:
        return 'short_flank'
    return 'ok'


def runs(q, a, b):
    """(longest homopolymer, most copies of a period 2-6 run) touching [a - 1, max(a, b) + 1)."""
    lo, hi = max(0, a - 1), min(len(q), max(a, b) + 1)
    hp, cp = 0, 0
    for p in range(1, 7):
        k = max(0, lo - p)
        while k < min(len(q) - p, hi):
            if q[k] != q[k + p]:
                k += 1
                continue
            s, e = k, k
            while s > 0 and q[s - 1] == q[s - 1 + p]:
                s -= 1
            while e + 1 + p < len(q) and q[e + 1] == q[e + 1 + p]:
                e += 1
            n = e + 1 + p - s
            if s < hi and e + 1 + p > lo:
                if p == 1:
                    hp = max(hp, n)
                else:
                    cp = max(cp, n // p)
            k = e + 1
    return hp, cp


def context(r):
    a, b = site_of(r)
    ra, rb = map(int, ref_site(r).split('-'))
    h1, c1 = runs(r['alt_q'], a, b)
    h2, c2 = runs(r['ref_q'], ra, rb)
    if max(h1, h2) >= 6:
        return 'hp'
    if max(c1, c2) >= 3:
        return 'str'
    return 'other'


def platform(name):
    """platform of a set name (HG008T_PacBio_d46) or a count label (T_PacBio, N_HiFi)."""
    p = name.split('_')[1]
    return {'HiFi': 'PacBio'}.get(p, p)

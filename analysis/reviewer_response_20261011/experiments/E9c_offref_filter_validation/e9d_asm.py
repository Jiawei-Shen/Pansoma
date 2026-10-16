"""E9d: assembly check of call contexts (HG008-T v3.2 tumor and HG008-N v6.2 normal assemblies; COLO829BL verkko
hapX / hapY normal assembly), with minimap2 -x sr like offref_negatives_20260929.

usage:
  e9d_asm.py queries SAMPLE         write DATA/asm/SAMPLE/queries.fa + queries.tsv from every DATA/runs/*/calls.tsv
                                    (and DATA/gold/*/calls.tsv) of that sample: alt_read, alt_graph, ref_graph contexts
  e9d_asm.py classify SAMPLE        parse DATA/asm/SAMPLE/<target>.paf -> DATA/asm/SAMPLE/hits.tsv
                                    (per query and target: passing hits)

A passing hit: alignment covering >= 90 % of the query, no inserted or deleted base anywhere in the alignment, at most
one mismatch (NM <= 1; an assembly or consensus error), and no mismatch within one base of the site interval (for a
deletion: the two bases around the junction). Indels are never allowed, so a homopolymer or repeat-length difference is
never a hit (minimap2 may place the gap anywhere in the run).
"""
import csv
import glob
import hashlib
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e9c_lib as L  # noqa: E402

CTX = ('alt_read', 'alt_graph', 'ref_graph', 'ref_read')
CS = re.compile(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})')


def qid(seq, site):
    return 'q' + hashlib.sha1(f'{seq}|{site}'.encode()).hexdigest()[:16]


def call_tables(sample):
    for p in sorted(glob.glob(f'{L.DATA}/runs/*/calls.tsv') + glob.glob(f'{L.DATA}/gold/*/calls.tsv')):
        rj = os.path.join(os.path.dirname(p), 'run.json')
        import json
        if json.load(open(rj)).get('sample') == sample:
            yield p


def queries(sample):
    out = f'{L.DATA}/asm/{sample}'
    os.makedirs(out, exist_ok=True)
    seen = {}
    csv.field_size_limit(sys.maxsize)
    for p in call_tables(sample):
        for r in csv.DictReader(open(p), delimiter='\t'):
            for c in CTX:
                s, site = r.get(c, ''), r.get(f'{c}_site', '')
                if s and site and len(s) >= 25:
                    seen.setdefault(qid(s, site), (s, site))
    with open(f'{out}/queries.fa.tmp', 'w') as fa, open(f'{out}/queries.tsv.tmp', 'w') as tsv:
        tsv.write('qid\tqlen\tsite\n')
        for k, (s, site) in sorted(seen.items()):
            fa.write(f'>{k}\n{s}\n')
            tsv.write(f'{k}\t{len(s)}\t{site}\n')
    os.replace(f'{out}/queries.fa.tmp', f'{out}/queries.fa')
    os.replace(f'{out}/queries.tsv.tmp', f'{out}/queries.tsv')
    print(f'{sample}: {len(seen):,} unique query contexts -> {out}/queries.fa')


def site_exact(cs, start, lo, hi):
    """cs walk from query position `start` (in the alignment's query orientation): True when [lo, hi) is inside the
    alignment and has no mismatch / insertion inside and no deletion strictly inside."""
    p = start
    if lo < start:
        return False
    for op in CS.findall(cs):
        t = op[0]
        if t == ':':
            p += int(op[1:])
        elif t == '*':
            if lo <= p < hi:
                return False
            p += 1
        elif t == '+':
            n = len(op) - 1
            if p < hi and p + n > lo:
                return False
            p += n
        elif t == '-':
            if lo < p < hi:
                return False
        elif t == '~':
            return False
        if p >= hi and t == ':':
            return True
    return p >= hi


def classify(sample):
    out = f'{L.DATA}/asm/{sample}'
    q = {r['qid']: r for r in csv.DictReader(open(f'{out}/queries.tsv'), delimiter='\t')}
    hits = defaultdict(set)   # (qid, target) -> {(contig, pos // 50)}
    for target, role, _ in L.ASSEMBLIES[sample]:
        paf = f'{out}/{target}.paf'
        if not os.path.exists(paf):
            raise SystemExit(f'missing {paf}')
        for line in open(paf):
            f = line.rstrip('\n').split('\t')
            name, qlen, qs, qe, strand = f[0], int(f[1]), int(f[2]), int(f[3]), f[4]
            tags = dict((t[:2], t[5:]) for t in f[12:])
            if qe - qs < 0.9 * qlen or int(tags.get('NM', 99)) > 1 or 'cs' not in tags:
                continue
            if '+' in tags['cs'] or '-' in tags['cs']:
                continue
            a, b = map(int, q[name]['site'].split('-'))
            lo, hi = max(0, a - 1), min(qlen, max(b, a) + 1)
            if strand == '-':
                lo, hi = qlen - hi, qlen - lo
                start = qlen - qe
            else:
                start = qs
            if site_exact(tags['cs'], start, lo, hi):
                hits[(name, target)].add((f[5], int(f[7]) // 50))
    with open(f'{out}/hits.tsv', 'w') as w:
        w.write('qid\ttarget\trole\tn_loci\n')
        for (name, target), loci in sorted(hits.items()):
            role = next(r for t, r, _ in L.ASSEMBLIES[sample] if t == target)
            w.write(f'{name}\t{target}\t{role}\t{len(loci)}\n')
    print(f'{sample}: {len(hits):,} (query, target) pairs with a passing hit; {len(q):,} queries')


if __name__ == '__main__':
    {'queries': queries, 'classify': classify}[sys.argv[1]](sys.argv[2])

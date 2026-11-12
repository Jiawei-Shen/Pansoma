"""s2: unique query sequences of one tumor (all three platforms) for the assembly and GRCh38 alignments.

usage: s2_queries.py TUMOR       (HG008T | COLO829T)

From every s1 row that common.usable() accepts (status ok, alt_q >= 30 bp, >= 10 read-consensus bases on each side of the
site and >= 5 beyond the event window on both queries): the ALT context alt_q (site [a, b)) and the REF-swap context
ref_q (site from ref_site(); written only when >= 30 bp). Query id = 'q' + sha1(seq|site)[:16], so identical contexts of
different tensors or platforms map once. win = the union, over the rows that use the query, of its event window
(common.event_window), used by s3 to rank normal alignments by their differences outside the window.
Writes DATA/s2/<TUMOR>/queries.fa, queries.tsv (qid, qlen, site, win).
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, MIN_LEN, qid, ref_site, s1_files, usable, windows  # noqa: E402,F401


def main():
    tumor = sys.argv[1]
    out = f'{DATA}/s2/{tumor}'
    os.makedirs(out, exist_ok=True)
    seen = {}
    n = {}
    for p in s1_files(tumor):
        for r in csv.DictReader(open(p), delimiter='\t'):
            u = usable(r)
            n[u] = n.get(u, 0) + 1
            if u != 'ok':
                continue
            (wl, wh), (rl, rh) = windows(r)
            ctx = [(r['alt_q'], r['site'], wl, wh)]
            if len(r['ref_q']) >= MIN_LEN:
                ctx.append((r['ref_q'], ref_site(r), rl, rh))
            for s, site, lo, hi in ctx:
                k = qid(s, site)
                if k in seen:
                    v = seen[k]
                    v[2], v[3] = min(v[2], lo), max(v[3], hi)
                else:
                    seen[k] = [s, site, lo, hi]
    with open(f'{out}/queries.fa.tmp', 'w') as fa, open(f'{out}/queries.tsv.tmp', 'w') as tsv:
        tsv.write('qid\tqlen\tsite\twin\n')
        for k, (s, site, lo, hi) in sorted(seen.items()):
            fa.write(f'>{k}\n{s}\n')
            tsv.write(f'{k}\t{len(s)}\t{site}\t{lo}-{hi}\n')
    os.replace(f'{out}/queries.tsv.tmp', f'{out}/queries.tsv')
    os.replace(f'{out}/queries.fa.tmp', f'{out}/queries.fa')
    print(f'{tumor}: s1 rows {n}, unique queries {len(seen):,} -> {out}/queries.fa')


if __name__ == '__main__':
    main()

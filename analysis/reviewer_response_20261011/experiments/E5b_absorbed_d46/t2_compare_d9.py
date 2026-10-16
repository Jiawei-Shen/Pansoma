"""E5b check t2: s4's d9 answer (vg-find subgraph) against t1's (complete d9 edge table) on the same HG008T loci.
Run after s4 d9 and t1: python t2_compare_d9.py"""
import csv, gzip
from collections import Counter
from e5b_common import DATA
t1 = {r['truth_id']: r for r in csv.DictReader(open(f'{DATA}/t1_spell_d9_sqlite.tsv'), delimiter='\t')}
s4 = {r['truth_id']: r for r in csv.DictReader(gzip.open(f'{DATA}/spell_d9.tsv.gz', 'rt'), delimiter='\t') if r['sample'] == 'HG008T'}
c = Counter((t1[k]['ref_ok'], s4[k]['ref_ok'], t1[k]['alt_spelled'], s4[k]['alt_spelled']) for k in t1)
print('(ref_ok sqlite, ref_ok vg-find, alt sqlite, alt vg-find): count over', len(t1), 'loci')
for k, v in sorted(c.items()):
    print(k, v)
print('agree on alt_spelled:', sum(t1[k]['alt_spelled'] == s4[k]['alt_spelled'] for k in t1), 'of', len(t1))

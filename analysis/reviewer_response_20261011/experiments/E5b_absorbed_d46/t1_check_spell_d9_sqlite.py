"""E5b check t1: the spelling search of s4 on the d9 graph with the complete d9 edge table (the path_index sqlite of
graph_absorbed_somatic_20261001 / classify.py) instead of a vg-find subgraph, against s2_graph_paths' graph_match of the
2,478 HG008T read-level absorbed loci (exact = a d9 path spells GRCh38 + truth; with_germline / closest = it does not).
Expect: exact -> spelled, closest -> not spelled; with_germline -> mostly not spelled (it needs germline alleles too).
Also writes the per-locus answer, so s4's vg-find result can be compared with it (t2).
Run: python t1_check_spell_d9_sqlite.py [n_loci]   (login node; HG008T loci of the read-level set)
"""
import csv, gzip, sqlite3, sys
from collections import Counter
import pysam
from e5b_common import DATA, FASTA, RefPath, spells, FLIP

IDX = '/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.GRCh38_CHM13.path_index.sqlite'
PV = '/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/per_variant.tsv'
csv.field_size_limit(1 << 30)


class SqliteGraph:
    """SubGraph interface over the complete d9 node / edge tables."""
    RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

    def __init__(self, path):
        self.db = sqlite3.connect(f'file:{path}?mode=ro', uri=True)
        self._seq, self._succ = {}, {}
        self.seq = self

    def __contains__(self, n):
        return self._get(n) is not None

    def _get(self, n):
        if n not in self._seq:
            r = self.db.execute('select seq from nodes where node_id=?', (str(n),)).fetchone()
            self._seq[n] = r[0] if r else None
        return self._seq[n]

    def oseq(self, n, o):
        s = self._get(n)
        return s if o == '+' else s.translate(self.RC)[::-1]

    def succ(self, n, o):
        if (n, o) not in self._succ:
            out = {(int(m), mo) for m, mo in self.db.execute('select to_id, to_orient from edges where from_id=? and from_orient=?', (str(n), o))}
            out |= {(int(m), FLIP[mo]) for m, mo in self.db.execute('select from_id, from_orient from edges where to_id=? and to_orient=?', (str(n), FLIP[o]))}
            self._succ[(n, o)] = sorted(out)
        return self._succ[(n, o)]


n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 9
R = RefPath('d9')
G = SqliteGraph(IDX)
fa = pysam.FastaFile(FASTA)
pv = {(r['chrom'], r['vcf_pos'], r['vcf_ref'], r['vcf_alt']): r['graph_match'] for r in csv.DictReader(open(PV), delimiter='\t')}
loci = [r for r in csv.DictReader(gzip.open(f'{DATA}/loci.tsv.gz', 'rt'), delimiter='\t')
        if r['sample'] == 'HG008T' and (r['chrom'], r['vcf_pos'], r['vcf_ref'], r['vcf_alt']) in pv][:n_max]
tab = Counter()
with open(f'{DATA}/t1_spell_d9_sqlite.tsv', 'w') as w:
    w.write('truth_id\tgraph_match\tref_ok\talt_spelled\tinfo\n')
    for r in loci:
        chrom, pos0, vref, valt = r['chrom'], int(r['vcf_pos']) - 1, r['vcf_ref'].upper(), r['vcf_alt'].upper()
        first, last, lo, hi = int(r['d9_first']), int(r['d9_last']), int(r['d9_win_start0']), int(r['d9_win_end0'])
        region = fa.fetch(chrom, lo, hi + 1).upper()
        off = pos0 - lo
        alt = region[:off] + valt + region[off + len(vref):]
        okr = spells(G, R, chrom, first, last, lo, hi, region)[0]
        oka, _, info = spells(G, R, chrom, first, last, lo, hi, alt)
        gm = pv[(r['chrom'], r['vcf_pos'], r['vcf_ref'], r['vcf_alt'])]
        tab[(gm, int(okr), int(oka))] += 1
        w.write(f"{r['truth_id']}\t{gm}\t{int(okr)}\t{int(oka)}\t{info}\n")
print('graph_match (s2), GRCh38 spelled, ALT spelled: count')
for k, v in sorted(tab.items()):
    print(k, v)

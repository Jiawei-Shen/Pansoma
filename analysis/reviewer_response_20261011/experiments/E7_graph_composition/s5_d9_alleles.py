"""E7 step 5: the alleles the real HPRC v1.1 d9 graph contains, from its own vg deconstruct VCF.

Input: /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz (deconstruct of the d9 graph Pansoma
aligns to: top-level snarls, same 45 columns). Not vcfwave-decomposed, so a multi-variant ALT allele stays one COMPLEX key
(after trimming); SNVs, simple indels and MNPs (split into SNVs) are comparable with the s1 / s2 keys. This gives a lower
bound on what d9 spells; used to (i) measure how sharp the d9 frequency filter is per HPRC carrier count and (ii) check
directly which sample alleles d9 contains.
`bcftools norm -m -any -f GRCh38 -c w` then common.primitives, as in s1.
Output: WORK/d9_alleles/CHROM.tsv.gz: key, kind, AC (d9 INFO AC of the record allele), complex parent (0/1).
Usage: s5_d9_alleles.py CHROM
"""
import gzip, json, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import BCFTOOLS, FASTA, WORK, key, primitives

D9 = '/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz'
chrom = sys.argv[1]
os.makedirs(f'{WORK}/d9_alleles', exist_ok=True)
cmd = (f'{BCFTOOLS} view -r {chrom} {D9} -Ou | {BCFTOOLS} norm -m -any -f {FASTA} -c w -Ou '
       f'2> {WORK}/d9_alleles/{chrom}.norm.log | {BCFTOOLS} query -f "%POS\\t%REF\\t%ALT\\t%INFO/AC\\n"')
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
seen = {}
stat = dict(records=0, primitives=0, complex=0)
for line in p.stdout:
    pos, ref, alt, ac = line.rstrip('\n').split('\t')
    stat['records'] += 1
    if alt in ('*', '.') or alt.startswith('<'):
        continue
    prims = primitives(int(pos), ref, alt)
    for pp, r, a, kind in prims:
        stat['primitives'] += 1
        stat['complex'] += kind == 'COMPLEX'
        k = key(chrom, pp, r, a)
        acv = int(ac) if ac not in ('.', '') else 0
        if k in seen:
            seen[k][1] += acv
        else:
            seen[k] = [kind, acv]
assert p.wait() == 0
with gzip.open(f'{WORK}/d9_alleles/{chrom}.tsv.gz.tmp', 'wt', compresslevel=3) as fo:
    fo.write('key\tkind\tAC\n')
    for k, (kind, ac) in seen.items():
        fo.write(f'{k}\t{kind}\t{ac}\n')
os.replace(f'{WORK}/d9_alleles/{chrom}.tsv.gz.tmp', f'{WORK}/d9_alleles/{chrom}.tsv.gz')
stat['keys'] = len(seen)
json.dump(stat, open(f'{WORK}/d9_alleles/{chrom}.json', 'w'), indent=1)
print(chrom, stat)

"""E7 step 1: per-haplotype carriers of every HPRC v1.1 primitive allele on one chromosome.

Inputs
  wave: hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz (vg deconstruct -> vcfbub -l 0 -a 100000 -> vcfwave; 45 columns:
        CHM13 haploid + 44 diploid samples = 89 non-GRCh38 haplotypes; GRCh38 is the reference and never a carrier).
  raw:  hprc-v1.1-mc-grch38.raw.vcf.gz (the vg deconstruct VCF before vcfbub / vcfwave; tabix index from
        analysis/graph_absorbed_somatic_20261001/s5_index_full_vcf.sh).
Why both: in the wave VCF the genotypes of the vcfwave-decomposed records (INFO ORIGIN) are incomplete. On chr1, 82,346
of 533,408 decomposed records carry fewer GT carriers than their INFO AC (36,985 have none; GT total 4.53 M vs AC 5.42 M),
while the 1.61 M records that vcfwave left alone agree with AC. vcfwave merges a primitive shared by several parent
alleles into one record and keeps the genotypes of only some of them. So, for every decomposed bubble (wave ID
'<raw ID>_<n>'), each ALT allele of the raw record is re-expressed as a set of that bubble's wave primitives
(min-size exact cover, `decompose`), and a primitive's carriers are the union of the raw carriers of every parent
allele that contains it. Records vcfwave did not touch keep their wave GT.
Then every allele goes through `bcftools norm -f GRCh38 -c w` (left-align) and common.primitives (MNP -> SNVs), and
identical keys are merged (union of carriers).
Haplotype bit order: bit 0 = CHM13, bit 2i-1 / 2i = hap 1 / hap 2 of the i-th sample column (i = 1..44).
Output: WORK/hprc_alleles/<tag>.tsv.gz: key, kind, rep_class, unit, tract, ilen, n_carriers, mask (hex);
        WORK/hprc_alleles/<tag>.json: counts, including the check of the repaired carriers against wave INFO AC.
Usage: s1_hprc_alleles.py CHROM [REGION] [TAG]
"""
import collections, gzip, json, os, subprocess, sys
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import BCFTOOLS, FASTA, HPRC_VCF, WORK, key, primitives, repeat_context

RAW = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/full_vcf_index/hprc-v1.1-mc-grch38.raw.vcf.gz'
chrom = sys.argv[1]
region = sys.argv[2] if len(sys.argv) > 2 else chrom
tag = sys.argv[3] if len(sys.argv) > 3 else chrom
sys.setrecursionlimit(100000)
stat = collections.Counter()


def query(vcf, fmt, region):
    p = subprocess.Popen([BCFTOOLS, 'query', '-r', region, '-f', fmt, vcf], stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 20)
    for line in p.stdout:
        yield line.rstrip('\n').split('\t')
    assert p.wait() == 0, (vcf, region)


def gt_masks(gts, n_alleles):
    """allele index -> haplotype mask, from the 45 GT strings."""
    m = [0] * (n_alleles + 1)
    for i, g in enumerate(gts):
        a = g.replace('/', '|').split('|')
        bits = [0] if i == 0 else [2 * i - 1, 2 * i]
        for b, x in zip(bits, a):
            if x not in ('.', '0'):
                m[int(x)] |= 1 << b
    return m


def decompose(ref, alt, prims):
    """Smallest set of primitives (offset, ref, alt) that turns ref into alt; None if there is none."""
    prims = sorted(prims)
    offs = [p[0] for p in prims]
    memo = {}

    def solve(rp, ap):
        k = (rp, ap)
        if k in memo:
            return memo[k]
        best = () if ref[rp:] == alt[ap:] else None
        j = 0
        while j < len(prims) and offs[j] < rp:
            j += 1
        for o, r, a in prims[j:]:
            gap = o - rp
            if ref[rp:o] != alt[ap:ap + gap]:
                break  # every later primitive needs this gap too
            if ref[o:o + len(r)] != r or alt[ap + gap:ap + gap + len(a)] != a:
                continue
            sub = solve(o + len(r), ap + gap + len(a))
            if sub is not None and (best is None or len(sub) + 1 < len(best)):
                best = ((o, r, a),) + sub
        memo[k] = best
        return best

    return solve(0, 0)


# 1. wave records
plain, groups = [], collections.defaultdict(list)
for f in query(HPRC_VCF, '%ID\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/ORIGIN[\t%GT]\n', region):
    vid, pos, ref, alts, acs, origin, gts = f[0], int(f[1]), f[2], f[3].split(','), f[4].split(','), f[5], f[6:]
    assert len(gts) == 45
    m = gt_masks(gts, len(alts)) if '1' in ''.join(gts) or len(alts) > 1 else [0] * (len(alts) + 1)
    stat['wave_records'] += 1
    if origin == '.':
        for j, a in enumerate(alts, start=1):
            plain.append((pos, ref, a, m[j]))
    else:
        stat['wave_records_decomposed'] += 1
        assert len(alts) == 1
        groups[vid.rsplit('_', 1)[0]].append((pos, ref, alts[0], int(acs[0]), m[1]))

# 2. raw records of the decomposed bubbles
raw = {}
for f in query(RAW, '%ID\t%POS\t%REF\t%ALT[\t%GT]\n', region):
    if f[0] in groups and f[0] not in raw:
        alts = f[3].split(',')
        raw[f[0]] = (int(f[1]), f[2], alts, gt_masks(f[4:], len(alts)))
stat['decomposed_bubbles'] = len(groups)
stat['decomposed_bubbles_without_raw'] = sum(1 for g in groups if g not in raw)

# 3. repaired carriers of the decomposed primitives
records = list(plain)
check = collections.Counter()
for gid, prims in groups.items():
    uniq = {}
    for pos, ref, alt, ac, m in prims:
        k = (pos, ref, alt)
        if k in uniq:
            uniq[k][0] = max(uniq[k][0], ac)
            uniq[k][1] |= m
        else:
            uniq[k] = [ac, m]
    rep = {k: 0 for k in uniq}
    if gid in raw:
        rpos, rref, ralts, rm = raw[gid]
        cand = [(p - rpos, r, a) for (p, r, a) in uniq if 0 <= p - rpos and p - rpos + len(r) <= len(rref)]
        stat['primitives_outside_raw_ref'] += len(uniq) - len(cand)
        for ai, a in enumerate(ralts, start=1):
            if not rm[ai] or a.startswith('<') or a == '*':
                continue
            d = decompose(rref, a, cand) if len(rref) + len(a) <= 200000 else None
            if d is None:
                stat['raw_alleles_not_decomposed'] += 1
                stat['raw_allele_carriers_not_decomposed'] += bin(rm[ai]).count('1')
                continue
            stat['raw_alleles_decomposed'] += 1
            for o, r, aa in d:
                rep[(o + rpos, r, aa)] |= rm[ai]
    for k, (ac, m) in uniq.items():
        final = rep[k] | m
        n = bin(final).count('1')
        check['n==AC' if n == ac else ('n<AC' if n < ac else 'n>AC')] += 1
        check['sum_AC'] += ac
        check['sum_n'] += n
        check['sum_wave_gt'] += bin(m).count('1')
        records.append((k[0], k[1], k[2], final))
stat.update({f'decomposed_primitive_{k}': v for k, v in check.items()})

# 4. left-align with bcftools norm (the mask rides in the ID column), split MNPs, merge identical keys
fai = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(FASTA + '.fai')}
records.sort(key=lambda r: r[0])
norm = subprocess.Popen(f'{BCFTOOLS} norm -f {FASTA} -c w -Ov - 2> {WORK}/hprc_alleles/{tag}.norm.log | '
                        f'{BCFTOOLS} query -f "%POS\\t%ID\\t%REF\\t%ALT\\n"', shell=True, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
os.makedirs(f'{WORK}/hprc_alleles', exist_ok=True)
alleles = {}


def reader():
    for line in norm.stdout:
        pos, mid, ref, alt = line.rstrip('\n').split('\t')
        stat['norm_records'] += 1
        for p, r, a, kind in primitives(int(pos), ref, alt):
            k = key(chrom, p, r, a)
            v = alleles.get(k)
            if v is None:
                alleles[k] = [int(mid, 16), kind, p, r, a]
            else:
                v[0] |= int(mid, 16)


import threading
t = threading.Thread(target=reader)
t.start()
norm.stdin.write(f'##fileformat=VCFv4.2\n##contig=<ID={chrom},length={fai[chrom]}>\n'
                 '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
for pos, ref, alt, m in records:
    if alt in ('*', '.') or alt.startswith('<') or ref == alt:
        continue
    norm.stdin.write(f'{chrom}\t{pos}\t{m:x}\t{ref}\t{alt}\t.\t.\t.\n')
norm.stdin.close()
t.join()
assert norm.wait() == 0

fa = pysam.FastaFile(FASTA)
stat['keys'] = len(alleles)
stat['keys_with_carrier'] = sum(1 for v in alleles.values() if v[0])
with gzip.open(f'{WORK}/hprc_alleles/{tag}.tsv.gz.tmp', 'wt', compresslevel=3) as fo:
    fo.write('key\tkind\trep_class\tunit\ttract\tilen\tn_carriers\tmask\n')
    for k, (mask, kind, p, r, a) in alleles.items():
        rc, u, tr = repeat_context(fa, chrom, p, r, a, kind)
        fo.write(f'{k}\t{kind}\t{rc}\t{u}\t{tr}\t{len(a) - len(r)}\t{bin(mask).count("1")}\t{mask:x}\n')
os.replace(f'{WORK}/hprc_alleles/{tag}.tsv.gz.tmp', f'{WORK}/hprc_alleles/{tag}.tsv.gz')
json.dump(dict(stat), open(f'{WORK}/hprc_alleles/{tag}.json', 'w'), indent=1)
print(chrom, region, dict(stat))

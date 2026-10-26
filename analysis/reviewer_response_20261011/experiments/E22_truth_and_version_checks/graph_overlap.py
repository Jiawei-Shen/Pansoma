"""E22 check 3: how many somatic truth alleles of each tumor are already alleles of the HPRC v2.1 pangenome (R1.m1).

Same method for all six tumors. Truth alleles come from E7 step 2 (`WORK/sample_alleles/<sample>.somatic.tsv.gz`:
left-aligned with bcftools norm, split into primitive keys; HCC1954 / H1437 / H2009 added on 2026-10-26 with
`s2_sample_alleles.py <LINE> <LINE>T_PacBio_d46 somatic`). HPRC carriers come from E7 step 1 on the v2.1 wave VCF
(`WORK/v2.1/hprc_alleles/chr*.tsv.gz`: per primitive key, the carrier count and mask over the 463 non-GRCh38
haplotypes). A truth allele with several primitives counts only when one haplotype set carries all of them (AND of the
masks), as in E7.
- in_hprc: at least 1 of the 463 haplotypes carries the allele.
- in_d46: at least 46 carry it, i.e. the allele is on a node the d46 graph keeps (E7 s4 'predicted_d46').
Scopes: bed = FILTER PASS and inside the truth's benchmark BED (as E7 'bed'), chr1-22; bed_chr1 = the same on chr1.
For the CASTLE lines and HCC1395 the alleles are also joined (chrom, VCF pos / ref / alt of the truth table) to E22's
per-allele table (`<set>.alleles.tsv.gz`: normal genotype, gnomAD / CoLoRSdb AF, tumor VAF, het class of the normal's
1-Mb bin) to split by loss of heterozygosity in the normal (het fraction < 0.2) and to list examples.
Outputs (this folder): graph_overlap.tsv, graph_overlap_examples.tsv; WORK copy of the per-allele join.
"""
import collections, csv, gzip, os, sys

E7 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'E7_graph_composition')
os.environ['E7_GRAPH'] = 'v2.1'
sys.path.insert(0, E7)
from common import AUTOSOMES, GWORK, WORK  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
E22W = '/scratch/jshen/data/reviewer_experiments_20261013/E22_truth_and_version_checks'
SETS = '/scratch/jshen/data/pansoma_v2_tensors'
SAMPLES = [  # sample, tensor set holding the truth table, E22 per-allele set (None: not in E22)
    ('HG008', 'HG008T_PacBio', None),
    ('COLO829', 'COLO829T_fiberseq', None),
    ('HCC1395', 'HCC1395T_PacBio_d46', 'HCC1395_SEQC2'),
    ('HCC1954', 'HCC1954T_PacBio_d46', 'HCC1954_CASTLE'),
    ('H1437', 'H1437T_PacBio_d46', 'H1437_CASTLE'),
    ('H2009', 'H2009T_PacBio_d46', 'H2009_CASTLE'),
]
LOH_CLASSES = {'het0.0-0.1', 'het0.1-0.2'}

# 1. truth alleles: keys (E7 s2) + raw VCF allele (truth table) + E22 annotations
rows = []
need = collections.defaultdict(set)
for s, tset, e22 in SAMPLES:
    raw = {}
    with open(f'{SETS}/{tset}/truth/somatic.graph.tsv') as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {c: hdr.index(c) for c in ('truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt')}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            raw[f[ix['truth_id']]] = (f[ix['chrom']], f[ix['vcf_pos']], f[ix['vcf_ref']], f[ix['vcf_alt']])
    ann = {}
    if e22:
        with gzip.open(f'{E22W}/{e22}.alleles.tsv.gz', 'rt') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                ann[(r['chrom'], r['pos'], r['ref'], r['alt'])] = r
    with gzip.open(f'{WORK}/sample_alleles/{s}.somatic.tsv.gz', 'rt') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['chrom'] not in AUTOSOMES:
                continue
            keys = r['keys'].split(';') if r['keys'] else []
            for k in keys:
                need[r['chrom']].add(k)
            a = ann.get(raw[r['truth_id']]) if e22 else None
            rows.append(dict(sample=s, truth_id=r['truth_id'], chrom=r['chrom'], raw=raw[r['truth_id']],
                             kind=r['nkind'], bed=r['passed'] == 'T' and r['in_bed'] == 'T', keys=keys,
                             rep_class=r['rep_class'], tract=r['tract'], ann=a, e22=e22 is not None))
print('truth alleles', len(rows), 'keys', sum(len(v) for v in need.values()), flush=True)

# 2. HPRC v2.1 carrier masks of the needed keys
masks = {}
for c in AUTOSOMES:
    nk = need.get(c, set())
    with gzip.open(f'{GWORK}/hprc_alleles/{c}.tsv.gz', 'rt') as fh:
        fh.readline()
        for line in fh:
            k, rest = line.split('\t', 1)
            if k in nk:
                masks[k] = int(line.rstrip('\n').rsplit('\t', 1)[1], 16)
    print(c, len(masks), flush=True)

# 3. carriers per truth allele
for r in rows:
    ms = [masks.get(k) for k in r['keys']]
    if r['keys'] and all(m is not None for m in ms):
        m = ms[0]
        for x in ms[1:]:
            m &= x
        r['carriers'] = bin(m).count('1')
    else:
        r['carriers'] = 0
    r['stratum'] = 'SNV' if r['kind'] == 'SNV' else ('INDEL' if r['kind'] in ('INS', 'DEL') else 'OTHER')
    a = r['ann']
    r['loh'] = None if a is None else (a['normal_bin_class'] in LOH_CLASSES)

# 4. summary
out = []
for s, _, e22 in SAMPLES:
    for scope in ('bed', 'bed_chr1', 'all'):
        sel0 = [r for r in rows if r['sample'] == s and (scope == 'all' or r['bed'])
                and (scope != 'bed_chr1' or r['chrom'] == 'chr1')]
        for stratum in ('SNV', 'INDEL'):
            for part in (('all', 'normal_not_LOH', 'normal_LOH') if e22 else ('all',)):
                sel = [r for r in sel0 if r['stratum'] == stratum and
                       (part == 'all' or (part == 'normal_LOH') == bool(r['loh']))]
                n = len(sel)
                n1 = sum(r['carriers'] >= 1 for r in sel)
                n46 = sum(r['carriers'] >= 46 for r in sel)
                ncommon = sum(r['ann'] is not None and r['ann']['gnomad_af'] not in ('', '.') and
                              float(r['ann']['gnomad_af']) >= 0.01 for r in sel) if e22 else ''
                out.append(dict(sample=s, scope=scope, stratum=stratum, part=part, n=n, in_hprc=n1,
                                in_d46=n46, pct_in_hprc=round(100 * n1 / n, 2) if n else '',
                                pct_in_d46=round(100 * n46 / n, 2) if n else '',
                                gnomad_ge_1pct=ncommon,
                                pct_gnomad_ge_1pct=round(100 * ncommon / n, 2) if (n and e22) else ''))
with open(f'{HERE}/graph_overlap.tsv', 'w') as fo:
    w = csv.DictWriter(fo, list(out[0]), delimiter='\t')
    w.writeheader()
    w.writerows(out)

# 5. per-allele join (for examples) and an examples table: in_d46 SNVs in the benchmark BED of the CASTLE lines,
# a fixed spread over chromosomes (every 1/8 of the sorted list), with E22's annotations
with gzip.open(f'{E22W}/graph_overlap.alleles.tsv.gz', 'wt') as fo:
    fo.write('sample\ttruth_id\tchrom\tpos\tref\talt\tstratum\tbed\tcarriers\tnormal_gt\tnormal_bin_class\tgnomad_af\t'
             'colorsdb_af\ttumor_vaf\tn_tech\n')
    for r in rows:
        a = r['ann'] or {}
        fo.write('\t'.join(map(str, [r['sample'], r['truth_id'], *r['raw'], r['stratum'], int(r['bed']), r['carriers'],
                                     a.get('normal_gt', ''), a.get('normal_bin_class', ''), a.get('gnomad_af', ''),
                                     a.get('colorsdb_af', ''), a.get('tumor_vaf', ''), a.get('n_tech', '')])) + '\n')
ex = []
for s, _, e22 in SAMPLES:
    if not e22 or s == 'HCC1395':
        continue
    for loh in (False, True):
        cand = sorted((r for r in rows if r['sample'] == s and r['bed'] and r['stratum'] == 'SNV' and
                       r['carriers'] >= 46 and r['loh'] is loh),
                      key=lambda r: (AUTOSOMES.index(r['chrom']), int(r['raw'][1])))
        for i in range(0, len(cand), max(1, len(cand) // 4))[:4]:
            r = cand[i]
            a = r['ann'] or {}
            ex.append(dict(sample=s, normal_LOH=int(loh), chrom=r['chrom'], pos=r['raw'][1], ref=r['raw'][2],
                           alt=r['raw'][3], hprc_carriers_of_463=r['carriers'], normal_gt=a.get('normal_gt', ''),
                           normal_bin_class=a.get('normal_bin_class', ''), gnomad_af=a.get('gnomad_af', ''),
                           colorsdb_af=a.get('colorsdb_af', ''), tumor_vaf=a.get('tumor_vaf', ''),
                           n_tech=a.get('n_tech', '')))
with open(f'{HERE}/graph_overlap_examples.tsv', 'w') as fo:
    w = csv.DictWriter(fo, list(ex[0]), delimiter='\t')
    w.writeheader()
    w.writerows(ex)
print('done', flush=True)

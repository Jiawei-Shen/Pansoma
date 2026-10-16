"""Shared definitions for E7 (in-silico graph composition).

Allele keys: after `bcftools norm -m -any -f GRCh38` (left-aligned, biallelic, one anchor base for indels), every
allele becomes one or more primitive alleles:
  SNV   ref and alt 1 bp;
  MNP   same length > 1: split into the SNVs at the differing positions (what `bcftools norm --atomize` does for MNPs;
        not used on the whole VCF because it also atomizes 50 bp - 100 kb complex alleles into thousands of SNVs);
  INS / DEL  one side is 1 bp and shares the first base with the other;
  COMPLEX    anything else (kept as one key, matched only exactly).
Long sequences are replaced by their md5 so that structural alleles do not blow up memory; the same function is used
for HPRC and for the sample alleles, so keys stay comparable.
"""
import csv, hashlib, os

# GRCh38 no-alt analysis set, the fasta of every labels manifest. Slurm jobs copy it to node-local /tmp and set E7_FASTA,
# because random access over BeeGFS stalled on busy nodes (2026-10-13).
FASTA_SRC = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
FASTA = os.environ.get('E7_FASTA', FASTA_SRC)
BCFTOOLS = '/wanglab/jshen/anaconda3/envs/sambcfenv/bin/bcftools'
WORK = '/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition'
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
LONG = 200

# Which HPRC release the haplotype axis comes from: E7_GRAPH=v1.1 (default; the 2026-10-13 results in this folder) or
# v2.1 (added 2026-10-16; results in ./v2.1, large outputs in WORK/v2.1). The sample alleles (s2) are shared.
GRAPH = os.environ.get('E7_GRAPH', 'v1.1')
GRAPHS = {
    'v1.1': dict(vcf='/scratch/jshen/data/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz',
                 raw='/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/full_vcf_index/'
                     'hprc-v1.1-mc-grch38.raw.vcf.gz',
                 nhap=89, real='d9', d_real=9, sub=''),
    # HPRC release 2 v2.1 Minigraph-Cactus GRCh38 graph, decomposed VCF (vcfbub -l 0 -a 100000, vcfwave -I 1000,
    # vcffixup; s3://human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38/,
    # downloaded 2026-10-16). Its README says the missing-genotype bug of the v1.1 wave VCF is resolved, so s1 uses its
    # genotypes as they are (raw None; checked against the raw VCF on a test region, see README v2.1 section).
    'v2.1': dict(vcf='/scratch/jshen/data/hprc_v2.1_vcf/hprc-v2.1-mc-grch38.wave.vcf.gz', raw=None,
                 nhap=463, real='d46', d_real=46, sub='v2.1'),
}
G = GRAPHS[GRAPH]
HPRC_VCF = G['vcf']
NHAP = G['nhap']  # non-GRCh38 haplotypes: CHM13 (haploid column 0) + 2 per diploid sample column
NW = (NHAP + 63) // 64  # uint64 words per carrier mask
GWORK = f"{WORK}/{G['sub']}" if G['sub'] else WORK  # graph-specific large outputs (hprc_alleles, joined, hist_draws)
GOUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), G['sub'])  # graph-specific result tables

# 1000 Genomes population -> superpopulation, for the HPRC release 2 metadata (population_abbreviation). ASL (African
# Americans in St. Louis) and MKK (HapMap Maasai) are not 1000G codes; both are African-ancestry populations.
SUPERPOP = {**{p: 'AFR' for p in ('YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB', 'MKK', 'ASL')},
            **{p: 'AMR' for p in ('MXL', 'PUR', 'CLM', 'PEL')},
            **{p: 'EAS' for p in ('CHB', 'JPT', 'CHS', 'CDX', 'KHV')},
            **{p: 'EUR' for p in ('CEU', 'TSI', 'FIN', 'GBR', 'IBS')},
            **{p: 'SAS' for p in ('GIH', 'PJL', 'BEB', 'STU', 'ITU')}}
# GIAB samples without a 1000G population code: HG002 (Ashkenazi Jewish) and HG005 (Han Chinese).
SUPERPOP_EXTRA = {'HG002': 'EUR', 'HG005': 'EAS'}
V2_META = '/scratch/jshen/data/hprc_v2.1_meta/hprc_release2_sample_metadata.csv'


def haplotypes():
    """One dict per carrier-mask bit (sample, hap, superpopulation), in the bit order of s1 (VCF column order)."""
    cols = [l.strip() for l in os.popen(f'{BCFTOOLS} query -l {HPRC_VCF}')]
    assert cols[0] == 'CHM13' and 1 + 2 * (len(cols) - 1) == NHAP, (cols[:2], len(cols))
    if GRAPH == 'v1.1':
        haps = list(csv.DictReader(open('/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/'
                                        'hprc_individuals.tsv'), delimiter='\t'))
    else:
        pop = {r['sample_id']: r['population_abbreviation'] for r in csv.DictReader(open(V2_META))}
        haps = [dict(sample='CHM13', hap='0', superpopulation='reference')]
        for s in cols[1:]:
            sp = SUPERPOP_EXTRA.get(s) or SUPERPOP[pop[s]]
            haps += [dict(sample=s, hap='1', superpopulation=sp), dict(sample=s, hap='2', superpopulation=sp)]
    assert len(haps) == NHAP and haps[0]['sample'] == 'CHM13'
    assert [h['sample'] for h in haps] == [cols[0]] + [s for s in cols[1:] for _ in (1, 2)]
    assert [h['hap'] for h in haps[1:]] == ['1', '2'] * (len(cols) - 1)
    return haps


def mask_words(m):
    """Python int carrier mask -> NW uint64 words (word 0 = bits 0-63)."""
    return [(m >> (64 * i)) & 0xFFFFFFFFFFFFFFFF for i in range(NW)]


def seqkey(s):
    return s if len(s) <= LONG else 'md5:' + hashlib.md5(s.encode()).hexdigest()


def primitives(pos, ref, alt):
    """[(pos, ref, alt, kind)] for one normalized biallelic record (pos 1-based)."""
    ref, alt = ref.upper(), alt.upper()
    if len(ref) == 1 and len(alt) == 1:
        return [(pos, ref, alt, 'SNV')] if ref != alt else []
    if len(ref) == len(alt):
        return [(pos + i, r, a, 'SNV') for i, (r, a) in enumerate(zip(ref, alt)) if r != a]
    if len(ref) == 1 and alt[0] == ref:
        return [(pos, ref, alt, 'INS')]
    if len(alt) == 1 and ref[0] == alt:
        return [(pos, ref, alt, 'DEL')]
    return [(pos, ref, alt, 'COMPLEX')]


def key(chrom, pos, ref, alt):
    return f'{chrom}:{pos}:{seqkey(ref)}:{seqkey(alt)}'


def min_period(s):
    n = len(s)
    for p in range(1, n + 1):
        if n % p == 0 and s[:p] * (n // p) == s:
            return p
    return n


def repeat_context(fa, chrom, pos, ref, alt, kind):
    """(class, unit length, GRCh38 tract length in bp) for a left-aligned INS / DEL; class HP>=7 / HP2-6 / STR / other.

    GRCh38 frame: the tract is the run of the event's minimal unit that starts right after the anchor base (left
    alignment puts the event at the left end of the run). HP>=7: unit 1 and run >= 7 bp; HP2-6: unit 1 and run 2-6 bp;
    STR: unit 2-6 and >= 3 copies in GRCh38; other: the rest. Mirrors the INDEL PoN rule's 'homopolymer >= 7 bp' (GRCh38 frame).
    """
    if kind not in ('INS', 'DEL'):
        return ('', 0, 0)
    seq = alt[1:] if kind == 'INS' else ref[1:]
    if len(seq) > 50:
        return ('SV', 0, 0)
    u = min_period(seq)
    unit = seq[:u]
    win = fa.fetch(chrom, pos, pos + 400).upper()  # 0-based pos = first base after the anchor
    n = 0
    while win.startswith(unit, n * u):
        n += 1
    tract = n * u
    if u == 1 and tract >= 7:
        return ('HP>=7', u, tract)
    if u == 1 and tract >= 2:
        return ('HP2-6', u, tract)
    if 2 <= u <= 6 and n >= 3:
        return ('STR', u, tract)
    return ('other', u, tract)


# strata shared by s3 / s4
SAMPLES = ['HG008', 'COLO829', 'HCC1395']
ROLES = ['germline', 'somatic']
FINE = ['SNV', 'INS_HP>=7', 'INS_HP2-6', 'INS_STR', 'INS_other', 'DEL_HP>=7', 'DEL_HP2-6', 'DEL_STR', 'DEL_other',
        'SV>50', 'other']
FI = {s: i for i, s in enumerate(FINE)}
GROUPS = {  # reported strata -> fine strata
    'SNV': ['SNV'],
    'INDEL': [f for f in FINE if f[:3] in ('INS', 'DEL')],
    'INS': [f for f in FINE if f.startswith('INS')],
    'DEL': [f for f in FINE if f.startswith('DEL')],
    'INDEL_HP>=7': ['INS_HP>=7', 'DEL_HP>=7'],
    'INDEL_not_HP>=7': ['INS_HP2-6', 'INS_STR', 'INS_other', 'DEL_HP2-6', 'DEL_STR', 'DEL_other'],
    'INDEL_STR': ['INS_STR', 'DEL_STR'],
    'INDEL_HP2-6_other': ['INS_HP2-6', 'INS_other', 'DEL_HP2-6', 'DEL_other'],
    'SV>50': ['SV>50'],
}

def fine_code(nkind, rep_class, ilen):
    if nkind == 'SNV':
        return FI['SNV']
    if nkind in ('INS', 'DEL'):
        if abs(int(ilen)) > 50:
            return FI['SV>50']
        return FI[f'{nkind}_{rep_class}']
    return FI['other']

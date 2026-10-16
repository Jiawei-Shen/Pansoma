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
import hashlib, os

# GRCh38 no-alt analysis set, the fasta of every labels manifest. Slurm jobs copy it to node-local /tmp and set E7_FASTA,
# because random access over BeeGFS stalled on busy nodes (2026-10-13).
FASTA_SRC = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
FASTA = os.environ.get('E7_FASTA', FASTA_SRC)
HPRC_VCF = '/scratch/jshen/data/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz'
BCFTOOLS = '/wanglab/jshen/anaconda3/envs/sambcfenv/bin/bcftools'
WORK = '/scratch/jshen/data/reviewer_experiments_20261013/E7_graph_composition'
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
LONG = 200


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

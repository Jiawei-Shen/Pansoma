"""E7 step 2: normalize a sample's somatic truth and germline (normal) alleles exactly like the HPRC alleles.

Input: the truth tables the labels were built from (tensors/*/labels.manifest.json -> truth.{somatic,germline}.table),
one row per (record, ALT allele) with truth_id, vcf_pos / vcf_ref / vcf_alt, kind, passed (FILTER PASS), in_bed (inside
that truth's BED). truth_id is the same in every tensor set built from the same truth VCF (checked for the somatic tables
of the d9 / linear / d46 sets: identical truth_id, chrom, pos, ref, alt, filters, passed, in_bed), so the output joins
to the per-allele recall tables (tensors/{somatic,germline}.recall.tsv).
Each allele -> VCF line (ID = truth_id, INFO carries kind / passed / in_bed) -> `bcftools norm -f GRCh38 -c w`
(left-align) -> common.primitives (MNP -> SNVs; a truth allele with several primitives keeps all their keys).
Output: WORK/sample_alleles/<sample>.<role>.tsv.gz: truth_id, chrom, kind, passed, in_bed, keys (';'), nkind,
        rep_class, unit, tract, ilen.
Usage: s2_sample_alleles.py SAMPLE TENSOR_SET [ROLES, default somatic,germline]
"""
import gzip, os, subprocess, sys, threading
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import AUTOSOMES, BCFTOOLS, FASTA, WORK, key, primitives, repeat_context

sample, tset = sys.argv[1], sys.argv[2]
SETS = '/scratch/jshen/data/pansoma_v2_tensors'
os.makedirs(f'{WORK}/sample_alleles', exist_ok=True)
fai = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(FASTA + '.fai')}
fa = pysam.FastaFile(FASTA)
order = {c: i for i, c in enumerate(AUTOSOMES)}
for role in (sys.argv[3].split(',') if len(sys.argv) > 3 else ('somatic', 'germline')):
    table = f'{SETS}/{tset}/truth/{role}.graph.tsv'
    rows = []
    with open(table) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {c: hdr.index(c) for c in ('truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind', 'passed', 'in_bed')}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            c = f[ix['chrom']]
            if c not in order or f[ix['kind']] == 'COMPLEX' and not f[ix['vcf_alt']]:
                continue
            rows.append((order[c], int(f[ix['vcf_pos']]), f[ix['truth_id']], c, f[ix['vcf_ref']], f[ix['vcf_alt']],
                         f[ix['kind']], f[ix['passed']][0], f[ix['in_bed']][0]))
    rows.sort()
    out = f'{WORK}/sample_alleles/{sample}.{role}.tsv.gz'
    norm = subprocess.Popen(f'{BCFTOOLS} norm -f {FASTA} -c w -Ov - 2> {out}.norm.log | {BCFTOOLS} query -f '
                            '"%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/K\\t%INFO/P\\t%INFO/B\\n"', shell=True,
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    res = {}

    def reader():
        for line in norm.stdout:
            tid, c, pos, ref, alt, k, p, b = line.rstrip('\n').split('\t')
            res.setdefault(tid, [c, k, p, b, []])[4].extend(primitives(int(pos), ref, alt))

    t = threading.Thread(target=reader)
    t.start()
    norm.stdin.write('##fileformat=VCFv4.2\n' + ''.join(f'##contig=<ID={c},length={fai[c]}>\n' for c in AUTOSOMES) +
                     '##INFO=<ID=K,Number=1,Type=String,Description="truth table kind">\n'
                     '##INFO=<ID=P,Number=1,Type=String,Description="passed">\n'
                     '##INFO=<ID=B,Number=1,Type=String,Description="in_bed">\n'
                     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for _, pos, tid, c, ref, alt, k, p, b in rows:
        norm.stdin.write(f'{c}\t{pos}\t{tid}\t{ref}\t{alt}\t.\t.\tK={k};P={p};B={b}\n')
    norm.stdin.close()
    t.join()
    assert norm.wait() == 0
    n_in = len(rows)
    del rows
    with gzip.open(out + '.tmp', 'wt', compresslevel=3) as fo:
        fo.write('truth_id\tchrom\tkind\tpassed\tin_bed\tkeys\tnkind\trep_class\tunit\ttract\tilen\n')
        for tid, (c, k, p, b, prims) in res.items():
            keys = ';'.join(key(c, *pr[:3]) for pr in prims)
            if len(prims) == 1:
                pos, ref, alt, nk = prims[0]
                rc, u, tr = repeat_context(fa, c, pos, ref, alt, nk)
                il = len(alt) - len(ref)
            else:
                nk, rc, u, tr, il = ('MULTI' if prims else 'NONE'), '', 0, 0, 0
            fo.write(f'{tid}\t{c}\t{k}\t{p}\t{b}\t{keys}\t{nk}\t{rc}\t{u}\t{tr}\t{il}\n')
    os.replace(out + '.tmp', out)
    print(sample, role, 'table alleles', n_in, 'normalized', len(res), flush=True)

"""Step 5: which HPRC v1.1 haplotypes carry each truth allele, from the vg deconstruct VCFs (by sequence, per locus).

Inputs: loci.tsv (window win_start0..win_end0, ref_hap, alt_hap); hprc_sample_metadata.tsv (s5b);
  d9    /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz (deconstruct of the graph Pansoma aligns to;
        top-level snarls only, no LV tag; AT = allele traversal in d9 node ids, comparable with loci.tsv nodes)
  full  /scratch/jshen/data/hprc-v1.1-mc-grch38.raw.vcf.gz (deconstruct of the unfiltered default graph, nested: LV/PS;
        AT in default-graph node ids, NOT comparable with d9 ids); tabix index built by s5_index_full_vcf.sh into
        $D/full_vcf_index/ (symlink + .tbi; the original has none)
Method (the i1_germline.py rule, extended to records crossing the window edge): every record overlapping
[win_start0, win_end0]; an ALT allele 'carries the truth ALT' if, applied alone to GRCh38, it gives alt_hap (when the record
reaches past the window, both sides are extended with the same GRCh38 flank). Primary record = lowest LV (d9: all 0),
then largest INFO AC. Its phased GT: haplotype k of the GT = hap k of the sample (vg deconstruct order; checked at start
on chrX 50-52 Mb: male samples, whose hap 2 = maternal carries chrX, must be '.|x'); CHM13 is haploid and written CHM13#0
(its GFA W-line hap); GRCh38 = allele 0 of every record, so never a carrier. carriers_any = union over every matching
record: in the full VCF a nested child snarl can carry the allele inside a parent allele with other variants, and in both
VCFs a repeat-unit INDEL can be matched by two records at shifted positions (each holds part of the carriers), so
carriers_any is the carrier set to use. info_check: INFO AC/AN vs the GT count of the primary record.
Limitation: a haplotype with the truth allele plus another variant inside the same record span has a different allele and
is not a carrier here (the GFA-walk step handles that); a locus with no matching record may still be spelled by a path.
Outputs: $D/hprc_vcf_alleles.tsv (one row per truth_id x source d9|full), $D/hprc_vcf_haplotypes.tsv.gz (truth_id, source,
sample, hap, allele_index, carries_alt, carries_alt_any; loci with a matching record only).
"""
import collections, csv, gzip, statistics, sys
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
VCF = {'d9': pysam.VariantFile('/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz'),
       'full': pysam.VariantFile(f'{D}/full_vcf_index/hprc-v1.1-mc-grch38.raw.vcf.gz')}
fa = pysam.FastaFile('/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta')
META = list(csv.DictReader(open(f'{D}/hprc_sample_metadata.tsv'), delimiter='\t'))
SUP, SEX = {r['sample']: r['superpopulation'] for r in META}, {r['sample']: r['sex'] for r in META}
SAMPLES = list(VCF['d9'].header.samples)
assert SAMPLES == list(VCF['full'].header.samples)
MALE = [s for s in SAMPLES if SEX.get(s) == 'male']
phase = collections.Counter(('hap1 missing' if g[0] is None else 'hap2 missing')
                            for r in VCF['d9'].fetch('chrX', 50000000, 52000000) for s in MALE
                            for g in [r.samples[s]['GT']] if len(g) == 2 and (g[0] is None) != (g[1] is None))
print('phase check, male samples on chrX 50-52 Mb (GT order = hap order expects hap1 missing):', dict(phase))
loci = list(csv.DictReader(open(f'{D}/loci.tsv'), delimiter='\t'))[:int(sys.argv[1]) if sys.argv[1:] else None]  # argv: test on the first N
short = lambda s: s if len(s) <= 30 else f'{s[:12]}..{s[-6:]}({len(s)}bp)'


def haps(rec):
    """[(sample, hap, allele index or None)] in GT order."""
    out = []
    for s in SAMPLES:
        gt = rec.samples[s]['GT'] or (None,)
        out += [(s, '0' if s == 'CHM13' else str(k + 1), g) for k, g in enumerate(gt)]
    return out


def matching(vf, chrom, w0, w1, alt_hap):
    """records overlapping the window: [(rec, {matching ALT indices})], n overlapping, n REF != GRCh38."""
    hits, n, bad = [], 0, 0
    for rec in vf.fetch(chrom, w0, w1 + 1):
        n += 1
        lo, hi = min(w0, rec.start), max(w1 + 1, rec.stop)
        region = fa.fetch(chrom, lo, hi).upper()
        if region[rec.start - lo:rec.stop - lo] != rec.ref.upper():
            bad += 1; continue
        target = region[:w0 - lo] + alt_hap + region[w1 + 1 - lo:]
        m = {i + 1 for i, a in enumerate(rec.alts or ()) if a and a[0] not in '<*'
             and region[:rec.start - lo] + a.upper() + region[rec.stop - lo:] == target}
        if m:
            hits.append((rec, m))
    return hits, n, bad


cols = ['truth_id', 'source', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'n_records_overlapping', 'n_records_matching',
        'record', 'lv', 'record_start0', 'record_end0', 'allele_index', 'at', 'AC', 'AN', 'AF', 'carriers', 'n_carriers',
        'n_other', 'n_missing', 'info_check', 'carrier_superpops', 'conflict', 'other_records', 'carriers_any', 'n_carriers_any']
fo, fh = open(f'{D}/hprc_vcf_alleles.tsv', 'w'), gzip.open(f'{D}/hprc_vcf_haplotypes.tsv.gz', 'wt')
fo.write('\t'.join(cols) + '\n')
fh.write('truth_id\tsource\tsample\thap\tallele_index\tcarries_alt\tcarries_alt_any\n')
stat = collections.Counter()
for i, L in enumerate(loci):
    chrom, w0, w1 = L['chrom'], int(L['win_start0']), int(L['win_end0'])
    for src, vf in VCF.items():
        hits, n, bad = matching(vf, chrom, w0, w1, L['alt_hap'])
        stat[(src, 'ref_mismatch')] += bad
        out = {c: L.get(c, '') for c in cols}
        out.update(source=src, n_records_overlapping=n, n_records_matching=len(hits))
        if hits:
            lv = lambda r: int(r.info.get('LV', 0)) if 'LV' in r.header.info else 0
            hits.sort(key=lambda h: (lv(h[0]), -sum(h[0].info['AC'][k - 1] for k in h[1]), h[0].pos))
            rec, m = hits[0]
            H = haps(rec)
            anym = {(s, h) for r, mm in hits for s, h, g in haps(r) if g in mm}
            car = [f'{s}#{h}' for s, h, g in H if g in m]
            ac, an = sum(rec.info['AC'][k - 1] for k in m), rec.info.get('AN')
            n_oth, n_mis = sum(g is not None and g not in m for _, _, g in H), sum(g is None for _, _, g in H)
            at = rec.info.get('AT')
            out.update(record=f'{chrom}:{rec.pos}:{short(rec.ref)}:{",".join(short(rec.alts[k - 1]) for k in sorted(m))}',
                       lv=lv(rec) if src == 'full' else '', record_start0=rec.start, record_end0=rec.stop - 1,
                       allele_index=','.join(map(str, sorted(m))), at=','.join(at[k] for k in sorted(m)) if at else '',
                       AC=ac, AN=an, AF=f"{sum(rec.info['AF'][k - 1] for k in m):.4f}", carriers=','.join(car),
                       n_carriers=len(car), n_other=n_oth, n_missing=n_mis,
                       info_check='agree' if (ac, an) == (len(car), len(car) + n_oth) else f'DIFF calc {len(car)}/{len(car) + n_oth}',
                       carrier_superpops=','.join(f'{k}:{v}' for k, v in sorted(collections.Counter(
                           'CHM13' if s == 'CHM13' else SUP.get(s, '?') for s, h, g in H if g in m).items())),
                       conflict=','.join(rec.info.get('CONFLICT', ())),
                       other_records=';'.join(f'{r.pos}:LV{lv(r)}:AC{sum(r.info["AC"][k - 1] for k in mm)}' for r, mm in hits[1:]),
                       carriers_any=','.join(f'{s}#{h}' for s, h, _ in H if (s, h) in anym), n_carriers_any=len(anym))
            stat[(src, 'info_' + out['info_check'][:4])] += 1
            for s, h, g in H:
                fh.write(f'{L["truth_id"]}\t{src}\t{s}\t{h}\t{"" if g is None else g}\t{g in m}\t{(s, h) in anym}\n')
        fo.write('\t'.join(str(out[c]) for c in cols) + '\n')
    if i % 250 == 0:
        print(i, len(loci), flush=True)
fo.close(); fh.close()
rows = list(csv.DictReader(open(f'{D}/hprc_vcf_alleles.tsv'), delimiter='\t'))
for src in VCF:
    for kd in ('SNV', 'INDEL', 'all'):
        s = [r for r in rows if r['source'] == src and (kd == 'all' or r['kind2'] == kd)]
        h = [r for r in s if int(r['n_records_matching'])]
        af = [float(r['AF']) for r in h]
        print(src, kd, len(s), 'loci;', len(h), 'with a matching record;', 'median AF', statistics.median(af) if af else '-',
              '; median carriers', statistics.median(int(r['n_carriers']) for r in h) if h else '-',
              '; >1 matching record', sum(int(r['n_records_matching']) > 1 for r in h))
print(dict(stat))

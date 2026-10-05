"""COLO829T step 6: which HPRC v1.1 haplotypes carry each truth allele, from the vg deconstruct VCFs (by sequence, per locus).

Port of analysis/graph_absorbed_somatic_20261001/s5_hprc_vcf.py (HG008T); same rule and output columns (+ carrier_superpops_any).
Inputs: $D/loci.tsv (c1b; window win_start0..win_end0, ref_hap, alt_hap); $D/hprc_sample_metadata.tsv (HPRC v1.1 sample ->
  1000G population / superpopulation / sex; copied by this script from the HG008 data dir, built there by s5b, not re-downloaded);
  d9    /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz (deconstruct of the graph Pansoma aligns to;
        top-level snarls only, no LV tag; AT = allele traversal in d9 node ids, comparable with loci.tsv / graph_paths.tsv nodes)
  full  /scratch/jshen/data/hprc-v1.1-mc-grch38.raw.vcf.gz (deconstruct of the unfiltered default graph, nested: LV/PS;
        AT in default-graph node ids, NOT comparable with d9 ids), read through the HG008 symlink + tabix index
        (graph_absorbed_somatic_20261001/full_vcf_index, s5_index_full_vcf.sh; read only, the original has no index)
Method (as s5): every record overlapping [win_start0, win_end0]; an ALT allele 'carries the truth ALT' if, applied alone to
GRCh38, it gives alt_hap (when the record reaches past the window, both sides are extended with the same GRCh38 flank).
Primary record = lowest LV (d9: all 0), then largest INFO AC, then position. Its phased GT: haplotype k of the GT = hap k of
the sample (vg deconstruct order; checked at start on chrX 50-52 Mb: male samples, whose hap 2 = maternal carries chrX, must be
'.|x'); CHM13 is haploid and written CHM13#0 (its GFA W-line hap); GRCh38 = allele 0 of every record, so never a carrier.
carriers_any = union over every matching record (a nested child snarl in the full VCF, or a repeat-unit INDEL matched by two
records at shifted positions, can each hold part of the carriers), so carriers_any is the carrier set to use.
info_check: INFO AC/AN vs the GT count of the primary record.
Limitation: a haplotype with the truth allele plus another variant inside the same record span has a different allele and is
not a carrier here (the GFA-walk step handles that); a locus with no matching record may still be spelled by a graph path.
Outputs: $D/hprc_vcf_alleles.tsv (one row per truth_id x source d9|full, loci.tsv order), $D/hprc_vcf_haplotypes.tsv.gz
(truth_id, source, sample, hap, allele_index, carries_alt, carries_alt_any; loci with a matching record only).
Run: python c6_hprc_vcf.py [N]  (login node; N = test on the first N loci); python c6_hprc_vcf.py report  re-prints the
summary and the checks (vs graph_paths.tsv: d9 matching record per c2 match class; d9 AT inside the c2 primary path; d9 vs
full carriers_any).
Assumptions: COLO829T has no chrX/chrY truth, so the male hap-order caveat does not touch any locus; AF is INFO AF of the
matching ALT(s) of the primary record (called haplotypes incl. CHM13; AN varies with missing GTs).
"""
import collections, csv, filecmp, gzip, os, re, shutil, statistics, sys
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005'
H = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
VCFP = {'d9': '/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz',
        'full': f'{H}/full_vcf_index/hprc-v1.1-mc-grch38.raw.vcf.gz'}
rd = lambda p: list(csv.DictReader(gzip.open(p, 'rt') if p.endswith('.gz') else open(p), delimiter='\t'))
if not os.path.exists(f'{D}/hprc_sample_metadata.tsv'):
    shutil.copyfile(f'{H}/hprc_sample_metadata.tsv', f'{D}/hprc_sample_metadata.tsv')
assert filecmp.cmp(f'{H}/hprc_sample_metadata.tsv', f'{D}/hprc_sample_metadata.tsv', shallow=False)
META = rd(f'{D}/hprc_sample_metadata.tsv')
SUP, SEX = {r['sample']: r['superpopulation'] for r in META}, {r['sample']: r['sex'] for r in META}
cols = ['truth_id', 'source', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'n_records_overlapping', 'n_records_matching',
        'record', 'lv', 'record_start0', 'record_end0', 'allele_index', 'at', 'AC', 'AN', 'AF', 'carriers', 'n_carriers',
        'n_other', 'n_missing', 'info_check', 'carrier_superpops', 'conflict', 'other_records', 'carriers_any', 'n_carriers_any',
        'carrier_superpops_any']
spc = lambda ks: ','.join(f'{k}:{v}' for k, v in sorted(collections.Counter(
    'CHM13' if s == 'CHM13' else SUP.get(s, '?') for s in (x.split('#')[0] for x in ks)).items()))
walk = lambda a, rev=False: (lambda w: [('<' if o == '>' else '>', n) for o, n in w[::-1]] if rev else w)(re.findall(r'([<>])(\d+)', a))
sub = lambda a, b: any(b[i:i + len(a)] == a for i in range(len(b) - len(a) + 1))
short = lambda s: s if len(s) <= 30 else f'{s[:12]}..{s[-6:]}({len(s)}bp)'


def run(n=None):
    VCF = {k: pysam.VariantFile(p) for k, p in VCFP.items()}
    fa = pysam.FastaFile('/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta')
    SAMPLES = list(VCF['d9'].header.samples)
    assert SAMPLES == list(VCF['full'].header.samples) and set(SAMPLES) <= set(SUP), 'sample sets differ'
    MALE = [s for s in SAMPLES if SEX.get(s) == 'male']
    phase = collections.Counter(('hap1 missing' if g[0] is None else 'hap2 missing')
                                for r in VCF['d9'].fetch('chrX', 50000000, 52000000) for s in MALE
                                for g in [r.samples[s]['GT']] if len(g) == 2 and (g[0] is None) != (g[1] is None))
    print('phase check, male samples on chrX 50-52 Mb (GT order = hap order expects hap1 missing):', dict(phase))
    assert phase and set(phase) == {'hap1 missing'}
    loci = rd(f'{D}/loci.tsv')[:n]

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

    fo, fh = open(f'{D}/hprc_vcf_alleles.tsv', 'w'), gzip.open(f'{D}/hprc_vcf_haplotypes.tsv.gz', 'wt')
    fo.write('\t'.join(cols) + '\n')
    fh.write('truth_id\tsource\tsample\thap\tallele_index\tcarries_alt\tcarries_alt_any\n')
    stat = collections.Counter()
    for i, L in enumerate(loci):
        assert L['status'] == 'ok', L['truth_id']
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
                Hp = haps(rec)
                anym = {(s, h) for r, mm in hits for s, h, g in haps(r) if g in mm}
                car = [f'{s}#{h}' for s, h, g in Hp if g in m]
                cany = [f'{s}#{h}' for s, h, _ in Hp if (s, h) in anym]
                ac, an = sum(rec.info['AC'][k - 1] for k in m), rec.info.get('AN')
                n_oth, n_mis = sum(g is not None and g not in m for _, _, g in Hp), sum(g is None for _, _, g in Hp)
                at = rec.info.get('AT')
                out.update(record=f'{chrom}:{rec.pos}:{short(rec.ref)}:{",".join(short(rec.alts[k - 1]) for k in sorted(m))}',
                           lv=lv(rec) if src == 'full' else '', record_start0=rec.start, record_end0=rec.stop - 1,
                           allele_index=','.join(map(str, sorted(m))), at=','.join(at[k] for k in sorted(m)) if at else '',
                           AC=ac, AN=an, AF=f"{sum(rec.info['AF'][k - 1] for k in m):.4f}", carriers=','.join(car),
                           n_carriers=len(car), n_other=n_oth, n_missing=n_mis,
                           info_check='agree' if (ac, an) == (len(car), len(car) + n_oth) else f'DIFF calc {len(car)}/{len(car) + n_oth}',
                           carrier_superpops=spc(car), conflict=','.join(rec.info.get('CONFLICT', ())),
                           other_records=';'.join(f'{r.pos}:LV{lv(r)}:AC{sum(r.info["AC"][k - 1] for k in mm)}' for r, mm in hits[1:]),
                           carriers_any=','.join(cany), n_carriers_any=len(cany), carrier_superpops_any=spc(cany))
                stat[(src, 'info_' + out['info_check'][:4])] += 1
                for s, h, g in Hp:
                    fh.write(f'{L["truth_id"]}\t{src}\t{s}\t{h}\t{"" if g is None else g}\t{g in m}\t{(s, h) in anym}\n')
            fo.write('\t'.join(str(out[c]) for c in cols) + '\n')
        if i % 100 == 0:
            print(i, len(loci), flush=True)
    fo.close(); fh.close()
    print(dict(stat))


def report():
    rows, G = rd(f'{D}/hprc_vcf_alleles.tsv'), {r['truth_id']: r for r in rd(f'{D}/graph_paths.tsv')}
    for src in ('d9', 'full'):
        for kd in ('SNV', 'INDEL', 'all'):
            s = [r for r in rows if r['source'] == src and (kd == 'all' or r['kind2'] == kd)]
            h = [r for r in s if int(r['n_records_matching'])]
            af = [float(r['AF']) for r in h]
            print(src, kd, len(s), 'loci;', len(h), 'with a matching record;', 'median AF', statistics.median(af) if af else '-',
                  '; median carriers', statistics.median(int(r['n_carriers']) for r in h) if h else '-',
                  '; median carriers_any', statistics.median(int(r['n_carriers_any']) for r in h) if h else '-',
                  '; >1 matching record', sum(int(r['n_records_matching']) > 1 for r in h),
                  '; info_check', dict(collections.Counter(r['info_check'][:4] for r in h)))
    R = {(r['truth_id'], r['source']): r for r in rows}
    ids = [r['truth_id'] for r in rows if r['source'] == 'd9']
    print('d9 matching record x c2 match class (kind2, match, d9 record, full record):')
    for k, v in sorted(collections.Counter((R[t, 'd9']['kind2'], G[t]['match'], R[t, 'd9']['n_records_matching'] != '0',
                                            R[t, 'full']['n_records_matching'] != '0') for t in ids).items()):
        print('  ', k, v)
    at = collections.Counter()
    for t in ids:
        r = R[t, 'd9']
        if G[t]['match'] == 'exact' and r['at']:
            pp = walk(G[t]['primary_path'])
            at[any(sub(walk(a), pp) or sub(walk(a, True), pp) for a in r['at'].split(','))] += 1
    print('exact loci with a d9 matching record: d9 AT (matching allele walk) inside the c2 primary path:', dict(at))
    cmp = collections.Counter()
    for t in ids:
        a, b = (set(R[t, s]['carriers_any'].split(',')) - {''} for s in ('d9', 'full'))
        if a or b:
            cmp[('identical' if a == b else 'd9 < full' if a < b else 'full < d9' if b < a else 'differ', R[t, 'd9']['kind2'])] += 1
    print('carriers_any d9 vs full (loci with a match in either):', sorted(cmp.items()))
    sp = collections.Counter()
    for t in ids:
        for x in R[t, 'd9']['carriers_any'].split(','):
            if x:
                sp['CHM13' if x.startswith('CHM13') else SUP.get(x.split('#')[0], '?')] += 1
    print('d9 carriers_any haplotype-loci per superpopulation:', dict(sp))


if sys.argv[1:] == ['report']:
    report()
else:
    run(int(sys.argv[1]) if sys.argv[1:] else None)
    report()

# Result (2026-10-05, login node, 2.8 min, 2.9 GB; log tmp/graph_absorbed_somatic_colo829t_20261005/c6_hprc_vcf.login.log):
#   phase check 48,921 'hap1 missing' (GT order = hap order); REF vs GRCh38 mismatches 0; info_check agree on every primary record
#   matching record: d9 395 / 497 (SNV 62/67, INDEL 333/430), full 410 / 497 (SNV 65, INDEL 345); median AF d9 0.205, full 0.180;
#     median carriers_any 17 (d9) / 16.5 (full) of 89 haplotypes (44 samples x 2 + CHM13#0); >1 matching record d9 6, full 19
#   vs c2: every d9 match is a c2 'exact' locus (none for with_germline / closest); d9 AT walk = c2 primary path in 392 / 395,
#     the other 3 (3232, 21304, 40338) are another c2 exact path (same repeat INDEL at a shifted copy, element in graph_elements)
#   exact loci without a d9 match: 63 (61 INDEL, 2 SNV; 13 of them matched in the full VCF; 35165/35166 chr14 ~19.2 Mb have no
#     d9 record at all); closest / with_germline loci matched only in the full VCF: 4 (allele in the unfiltered graph only)
#   carriers_any d9 vs full (412 loci with a match in either): identical 344, d9 subset of full 59, full subset of d9 7, differ 2;
#     d9-only carriers = haplotypes with another variant inside the full record span (different allele there) - representation
#     difference, resolved at allele level by the GFA-walk step; d9 carriers_any include CHM13#0 at 89 loci
